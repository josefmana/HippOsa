# Computes regressions evaluating the effect of OSA on brain structures' volume:
#
# RQ1) the moderating effect OSA on subcortical brain structures' volume in de novo PD
# RQ2) the moderating effect OSA on hippocampal substructures' volume in de novo PD


rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

# load packages
library(here)
library(tidyverse)
library(marginaleffects)
library(performance)
library(ggh4x)

theme_set( theme_bw() )

# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# read helper files with variable names
subco <- read.csv( here("helpers","subcortical.csv"), sep = ",") # subcortical structures
hippo <- read.csv( here("helpers","hippocampus.csv"), sep = ",") %>% filter( complete.cases(name) ) # hippocampal structures

# in-house functions
source( here("scripts","utils.R") )


# PRE-PROCESSING  ----

# read the data
d0 <-

  read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>% # read data
  filter(event == "enrollment") %>% # keep enrollment only
  
  # re-calculate hippocampal fields according to the legend in hippo
  mutate(
    !!!setNames( rep(NA, length( unique(hippo$name) ) ), unique(hippo$name) ),
    across(
      .cols = unique(hippo$name),
      .fns = ~ rowSums( across( all_of( with( hippo, var[name == cur_column()] ) ) ) )
    )
  )

# format it for analyses
df <- d0 %>%
  
  # pre-process brain and demography variables
  mutate_if( is.character, as.factor ) %>%
  mutate(
    PD = if_else(SUBJ == "PD", 1, 0),
    GENDER = as.factor(GENDER),
    across( c("EDU.Y","AGE","BMI"), ~ as.numeric( scale(.x) ) )
  ) %>%
  mutate(
    SBTIV = as.numeric( scale(SBTIV) ),
    across( all_of(subco$name), ~ as.numeric( scale(.x) ) ),
    across( all_of( unique(hippo$name) ), ~ as.numeric( scale(.x) ) )
  ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
    contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
    contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  } )


# INTERACTION BOXPLOTS ----

d0 %>%
  
  select(SUBJ, AHI.F, all_of(subco$scaled) ) %>%
  pivot_longer(
    cols = all_of(subco$scaled),
    values_to = "volume",
    names_to = "struct"
  ) %>%
  mutate(
    side = unlist(
      sapply( 1:nrow(.), function(i) with( subco, side[scaled == struct[i]] ) ),
      use.names = F
    ),
    structure = factor(
      unlist(
        sapply( 1:nrow(.), function(i) with( subco, structure[scaled == struct[i]] ) ),
        use.names = F
      ),
      levels = unique(subco$structure),
      ordered = T
    ),
    Diagnosis = if_else(SUBJ == "PD", "PD", "HC"),
    `OSA: ` = factor(
      if_else(AHI.F == "H", "OSA+", "OSA-"),
      levels = c("OSA-","OSA+"),
      ordered = T
    )
  ) %>%
  
  ggplot() +
  aes(y = volume, x = Diagnosis) +
  geom_boxplot( aes(fill = `OSA: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
  geom_dotplot( aes(fill = `OSA: `), binaxis = "y", stackdir = "center", position = position_dodge(.7), dotsize = 1.5 ) +
  labs( y = bquote("Standardized volume "(mm^3) ) ) +
  facet_grid2( structure ~ side, scales = "free_y", independent = "y" ) +
  scale_fill_manual( values = c("#64CDCC","#F9A729") ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","subcort_boxplots.jpg"),
  dpi = 300,
  width = 6,
  height = 10
)


# LINEAR REGRESSIONS ----

# set-up formulas
preds <- data.frame(
  
  model = c("base", "interaction","bmi"),
  X = c("SUBJ * AHI.F + AGE + GENDER + SBTIV", "SUBJ * AHI.F * AGE + GENDER + SBTIV","SUBJ * AHI.F * (AGE + BMI) + GENDER + SBTIV")
  
)

# loop through types of regressions (base vs interaction), and structures
lapply(
  
  1:nrow(preds),
  function(r)
    
    lapply(
      
      c("subco","hippo"),
      function(i) {
        
        with(
          
          get(i), {
            
            
            # fit the models
            fit <- fit_reg(df, unique(name), X = preds[r, "X"], w = F)
            
            # regression coefficients with threshold-based decisions
            # do full analysis for subcortical structures and interaction only for hippocampal substructures
            if (i == "subco") write.table(
              
              # extract and write main effects and interactions for subcortical structures
              x = left_join(
                
                rbind.data.frame(
                  lm_coeff(fit, term = "SUBJ1"),
                  lm_coeff(fit, term = "AHI.F1"),
                  lm_coeff(fit, term = "SUBJ1:AHI.F1")
                ),
                lm_dia(fit),
                by = c("y","X")
                
              ) %>% mutate( sig_FDR = bh_adjust(`p value`) ), # re-calculate Benjamini-Hochberg adjusted significance statements
              
              file = here( "tables", paste0(i, "_", preds[r, "model"], "_regression_coefficients.csv") ),
              sep = ",", row.names = F, quote = F
              
            ) else if (i == "hippo") write.table(
              
              # extract and write 'main effects and interactions only for hippocampal substructures
              x = left_join( lm_coeff(fit, term = "SUBJ1:AHI.F1") ,lm_dia(fit), by = c("y","X") ),
              file = here( "tables", paste0(i, "_" , preds[r, "model"], "_regression_coefficients.csv") ),
              sep = ",", row.names = F, quote = F
            )
            
            
            # marginal effects
            write.table(
              
              # marginal effects
              x = left_join(
                meff( fit = fit, fit0 = fit, type = ifelse(i == "subco", "full", "moderation") ), # marginalised effects
                lm_dia(fit), # diagnostics
                by = c("y","X")
              ),
              
              file = here( "tables", paste0(i, "_", preds[r, "model"], "_marginal_effects.csv") ),
              sep = ",", row.names = F, quote = F
            )
            
          }
        )
        
      }
    )
  
)


# FOREST PLOTS ----

## ---- subcortical structures regressions ----

left_join(
  
  # the marginal effects
  read.csv(here("tables","subco_base_marginal_effects.csv"), sep = ","),
  
  # add BH adjusted significance statements
  read.csv(here("tables","subco_base_marginal_effects.csv"), sep = ",") %>%
    filter( is.na(SUBJ) ) %>%
    mutate( sig = bh_adjust(p.value) ) %>%
    filter(term == "PD - CON") %>%
    select(y, sig)
  
) %>%
  
  # prepare variables
  mutate(
    
    # group the estimated is conditioned on
    `Group: ` = case_when(
      contrast == "mean(L) - mean(H)" & SUBJ == "CON" ~ "HC",
      contrast == "mean(L) - mean(H)" & SUBJ == "PD" ~ "PD",
      term == "PD - CON" ~ "PD - HC"
    ),
    
    # significance statement
    `Sig. (5% FDR):` = if_else(sig == "*", T, F),
    
    # hemisphere of the outcome variable
    side = unlist(
      sapply( 1:nrow(.), function(i) with( subco, side[name == y[i]] ) ),
      use.names = F
    ),
    
    # outcome variable brain structure
    structure = factor(
      unlist(
        sapply( 1:nrow(.), function(i) with( subco, structure[name == y[i]] ) ),
        use.names = F
      ),
      levels = rev( unique(subco$structure) ),
      ordered = T
    )

  ) %>%
  
  # keep only rows and columns of interest/use
  filter( complete.cases(`Group: `) ) %>%
  select(y, side, structure, `Group: `, `Sig. (5% FDR):`, estimate, conf.low, conf.high) %>%
  
  # plotting proper
  ggplot() +
  aes(y = estimate, ymin = conf.low, ymax = conf.high, x = structure, shape = `Group: `, colour = `Sig. (5% FDR):`) +
  geom_point(position = position_dodge(width = .3), size = 3.3) + # point estimates
  geom_linerange(position = position_dodge(width = .3), linewidth = 1) + # 95% CIs
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") + # zero for reference
  coord_flip() + # use flip if decided to cut the range of values shown 
  facet_wrap( ~ side, ncol = 2 ) + # a column per hemisphere
  scale_colour_manual( values = c("grey","orange") ) + # simple effects blue, interaction red
  labs(x = NULL, y = "mean(OSA-) - mean(OSA+)") + # estimate is that of difference between OSA- and OSA+ expected means
  theme(legend.position = "right")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","subcort_forest.jpg"),
  dpi = 300,
  width = 8,
  height = 8
)


## ---- hippocampal structures regressions ----

# prepare order of structures
ord <- hippo[ c("structure","order") ] %>% unique() %>% arrange(order) %>% select(structure) %>% unlist(use.names = F)

# plot it
left_join(
  
  # marginal effects
  read.csv(here("tables","hippo_base_marginal_effects.csv"), sep = ","),
  
  # add BH adjusted significance statements for interactions
  read.csv(here("tables","hippo_base_marginal_effects.csv"), sep = ",") %>%
    filter(term == "PD - CON") %>%
    mutate( sig = bh_adjust(p.value) ) %>%
    select(y, sig)
  
) %>%
  
  # prepare variables
  mutate(
    
    # group the estimated is conditioned on
    `Group: ` = case_when(
      contrast == "mean(L) - mean(H)" & SUBJ == "CON" ~ "HC",
      contrast == "mean(L) - mean(H)" & SUBJ == "PD" ~ "PD",
      term == "PD - CON" ~ "PD - HC"
    ),
    
    # significance statement
    `Sig. (5% FDR):` = if_else(sig == "*", T, F),
    
    # hemisphere of the outcome variable
    side = case_when(
      grepl("_left", y) ~ "Left",
      grepl("_right", y) ~ "Right"
    ),
    
    # outcome variable brain structure
    struct = factor(
      unlist( sapply( 1:nrow(.), function(i) with( hippo, unique( structure[name == y[[i]]] ) ) ), use.names = F),
      levels = rev(ord),
      ordered = T
    ),
    
    # block
    block = factor(
      unlist( sapply( 1:nrow(.), function(i) with( hippo, unique( block[name == y[[i]]] ) ) ), use.names = F),
      levels = unique(hippo$block),
      ordered = T
    )
    
  ) %>%
  
  # keep only rows and columns of interest/use
  filter( complete.cases(`Group: `) ) %>%
  select(y, side, struct, block, `Group: `, `Sig. (5% FDR):`, estimate, conf.low, conf.high) %>%
  
  # plotting proper
  ggplot() +
  aes(y = estimate, ymin = conf.low, ymax = conf.high, x = struct, shape = `Group: `, colour = `Sig. (5% FDR):`) +
  geom_point(position = position_dodge(width = .5), size = 3) + # point estimates
  geom_linerange(position = position_dodge(width = .5), linewidth = 1) + # 95% CIs
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") + # zero for reference
  coord_flip() + # use flip if decided to cut the range of values shown 
  facet_grid(block ~ side, scale = "free_y", space = "free") + # a column per hemisphere
  scale_colour_manual( values = c("grey","orange") ) + # simple effects blue, interaction red
  labs(x = NULL, y = "mean(OSA-) - mean(OSA+)") + # estimate is that of difference between OSA- and OSA+ expected means
  theme(legend.position = "right")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","hippo_forest.jpg"),
  dpi = 300,
  width = 8,
  height = 12.5
)

