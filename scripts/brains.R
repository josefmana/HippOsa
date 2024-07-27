# Computes regressions evaluating the effect of OSA on brain structures' volume:
#
# RQ1) the moderating effect OSA on subcortical brain structures' volume in de novo PD
# RQ2) the moderating effect OSA on hippocampal substructures' volume in de novo PD
#
# Both RQs are addressed from two point of views:
#
# 1) 'classical' multiple regressions adjusting estimates for following covariates - demographics and TIV - assuming effect
#     invariance of these covariates across levels of diagnosis and decisions based on adjusted thresholds for p-values
#     within Neymar-Pearson framework for controlling type I error rates
#
# 2) propensity scores matching based weighted least square multiple regressions adjusting estimates for following
#     covariates - demographics and TIV - assuming effect invariance of these covariates across levels of diagnosis and
#     inference based on g-computation with cluster-robust SEs


rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

# load packages
library(here)
library(tidyverse)
library(MatchIt)
library(marginaleffects) # for inference using g-computation
library(performance)

theme_set( theme_bw() )

# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# read helper files with variable names
subco <- read.csv( here("helpers","subcortical.csv"), sep = ",") # subcortical structures
hippo <- read.csv( here("helpers","hippocampus.csv"), sep = ",") # hippocampal structures

# in-house functions
source( here("scripts","utils.R") )


# PRE-PROCESSING  ----

# read the data
d0 <-
  read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>%
  filter(event == "enrollment")

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
    across( all_of(hippo$name), ~ as.numeric( scale(.x) ) )
  ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
    contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
    contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  } )


## ---- Propensity score matching -----

# fit a logistic regression based propensity score model based on sex and age for base models
# and age, gender and BMI for 'robustness check'
m.out <- lapply(
  
  setNames( c("AGE + GENDER","AGE + GENDER + BMI"), c("AGE + GENDER","AGE + GENDER + BMI") ),
  function(x)
    
    matchit(
      formula = as.formula( paste0("PD ~ ",x) ),
      data = df,
      method = "full",
      distance = "glm",
      link = "logit"
    )
  
)

# matching looks to have worked well with these data and variables
plot(m.out$`AGE + GENDER`, type = "density", interactive = F)
plot(m.out$`AGE + GENDER + BMI`, type = "density", interactive = F)

# prepare a data set with propensity scores weights
mdata <- lapply( setNames( names(m.out), names(m.out) ), function(x) match.data(m.out[[x]]) )


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
      if_else(AHI.F == "H", "AHI ≥ 15", "AHI < 15"),
      levels = c("AHI < 15","AHI ≥ 15"),
      ordered = T
    )
  ) %>%
  
  ggplot() +
  aes(y = volume, x = Diagnosis) +
  geom_boxplot( aes(fill = `OSA: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
  geom_dotplot( aes(fill = `OSA: `), binaxis = "y", stackdir = "center", position = position_dodge(.7), dotsize = 1.5 ) +
  labs( y = bquote("Standardized volume "(mm^3) ) ) +
  ggh4x::facet_grid2( structure ~ side, scales = "free_y", independent = "y" ) +
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
  predictor0 = c("SUBJ + AGE + GENDER + SBTIV", "SUBJ * AGE + GENDER + SBTIV","SUBJ * (AGE + BMI) + GENDER + SBTIV"),
  predictor1 = c("SUBJ * AHI.F + AGE + GENDER + SBTIV", "SUBJ * AHI.F * AGE + GENDER + SBTIV","SUBJ * AHI.F * (AGE + BMI) + GENDER + SBTIV"),
  weight = c("AGE + GENDER","AGE + GENDER","AGE + GENDER + BMI")
  
)

# loop through types of regressions (base vs interaction), structures and types of analysis
lapply(
  
  1:nrow(preds),
  function(r)
    
    lapply(
      
      c("subco","hippo"),
      function(i) {
        
        with(
          
          get(i), {
            
            ### ---- classical regressions with threshold-based decisions ----
            
            fit0 <- fit_reg(df, name, X = preds[r, "predictor0"], w = F) # main effect of disease
            fit1 <- fit_reg(df, name, X = preds[r, "predictor1"], w = F) # main effect of OSA & disease/OSA interaction
            
            # do full analysis for subcortical structures and interaction only for hippocampal substructures
            if (i == "subco") write.table(
              
              # extract and write main effects and interactions for subcortical structures
              x = left_join(
                
                rbind.data.frame( lm_coeff(fit0, term = "SUBJ1"), lm_coeff(fit1, term = "AHI.F1"), lm_coeff(fit1, term = "SUBJ1:AHI.F1") ),
                rbind.data.frame( lm_dia(fit0), lm_dia(fit1) ),
                by = c("y","X")
                
              ) %>% mutate( sig_FDR = bh_adjust(`p value`) ), # re-calculate Benjamini-Hochberg adjusted significance statements
              
              file = here( "tables", paste0(i, "_", preds[r, "model"], "_classical_regressions.csv") ),
              sep = ",", row.names = F, quote = F
              
            ) else if (i == "hippo") write.table(
              
              # extract and write 'main effects and interactions only for hippocampal substructures
              x = left_join( lm_coeff(fit1, term = "SUBJ1:AHI.F1") ,lm_dia(fit1), by = c("y","X") ),
              file = here( "tables", paste0(i, "_" , preds[r, "model"], "_classical_regressions.csv") ),
              sep = ",", row.names = F, quote = F
            )
            
            
            ### ---- weighted regressions with g-computation ----
            
            fit0 <- fit_reg(mdata[[preds[r,"weight"]]], name, X = preds[r, "predictor0"], w = T) # main effect of disease
            fit1 <- fit_reg(mdata[[preds[r,"weight"]]], name, X = preds[r, "predictor1"], w = T) # main effect of OSA & disease/OSA interaction
            
            # weighted regressions with g-computation
            write.table(
              
              # g-computation results
              x = left_join(
                meff( fit = fit1, fit0 = fit0, type = ifelse(i == "subco", "full", "moderation") ), # marginalised effects
                rbind.data.frame( lm_dia(fit0), lm_dia(fit1) ), # diagnostics
                by = c("y","X")
              ),
              
              file = here( "tables", paste0(i, "_", preds[r, "model"], "_weighted_regressions.csv") ),
              sep = ",", row.names = F, quote = F
            )
            
          }
        )
        
      }
    )
  
)


# FOREST PLOTS ----

## ---- subcortical structures regressions ----

read.csv(here("tables","subco_base_weighted_regressions.csv"), sep = ",") %>%
  
  # prepare variables
  mutate(
    
    # group the estimated is conditioned on
    `Group: ` = case_when(
      contrast == "mean(L) - mean(H)" & SUBJ == "CON" ~ "HC",
      contrast == "mean(L) - mean(H)" & SUBJ == "PD" ~ "PD",
      term == "PD - CON" ~ "PD - HC"
    ),
    
    # estimand of interest
    `Estimand: ` = case_when(
      `Group: ` %in% c("HC","PD") ~ "Simple main effect",
      `Group: ` == "PD - HC" ~ "Interaction"
    ),
    
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
  select(y, side, structure, `Group: `, `Estimand: `, estimate, conf.low, conf.high) %>%
  
  # plotting proper
  ggplot() +
  aes(y = estimate, ymin = conf.low, ymax = conf.high, x = structure, shape = `Group: `, colour = `Estimand: `) +
  geom_point(position = position_dodge(width = .3), size = 3.3) + # point estimates
  geom_linerange(position = position_dodge(width = .3), linewidth = 1) + # 95% CIs
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") + # zero for reference
  coord_flip() + # use flip if decided to cut the range of values shown 
  facet_wrap( ~ side, ncol = 2 ) + # a column per hemisphere
  scale_colour_manual( values = c("red","navyblue") ) + # simple effects blue, interaction red
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

# extract structures ordered by s-value (from largest to smallest) on the PD * OSA interaction
ord <- read.csv(here("tables","hippo_base_weighted_regressions.csv"), sep = ",") %>%
  filter(term == "PD - CON") %>%
  arrange( desc(s.value) ) %>%
  mutate( y = sub("_[^_]+$", "", y) ) %>%
  select(y) %>%
  unlist(use.names = F) %>%
  unique()

# do the plotting
read.csv(here("tables","hippo_base_weighted_regressions.csv"), sep = ",") %>%
  
  # prepare variables
  mutate(
    
    # group the estimated is conditioned on
    `Group: ` = case_when(
      contrast == "mean(L) - mean(H)" & SUBJ == "CON" ~ "HC",
      contrast == "mean(L) - mean(H)" & SUBJ == "PD" ~ "PD",
      term == "PD - CON" ~ "PD - HC"
    ),
    
    # estimand of interest
    `Estimand: ` = case_when(
      `Group: ` %in% c("HC","PD") ~ "Simple main effect",
      `Group: ` == "PD - HC" ~ "Interaction"
    ),
    
    # hemisphere of the outcome variable
    side = case_when(
      grepl("_lhx", y) ~ "Left",
      grepl("_rhx", y) ~ "Right"
    ),
    
    # outcome variable brain structure
    structure = factor(
      sub("_[^_]+$", "", y),
      levels = rev(ord),
      ordered = T
    )
    
  ) %>%
  
  # keep only rows and columns of interest/use
  filter( complete.cases(`Group: `) ) %>%
  select(y, side, structure, `Group: `, `Estimand: `, estimate, conf.low, conf.high) %>%
  
  # plotting proper
  ggplot() +
  aes(y = estimate, ymin = conf.low, ymax = conf.high, x = structure, shape = `Group: `, colour = `Estimand: `) +
  geom_point(position = position_dodge(width = .5), size = 3) + # point estimates
  geom_linerange(position = position_dodge(width = .5), linewidth = 1) + # 95% CIs
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") + # zero for reference
  coord_flip() + # use flip if decided to cut the range of values shown 
  facet_wrap( ~ side, ncol = 2 ) + # a column per hemisphere
  scale_colour_manual( values = c("red","navyblue") ) + # simple effects blue, interaction red
  labs(x = NULL, y = "mean(OSA-) - mean(OSA+)") + # estimate is that of difference between OSA- and OSA+ expected means
  theme(legend.position = "right")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","hippo_forest.jpg"),
  dpi = 300,
  width = 10,
  height = 13
)




# table for the classical analysis
fit_reg(df, subco$name, X = "SUBJ * AHI.F * (AGE + BMI) + GENDER + SBTIV", w = F) %>%
  lm_coeff() %>%
  mutate(
    Side = unlist( sapply( 1:nrow(.), function(i) with(subco, side[name == y[i]] ) ) , use.names = F),
    Structure = unlist( sapply( 1:nrow(.), function(i) with(subco, structure[name == y[i]] ) ) , use.names = F),
    Estimate = rprint(Estimate),
    SE = rprint(`Std. Error`),
    `95% CI` = paste0("[", rprint(`2.5 %`),", ",rprint(`97.5 %`), "]"),
    t = rprint(`t value`, 3),
    p = zerolead(`p value`),
    s = rprint(`s value`),
    sig = sig_FDR
  ) %>%
  select(Side, Structure, Estimate, SE, `95% CI`, t, p, s, sig) %>%
  pivot_wider( names_from = Side, values_from = c("Estimate","SE","95% CI","t","p","s","sig") ) %>%
  gt() %>%
  tab_spanner( label = "Left hemisphere", columns = ends_with("Left") ) %>%
  tab_spanner( label = "Right hemisphere", columns = ends_with("Right") ) %>%
  cols_align( align = "center", columns = -1) %>%
  cols_label(
    starts_with("Estimate_") ~ "{{:beta:}}",
    starts_with("SE_") ~ "SE",
    starts_with("95%") ~ "95% CI",
    starts_with("t_") ~ "t value",
    starts_with("p_") ~ "p-value",
    starts_with("s_") ~ "s-value",
    starts_with("sig") ~ "sig."
  )

# table for the g-comparisons
fit_reg(df.matched, subco$name, X = "SUBJ * AHI.F + AGE + GENDER + SBTIV", w = T) %>%
  meff() %>%
  mutate(
    Contrast = "OSA- vs OSA+",
    Group = case_when(term == "PD - CON" ~ "PD vs HC", SUBJ == "PD" ~ "PD", SUBJ == "CON" ~ "HC"),
    Side = unlist( sapply( 1:nrow(.), function(i) with(subco, side[name == y[i]] ) ) , use.names = F),
    Structure = unlist( sapply( 1:nrow(.), function(i) with(subco, structure[name == y[i]] ) ) , use.names = F),
    Estimate = rprint(estimate),
    SE = rprint(std.error),
    `95% CI` = paste0("[", rprint(conf.low),", ",rprint(conf.high), "]"),
    Z = rprint(statistic, 3),
    p = zerolead(p.value),
    s = rprint(s.value)
  ) %>%
  select(Contrast, Group, Side, Structure, Estimate, SE, `95% CI`, Z, p, s) %>%
  pivot_wider( names_from = Side, values_from = -c("Contrast","Group","Structure","Side") ) %>%
  gt(groupname_col = "Structure") %>%
  tab_spanner( label = "Left hemisphere", columns = ends_with("Left") ) %>%
  tab_spanner( label = "Right hemisphere", columns = ends_with("Right") ) %>%
  cols_align( align = "center", columns = -1) %>%
  cols_label(
    starts_with("Estimate_") ~ "{{:beta:}}",
    starts_with("SE_") ~ "SE",
    starts_with("95%") ~ "95% CI",
    starts_with("Z_") ~ "Z value",
    starts_with("p_") ~ "p-value",
    starts_with("s_") ~ "s-value"
  )

