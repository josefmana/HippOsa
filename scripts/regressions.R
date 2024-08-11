# Computes regressions evaluating the association between OSA and brain structures' volume/cognition:
#
# RQ1) the distribution of subcortical brain structures' volume in de novo PD conditional on OSA
# RQ2) the distribution of  hippocampal substructures' volume in de novo PD conditional on OSA
# RQ3) the distribution of cognitive performance in de novo PD conditional on OSA


rm( list = ls() ) # clear environment

# load packages
library(here)
library(tidyverse)
library(marginaleffects)
library(performance)
library(ggh4x)
library(brms)
#library(priorsense)
library(bayesplot)
library(patchwork)

theme_set( theme_bw() )

# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# read helper files with variable names
subco <- read.csv( here("helpers","subcortical.csv"), sep = ",") # subcortical structures
hippo <- read.csv( here("helpers","hippocampus.csv"), sep = ",") %>% filter( complete.cases(name) ) # hippocampal structures
psych <- read.csv( here("helpers","psychs.csv"), sep = ";") # psychologic variables

# extract response time variable names to be log-transformed
rt_vars <- subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable

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
    across( all_of(rt_vars), ~ -log(.x) ) # cognition
  ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
    contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
    contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  } )

# extract scaling values, i.e.,
# enrollment full sample means and SDs
scl <- sapply(
  
  c("AGE", "EDU.Y", "BMI", "SBTIV", subco$name, unique(hippo$name), psych$variable),
  function(i)
    with(df, c(M = mean( get(i), na.rm = T ), SD = sd( get(i), na.rm = T ) ) )
  
) %>% t()

# scale the variables
df <- df %>% mutate( across( all_of( rownames(scl) ), ~ (.x - scl[cur_column(), "M"]) / scl[cur_column(), "SD"] ) )


# INTERACTION BOXPLOTS ----

## ---- BRAINS ----

# prepare results for p-values reported in the boxplot
p <-
  fit_reg(d = df, outcomes = subco$scaled, X = "SUBJ * AHI.F", w = F) %>%
  mass() %>%
  mutate( sig = bh_adjust(p.value) )

# plot it
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


## ---- COGNTION ----

d0 %>%
  
  select(SUBJ, AHI.F, all_of(psych$variable) ) %>%
  pivot_longer(
    cols = all_of(psych$variable),
    values_to = "score",
    names_to = "test"
  ) %>%
  mutate(
    test = factor(
      unlist(
        sapply( 1:nrow(.), function(i) with( psych, label[variable == test[i]] ) ),
        use.names = F
      ),
      levels = psych$label,
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
  aes(y = score, x = Diagnosis) +
  geom_boxplot( aes(fill = `OSA: `), width = .6, position = position_dodge(.7), linewidth = .75 ) +
  geom_dotplot( aes(fill = `OSA: `), binaxis = "y", stackdir = "center", position = position_dodge(.7), dotsize = 1.5 ) +
  labs( y = "Test score" ) +
  facet_wrap( ~test, scales = "free", nrow = 5 ) +
  scale_fill_manual( values = c("#64CDCC","#F9A729") ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","cognition_boxplots.jpg"),
  dpi = 300,
  width = 9,
  height = 9
)


# LINEAR REGRESSIONS ----

# set-up formulas
forms <- data.frame(
  
  object = c("subco","hippo","psych"),
  y = c("subcortical", "hippocampi","cognition"),
  X = c( rep("SUBJ * AHI.F + AGE + GENDER + SBTIV", 2), "SUBJ * AHI.F + AGE + GENDER + EDU.Y")
  
)

# loop through types of regressions (base vs interaction), and structures
lapply(
  
  1:nrow(forms),
  function(r) with(
      
      get( forms[r, "object"] ), {
        
        # extract outcome
        y <- forms[r, "y"]
        
        # extract variables
        if(y == "cognition") vars <- variable else vars <- unique(name)
        
        # fit the models
        fit <- fit_reg(df, vars, X = forms[r, "X"], w = F)
        
        # regression coefficients with threshold-based decisions
        # do full analysis for subcortical structures and cognition, interaction only for hippocampal substructures
        if (y == "hippocampi") write.table(
          
          # extract and write 'main effects and interactions only for hippocampal substructures
          x = left_join( lm_coeff(fit, term = "SUBJ1:AHI.F1") ,lm_dia(fit), by = c("y","X") ),
          file = here( "tables", paste0(y,"_regression_coefficients.csv") ),
          sep = ",",
          row.names = F,
          quote = F

        ) else  write.table(
          
          # extract and write main effect coefficients and interactions for subcortical structures
          x = left_join(
            
            rbind.data.frame(
              lm_coeff(fit, term = "SUBJ1"),
              lm_coeff(fit, term = "AHI.F1"),
              lm_coeff(fit, term = "SUBJ1:AHI.F1")
            ),
            lm_dia(fit),
            by = c("y","X")
            
          ) %>% mutate( sig_FDR = bh_adjust(`p value`) ), # re-calculate Benjamini-Hochberg adjusted significance statements
          
          file = here( "tables", paste0(y,"_regression_coefficients.csv") ),
          sep = ",",
          row.names = F,
          quote = F
          
        )

        # marginal associations
        write.table(
          
          x = left_join(
            mass( fit = fit, fit0 = fit, type = ifelse(y == "hippocampi", "moderation", "fullest") ), # note that the "fullest" option contains some redundancy in its interactions, it's there for completeness
            lm_dia(fit), # diagnostics
            by = c("y","X")
          ),
          
          file = here( "tables", paste0(y,"_marginals.csv") ),
          sep = ",",
          row.names = F,
          quote = F

        )
      }
  )
)


# FOREST PLOTS ----

## ---- SUBCORTICAL STRUCTURES REGRESSIONS ----

left_join(
  
  # the marginal associations
  read.csv(here("tables","subcortical_marginals.csv"), sep = ","),
  
  # add BH adjusted significance statements
  read.csv(here("tables","subcortical_marginals.csv"), sep = ",") %>%
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
  scale_colour_manual( values = c("grey","orange") ) + # mark p < threshold
  labs(x = NULL, y = "mean(OSA-) - mean(OSA+)") + # estimate is that of difference between OSA- and OSA+ expected means
  theme(legend.position = "right")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","subcortical_forest.jpg"),
  dpi = 300,
  width = 8,
  height = 8
)


## ---- HIPPOCAMPAL STRUCTURES REGRESSIONS ----

# prepare order of structures
ord <-
  hippo[ c("structure","order") ] %>%
  unique() %>%
  arrange(order) %>%
  select(structure) %>%
  unlist(use.names = F)

# plot it
left_join(
  
  # marginal associations
  read.csv(here("tables","hippocampi_marginals.csv"), sep = ","),
  
  # add BH adjusted significance statements for interactions
  read.csv(here("tables","hippocampi_marginals.csv"), sep = ",") %>%
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
  scale_colour_manual( values = c("grey","orange") ) + # mark p < threshold
  labs(x = NULL, y = "mean(OSA-) - mean(OSA+)") + # estimate is that of difference between OSA- and OSA+ expected means
  theme(legend.position = "right")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","hippocampi_forest.jpg"),
  dpi = 300,
  width = 8,
  height = 12.5
)


## ---- COGNITIVE INDEXES REGRESSIONS ----

# use different colour because these are different comparisons
left_join(
  
  # marginal associations
  read.csv(here("tables","cognition_marginals.csv"), sep = ","),
  
  # BH adjusted significance statements for the main effects
  read.csv(here("tables","cognition_marginals.csv"), sep = ",") %>%
    filter(term == "SUBJ" & is.na(AHI.F) ) %>%
    mutate( sig = bh_adjust(p.value) ) %>%
    select(y, sig)
  
) %>%
  
  # prepare variables
  mutate(
    
    # group the estimated is conditioned on
    `OSA: ` = case_when(
      contrast == "mean(PD) - mean(CON)" & AHI.F == "L" ~ "OSA-",
      contrast == "mean(PD) - mean(CON)" & AHI.F == "H" ~ "OSA+",
      term == "H - L" ~ "OSA+ - OSA-"
    ),
    
    # significance statement
    `Sig. (5% FDR):` = if_else(sig == "*", T, F),
    
    # hemisphere of the outcome variable
    domain = factor(
      unlist( sapply( 1:nrow(.), function(i) with( psych, domain[variable == y[[i]]] ) ), use.names = F),
      levels = unique(psych$domain),
      ordered = T
    ),
    
    # outcome variable brain structure
    index = factor(
      unlist( sapply( 1:nrow(.), function(i) with( psych, label[variable == y[[i]]] ) ), use.names = F),
      levels = rev(psych$label),
      ordered = T
    )
    
  ) %>%
  
  # keep only rows and columns of interest/use
  filter( complete.cases(`OSA: `) ) %>%
  select(y, index, domain, `OSA: `, `Sig. (5% FDR):`, estimate, conf.low, conf.high) %>%
  
  # plotting proper
  ggplot() +
  aes(y = estimate, ymin = conf.low, ymax = conf.high, x = index, shape = `OSA: `, colour = `Sig. (5% FDR):`) +
  geom_point(position = position_dodge(width = .5), size = 3) + # point estimates
  geom_linerange(position = position_dodge(width = .5), linewidth = 1) + # 95% CIs
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") + # zero for reference
  coord_flip() + # use flip if decided to cut the range of values shown 
  facet_wrap( ~ domain, scale = "free_y", ncol = 2) + # a column per hemisphere
  scale_colour_manual( values = c("navyblue","red") ) + # mark p < threhold for main effects red
  labs(x = NULL, y = "mean(PD) - mean(HC)") + # estimate is that of difference between PD and HC expected means
  theme(legend.position = "right")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","cognition_forest.jpg"),
  dpi = 300,
  width = 8,
  height = 10
)



# BAYESIAN REGRESSIONS WITH HETEROSCEDASTICITY ----

# prepare formulas
formulas <- list(
  
  # formulas for base models (assuming homoscedasticity)
  varequal = lapply(
    
    setNames( c("tmt_b","gpt_phk","gpt_lhk"), c("tmt_b","gpt_phk","gpt_lhk") ),
    function(i)
      paste0(i," ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y") %>%
      as.formula() %>%
      bf()

  ),
  
  # formulas for variance adjusted models (allowing heteroscedasticity)
  # listing them one by one because as.formula do not want work with commas
  heteroscedastic = list(
    
    tmt_b = bf(tmt_b ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y, sigma ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y),
    gpt_phk = bf(gpt_phk ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y, sigma ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y),
    gpt_lhk = bf(gpt_lhk ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y, sigma ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y)

  )
  
)

# re-fit the models for TMT-B and GPT via brms with differing variance for the groups
priors <- list(
  
  # priors for base models (assuming homoscedasticity)
  varequal = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(1), class = sigma)
  ),
  
  # priors for variance adjusted models (allowing heteroscedasticity)
  heteroscedastic = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 1), class = b),
    prior(normal(0, 1), class = Intercept, dpar = sigma),
    prior(normal(0, 1), class = b, dpar = sigma)
  )
  
)

# fit the models
fit <- lapply(
  
  setNames( names(formulas), names(formulas) ),
  function(i)
    
    lapply(
      
      setNames( names(formulas[[i]]), names(formulas[[i]]) ),
      function(j)
        
        brm( formula = formulas[[i]][[j]],
             data = df,
             family = gaussian(link = "identity"),
             prior = priors[[i]],
             seed = 87542
        )
        
    )
)


## ---- POSTERIOR PREDICTIVE CHECKS ----

# check SDs
sapply( names(fit), function(i) sapply( c("group", "osa", "group_osa"), function(x) plot_ppc_stat(fit[[i]], df, var = x, sleep = 3) ) )

# check prior-sensitivity
sapply( names(fit), function(i) sapply( names(fit[[i]]), function(y) print( powerscale_sensitivity(fit[[i]][[y]]) ) ) ) # a-ok

# save some
lapply(
  
  c("dens_overlay_grouped","stat_grouped"),
  function(k) {
    
    ppc <- lapply(
      
      setNames( names(fit), names(fit) ),
      function(i)
        
        lapply(
          
          setNames( names(fit[[i]]), names(fit[[i]]) ),
          function(y) {
            
            # set distinct colour palletes for different types of models
            if (i == "varequal") bayesplot::color_scheme_set( scheme = "red" )
            else if (i == "heteroscedastic") bayesplot::color_scheme_set( scheme = "blue" )
            
            set.seed(87542) # seed for reproducibility
                
            # plotting proper
            pp_check(fit[[i]][[y]], ndraws = 100, type = k, group = "SUBJ", stat = "sd", bins = 12) +
              theme_bw(base_size = 14) +
              labs(
                x = paste0( with(psych, label[variable == y]), " (-log seconds)" ),
                title = case_when(
                  i == "varequal" ~ paste0( as.character(formulas[[i]][[y]])[1], ", sigma ~ 1" ),
                  i == "heteroscedastic" ~ paste0(
                    as.character(formulas[[i]][[y]])[1] , ", ",
                    sub( ")", "", sub( "list(sigma = ", "", as.character(formulas[[i]][[y]])[2], fixed = T ) )
                  )
                )
              ) +
              theme(
                legend.position = "none",
                panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 10)
              )
            }
        )
    )
    
    # put them to a single plot
    with(
      ppc,
      ( varequal$tmt_b | heteroscedastic$tmt_b ) /
        ( varequal$gpt_phk | heteroscedastic$gpt_phk ) /
        ( varequal$gpt_lhk | heteroscedastic$gpt_lhk )
    )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("cognition_ppc_", sub("_.*","",k), ".jpg") ),
      dpi = 300,
      width = 15,
      height = 15
    )
    
  }
)


## ---- COEFFICIENTS EXTRACTION ----

write.table(
  
  # extract model coefficients
  x = lapply(
    
    names(fit),
    function(i) do.call(
      
      rbind.data.frame,
      list(
        lm_coeff(fit[[i]], term = "SUBJ1", type = "Bayesian"),
        lm_coeff(fit[[i]], term = "AHI.F1", type = "Bayesian"),
        lm_coeff(fit[[i]], term = "SUBJ1:AHI.F1", type = "Bayesian")
      )
      
    )
    
  ) %>% do.call( rbind.data.frame, . ),
  
  file = here("tables","cognition_bayesian_regressions.csv"),
  sep = ",",
  row.names = F,
  quote = F
  
)
