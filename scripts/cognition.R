# Computes regressions evaluating:
#
# RQ3) the moderating effect OSA on cognitive performance in de novo PD
#
# The RQ is addressed via a series Bayesian test/re-test GLMMs with marginal means postprocessing of posteriors


rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

# load packages
library(here)
library(tidyverse)
library(brms)
library(priorsense)
library(bayestestR)
library(bayesplot)
library(patchwork)

theme_set( theme_bw() )
color_scheme_set("viridisA") # colour scheme

sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders
psych <- read.csv( here("helpers","psychs.csv"), sep = ";") # read helper files with variable names

# in-house functions
source( here("scripts","utils.R") )


# PRE-PROCESSING  ----

# read the data
d0 <-
  
  # read and pre-process variables
  read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>%
  mutate(
    id = Study.ID,
    sex = factor( case_when(GENDER == 1 ~ "male", GENDER == 0 ~ "female") ),
    group = factor( case_when(SUBJ == "PD" ~ "PD", SUBJ == "CON" ~ "HC") ),
    osa = factor( case_when(AHI.F == "H" ~ "+", AHI.F == "L" ~ "-") ),
    event = factor(
      case_when(event == "enrollment" ~ "enrollment", event == "followup_4" ~ "retest"),
      levels = c("enrollment", "retest"),
      ordered = T
    ),
    age = AGE,
    edu = EDU.Y
  ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(sex) <- contr.equalprior_pairs(2)
    contrasts(group) <- contr.equalprior_pairs(2)
    contrasts(osa) <- contr.equalprior_pairs(2)
    contrasts(event) <- contr.equalprior_pairs(2)
  } )


# prepare data set for analyses
df <- d0 %>%
  
  select(id, event, sex, group, osa, event, age, edu, all_of(psych$variable[-13]) ) %>%
  filter( complete.cases(event) ) %>%
  mutate(
    across(
      all_of( subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable ),
      ~ -log(.x)
    )
  )


# DO OSA+/- PATIENTS SHOW DISTINCTIVE RE-TEST COGNITIVE PROFILE COMPARED TO HEALTHY CONTROLS? ----

# amount of data points available
lapply(
  
  setNames( psych$variable[-13], psych$variable[-13] ),
  function(y)
    
    table( subset(df, complete.cases( get(y) ) )$id) %>%
    as.data.frame() %>%
    mutate( Group =  unlist( sapply( 1:nrow(.), function(j) with(df, unique( grp[id == Var1[j]] ) ) ) ) ) %>%
    mutate( OSA = unlist( sapply( 1:nrow(.), function(j) with(df, unique( osa[id == Var1[j]] ) ) ) ) ) %>%
    select(OSA, Group, Freq) %>%
    table()
  
)

# extract scaling values, i.e.,
# enrollment full sample means and SDs
scl <- sapply(
  
  c(psych$variable[-13], "age", "edu"),
  function(i)
    
    with(subset(df, event == "enrollment"), c(M = mean( get(i), na.rm = T ), SD = sd( get(i), na.rm = T ) ) )
  
) %>% t()

# scale the variables
df <- df %>%
  
  mutate(
    across(
      all_of( rownames(scl) ),
      ~ (.x - scl[cur_column(), "M"]) / scl[cur_column(), "SD"],
      .names = "{.col}_sc"
    )
  )

# set-up priors
prior1 <- c(
  
  prior("normal(0, 1)", class = "Intercept"),
  prior("normal(0, 1)", class = "b"),
  prior("exponential(1)", class = "sd", group = "id"),
  prior("exponential(1)", class = "sigma")
  
)

# fit the models
fit1 <- lapply(
  
  with( psych, setNames( paste0(psych$variable[-13],"_sc"), psych$variable[-13]) ),
  function(y) {
    
    data <- df %>% mutate( y = get(y) )
    
    brm(
      formula = bf( y ~ 1 + event * grp * osa + (1 | id) ),
      data = data,
      prior = prior1,
      family = gaussian(link = "identity"),
      seed = 87542
    )
    
  }
)

# do some posterior predictive checks
plot_ppc_stat(fit1, d2, stat = "mean")
plot_ppc_stat(fit1, d2, stat = "median")
plot_ppc_stat(fit1, d2, stat = "sd")

# save the results
write.table(
  
  x = do.call(
    
    rbind.data.frame,
    lapply(
      
      names(fit2),
      function(y)
        
        left_join(
          
          describe_posterior(fit2[[y]], test = "p_direction") %>%
            as.data.frame() %>%
            mutate(Outcome = y, .before = 1),
          
          powerscale_sensitivity(fit2[[y]]) %>%
            as.data.frame() %>%
            rename(
              "Parameter" = "variable",
              "prior_sensitivity" = "prior",
              "likelihood_sensitivity" = "likelihood",
              "powerscale_diagnosis" = "diagnosis"
            ),
          
          by = "Parameter"
          
        )
      
    )
  ),
  
  file = here("tables","cognition_retest.csv"),
  sep = ",",
  row.names = F,
  quote = F
  
)



# RQ3: DOES OSA AFFECT COGNITION DIFFERENTLY IN PD PATIENTS COMPARED TO HEALTHY CONTROLS? ----

# extract non-MoCA variables
psychvar <- with( psych, variable[variable != "moca"] )
psychlab <- with( psych, label[variable != "moca"] )

## Interaction boxplots ----

d0 %>%
  
  select(SUBJ, AHI.F, all_of(psychvar) ) %>%
  pivot_longer(
    cols = all_of(psychvar),
    values_to = "score",
    names_to = "test"
  ) %>%
  mutate(
    test = factor(
      unlist(
        sapply( 1:nrow(.), function(i) with( psych, label[variable == test[i]] ) ),
        use.names = F
      ),
      levels = psychlab,
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
  width = 8,
  height = 10
)


## ---- Propensity scores matching ----

# fit a logistic regression based propensity score model based on sex and age
m.out <- matchit(
  formula = PD ~ EDU.Y + GENDER, # this one works better than AGE + EDU.Y + GENDER
  data = df,
  method = "full",
  distance = "glm"
)

# matching looks to have worked well with these data and variables
plot(m.out, type = "density", interactive = F)

# prepare a data set with propensity scores weights
df.matched <- match.data(m.out)


## ---- Linear regressions ----

### ---- classical regressions with threshold-based decisions ----

write.table(
  
  # extract and write main effects and interactions for subcortical structures
  x = left_join(
    
    rbind.data.frame(
      lm_coeff(fit_reg(df, psychvar, X = "SUBJ * AGE + GENDER + EDU.Y", w = F), term = "SUBJ1"), # main effect of disease
      lm_coeff(fit_reg(df, psychvar, X = "SUBJ * AHI.F * AGE + GENDER + EDU.Y", w = F), term = "AHI.F1"), # main effect of OSA
      lm_coeff(fit_reg(df, psychvar, X = "SUBJ * AHI.F * AGE + GENDER + EDU.Y", w = F), term = "SUBJ1:AHI.F1") # disease/OSA interaction
    ),
    rbind.data.frame(
      lm_dia( fit_reg(df, psychvar, X = "SUBJ * AGE + GENDER + EDU.Y", w = F) ),
      lm_dia( fit_reg(df, psychvar, X = "SUBJ * AHI.F * AGE + GENDER + EDU.Y", w = F) )
    ),
    by = c("y","X")
    
  ) %>% mutate( sig_FDR = bh_adjust(`p value`) ), # re-calculate Benjamini-Hochberg adjusted significance statements
  
  file = here("tables","cognition_classical_regressions.csv"),
  sep = ",", row.names = F, quote = F
  
)


### ---- weighted regressions with g-computation ----

write.table(
  
  # g-computation results
  x = left_join(
    
    # marginalised effects
    meff(
      fit = fit_reg(df.matched, psychvar, X = "SUBJ * AHI.F * AGE + GENDER + EDU.Y", w = T), # main effect of disease
      fit0 = fit_reg(df.matched, psychvar, X = "SUBJ * AGE + GENDER + EDU.Y", w = T), # main effect of OSA & disease/OSA interaction
      type = "full"
    ),
    
    # diagnostics
    rbind.data.frame(
      lm_dia( fit_reg(df.matched, psychvar, X = "SUBJ * AGE + GENDER + EDU.Y", w = T) ),
      lm_dia( fit_reg(df.matched, psychvar, X = "SUBJ * AHI.F * AGE + GENDER + EDU.Y", w = T) )
    ),
    
    by = c("y","X")
    
  ),
  
  file = here("tables","cognition_weighted_regressions.csv"),
  sep = ",", row.names = F, quote = F
  
)


### ---- Bayesian regression with heteroscedasticity ----

# prepare formulas
formulas <- list(
  
  # formulas for base models (assuming homoscedasticity)
  varequal = lapply(
    setNames( c("tmt_b","gpt_phk","gpt_lhk"), c("tmt_b","gpt_phk","gpt_lhk") ),
    function(i)
      list(
        SUBJ = paste0(i," ~ SUBJ * AGE + GENDER + EDU.Y") %>% as.formula() %>% bf(),
        AHI = paste0(i," ~ SUBJ * AHI.F * AGE + GENDER + EDU.Y") %>% as.formula() %>% bf()
      )
  ),
  
  # formulas for variance adjusted models (allowing heteroscedasticity)
  # listing them one by one because as.formula do not want work with commas
  heteroscedastic = list(
    
    tmt_b = list(
      SUBJ = bf(tmt_b ~ SUBJ * AGE + GENDER + EDU.Y, sigma ~ SUBJ),
      AHI = bf(tmt_b ~ SUBJ * AHI.F * AGE + GENDER + EDU.Y, sigma ~ SUBJ)
    ),
    
    gpt_phk = list(
      SUBJ = bf(gpt_phk ~ SUBJ * AGE + GENDER + EDU.Y, sigma ~ SUBJ * AGE ),
      AHI = bf(gpt_phk ~ SUBJ * AHI.F * AGE + GENDER + EDU.Y, sigma ~ SUBJ * AGE )
    ),
    
    gpt_lhk = list(
      SUBJ = bf(gpt_lhk ~ SUBJ * AGE + GENDER + EDU.Y, sigma ~ SUBJ + AGE ),
      AHI = bf(gpt_lhk ~ SUBJ * AHI.F * AGE + GENDER + EDU.Y, sigma ~ SUBJ + AGE )
    )
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
        
        lapply(
          
          setNames( names(formulas[[i]][[j]]), names(formulas[[i]][[j]]) ),
          function(k)
            
            brm( formula = formulas[[i]][[j]][[k]],
                 data = df,
                 family = gaussian(link = "identity"),
                 prior = priors[[i]],
                 seed = 87542
            )
          
        )
    )
)


#### ----- posterior predictive checks ----

lapply(
  
  c("dens_overlay_grouped","stat_grouped"),
  function(k) {
    
    ppc <- lapply(
      
      setNames( names(fit), names(fit) ),
      function(i)
        
        lapply(
          
          setNames( names(fit[[i]]), names(fit[[i]]) ),
          function(y)
            
            lapply(
              
              setNames( names(fit[[i]][[y]]), names(fit[[i]][[y]]) ),
              function(x) {
                
                # set distinct colour palletes for different types of models
                if (i == "varequal") bayesplot::color_scheme_set( scheme = "red" )
                else if (i == "heteroscedastic") bayesplot::color_scheme_set( scheme = "blue" )
                
                set.seed(87542) # seed for reproducibility
                
                # plotting proper
                pp_check( fit[[i]][[y]][[x]], ndraws = 100, type = k, group = "SUBJ", stat = "sd", bins = 15 ) +
                  theme_bw(base_size = 14) +
                  labs(
                    x = paste0( with(psych, label[variable == y]), " (-log seconds)" ),
                    title = case_when(
                      i == "varequal" ~ paste0( as.character(formulas[[i]][[y]][[x]])[1], ", sigma ~ 1" ),
                      i == "heteroscedastic" ~ paste0(
                        as.character(formulas[[i]][[y]][[x]])[1] , ", ",
                        sub( ")", "", sub( "list(sigma = ", "", as.character(formulas[[i]][[y]][[x]])[2], fixed = T ) )
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
    )
    
    # put them to a single plot
    with(
      ppc,
      ( varequal$tmt_b$SUBJ | heteroscedastic$tmt_b$SUBJ | varequal$tmt_b$AHI | heteroscedastic$tmt_b$AHI ) /
        ( varequal$gpt_phk$SUBJ | heteroscedastic$gpt_phk$SUBJ | varequal$gpt_phk$AHI | heteroscedastic$gpt_phk$AHI ) /
        ( varequal$gpt_lhk$SUBJ | heteroscedastic$gpt_lhk$SUBJ | varequal$gpt_lhk$AHI | heteroscedastic$gpt_lhk$AHI )
    )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("cognition_ppc_", sub("_.*","",k), ".jpg") ),
      dpi = 300,
      width = 20,
      height = 15
    )
    
  }
)


#### ----- coefficients extraction ----

write.table(
  
  x = left_join(
    
    # extract model coefficients
    lapply(
      
      names(fit),
      function(i) {
        
        fit0 <- list(tmt_b = fit[[i]]$tmt_b$SUBJ, gpt_phk = fit[[i]]$gpt_phk$SUBJ, gpt_lhk = fit[[i]]$gpt_lhk$SUBJ)
        fit1 <- list(tmt_b = fit[[i]]$tmt_b$AHI, gpt_phk = fit[[i]]$gpt_phk$AHI, gpt_lhk = fit[[i]]$gpt_lhk$AHI)
        
        do.call(
          
          rbind.data.frame,
          list(
            lm_coeff(fit0, term = "SUBJ1", type = "Bayesian"),
            lm_coeff(fit1, term = "AHI.F1", type = "Bayesian"),
            lm_coeff(fit1, term = "SUBJ1:AHI.F1", type = "Bayesian")
          )
          
        )
      }
    ) %>% do.call( rbind.data.frame, . ),
    
    # extract diagnostics
    lapply(
      
      names(fit),
      function(k)
        
        lapply(
          
          names(fit[[k]]),
          function(y) {
            
            # print prior sensitivity diagnostics
            print( powerscale_sensitivity(fit[[k]][[y]][[1]]) )
            print( powerscale_sensitivity(fit[[k]][[y]][[2]]) )
            
            # return classical diagnostics
            return(
              lm_dia(fit[[k]][[y]], type = "Bayesian") %>%
                mutate(y = y, .before = 1) %>%
                mutate( p_breusch_pagan = zerolead(p_breusch_pagan) )
            )
            
          }
          
        ) %>% do.call( rbind.data.frame, . )
      
    ) %>% do.call( rbind.data.frame, . ),
    
    # join by outcome, predictor, and variance terms
    by = c("y","X","sigma")
    
  ),
  
  file = here("tables","cognition_bayesian_regressions.csv"),
  sep = ",", row.names = F, quote = F
  
)