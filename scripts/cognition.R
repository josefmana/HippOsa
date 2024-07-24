# Computes regressions evaluating:
#
# RQ4) the moderating effect OSA on longitudinal general cognitive performance in de novo PD
# RQ5) the moderating effect OSA on re-test cognitive profile in de novo PD
#
# Both RQs are addressed via Bayesian GLMMs


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


# UTILS ----

# plot posterior predictive stats for a set of models
plot_ppc_stat <- function(fit, data, labs = psych, sleep = 3, stat = "sd") lapply(
  
  names(fit),
  function(i) {
    
    y <- paste0(i,"_sc")
    d <- subset( data, complete.cases( get(y) ) ) %>% mutate(x = paste0(pd,"_",event) )
    
    print(
      pp_check(
        object = d[ , y],
        yrep = posterior_predict(fit[[i]], newdata = d),
        fun = ppc_stat_grouped,
        stat = stat,
        group = d$x
      ) +
        labs(
          title = with(labs , label[variable == i] ),
          subtitle = paste0("Observed (thick line) vs predicted (histogram) ",stat," for pre/post assessment of PD/HC participants")
        ) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5)
        )
    )
    
    Sys.sleep(sleep)
    
  }
)

# PRE-PROCESSING  ----

# read the data
d0 <-
  read.csv( here("_data","primary_dataset.csv"), sep = "," ) %>%
  mutate(
    id = Study.ID,
    sex = factor( case_when(GENDER == 1 ~ "male", GENDER == 0 ~ "female") ),
    pd = factor( case_when(SUBJ == "PD" ~ 1, SUBJ == "CON" ~ 0) ),
    osa = factor( case_when(AHI.F == "H" ~ 1, AHI.F == "L" ~ 0) ),
    age = AGE,
    edu = EDU.Y
  )

# prepare MoCA data set for RQ4
d1 <- d0 %>%
  select(id, event, sex, pd, osa, age, edu, moca) %>%
  filter( complete.cases(moca) ) %>%
  mutate( event_num = as.numeric( as.factor(event) ) - 1 )

# prepare battery data set for RQ5
d2 <- d0 %>%
  select(id, event, sex, pd, osa, age, edu, all_of(psych$variable[-13]) ) %>%
  filter( event %in% c("enrollment","followup_4") ) %>%
  mutate(
    across(
      all_of( subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable ),
      ~ -log(.x)
    )
  )

# check for missing values
sapply(

  1:nrow(df),
  function(i)
    if( sum( is.na(df[i, ]) ) > 0 ) df[i, ]

) %>%
  do.call( rbind.data.frame , . ) %>%
  mutate_if( is.numeric, round, 2 )


# RQ4: DO OSA+/- PATIENTS SHOW DISTINCTIVE LONGITUDINAL COGNITIVE TRAJECTORY COMPARED TO HEALTHY CONTROLS? ----

# plot number of assessments per participant
table(d1$id) %>%
  
  as.data.frame() %>%
  mutate( group =  if_else( grepl("CON",Var1), "HC", "PD" ) ) %>%
  
  ggplot() +
  aes(y = reorder(Var1, Freq), x = Freq, fill = group) +
  geom_bar( stat = "identity" ) +
  labs(y = NULL, x = "Number of assessments") +
  scale_x_continuous( breaks = seq(0,9,1), labels = seq(0,9,1) ) +
  facet_wrap( ~group, scales = "free_y" ) +
  theme_bw() +
  theme(legend.position = "none")

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","cognition_longitudinal_assessments.jpg"),
  dpi = 300,
  width = 8,
  height = 8
)


## ---- longitudinal GLM of MoCA ----

# set-up priors
prior1 <- c(
  
  prior("normal(0, 3)", class = "Intercept"),
  prior("normal(0, 3)", class = "b"),
  prior("normal(0, 3)", class = "sd", group = "id"),
  prior("lkj(2)", class = "cor")
  
)

# fit it
fit1 <- brm(
  
  moca | trials(30) ~ 1 + (1 + event_num | id) + event_num * pd * osa,
  data = d1,
  prior = prior1,
  family = binomial()
  
)

# check for prior sensitivity
powerscale_sensitivity(fit1)

# do some posterior checks
pp_check(fit1, type = "stat_2d", ndraws = NULL )
pp_check(fit1, type = "bars_grouped", ndraws = NULL, group = "event_num" )
pp_check(fit1, type = "bars_grouped", ndraws = NULL, group = "pd" )
pp_check(fit1, type = "bars_grouped", ndraws = NULL, group = "osa" )

# extract and save posterior summaries
describe_posterior(fit1, test = "p_direction") %>%
  as.data.frame() %>%
  write.table(file = here("tables","cognition_longitudinal.csv"), sep = ",", row.names = F, quote = F)


# RQ5: DO OSA+/- PATIENTS SHOW DISTINCTIVE RE-TEST COGNITIVE PROFILE COMPARED TO HEALTHY CONTROLS? ----

# amount of data points available
lapply(
  
  setNames( psych$variable[-13], psych$variable[-13] ),
  function(y)
    
    table( subset(d2, complete.cases( get(y) ) )$id) %>%
    as.data.frame() %>%
    mutate( Group =  if_else( grepl("CON",Var1), 0, 1 ) ) %>%
    mutate( OSA = unlist( sapply( 1:nrow(.), function(j) with(d1, unique( osa[id == Var1[j]] ) ) ) ) ) %>%
    select(OSA, Group, Freq) %>%
    table()
  
)

# extract scaling values, i.e.,
# enrollment full sample means and SDs
scl <- sapply(
  
  c(psych$variable[-13], "age", "edu"),
  function(i)
    
    with(subset(d2, event == "enrollment"), c(M = mean( get(i), na.rm = T ), SD = sd( get(i), na.rm = T ) ) )
  
) %>% t()

# scale the variables
d2 <- d2 %>%
  
  mutate(
    across(
      all_of( rownames(scl) ),
      ~ (.x - scl[cur_column(), "M"]) / scl[cur_column(), "SD"],
      .names = "{.col}_sc"
    )
  )

# set-up priors
prior2 <- c(
  
  prior("normal(0, 1)", class = "Intercept"),
  prior("normal(0, 1)", class = "b"),
  prior("exponential(1)", class = "sd", group = "id"),
  prior("exponential(1)", class = "sigma")
  
)

# fit the models
fit2 <- lapply(
  
  with( psych, setNames( paste0(psych$variable[-13],"_sc"), psych$variable[-13]) ),
  function(y) {
    
    data <- d2 %>% mutate( y = get(y) )
    
    brm(
      formula = bf( y ~ 1 + event * pd * osa + (1 | id) ),
      data = data,
      prior = prior2,
      family = gaussian(link = "identity"),
      seed = 87542
    )
    
  }
)

# do some posterior predictive checks
plot_ppc_stat(fit2, d2, stat = "mean")
plot_ppc_stat(fit2, d2, stat = "median")
plot_ppc_stat(fit2, d2, stat = "sd")

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
          
          powerscale_sensitivity(fit2$avlt_1_5) %>%
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