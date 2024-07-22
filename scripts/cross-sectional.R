# Computes regressions evaluating:
#
# RQ1) the moderating effect OSA on subcortical brain structures' volume in PD
# RQ2) the moderating effect OSA on hippocampal substructures' volume in PD
# RQ3) the moderating effect OSA on cognitive performance in PD
#
# All RQs are addressed from two point of views:
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
library(brms)
library(priorsense)
library(patchwork)

theme_set( theme_bw() )

# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# read helper files with variable names
psych <- read.csv( here("helpers","psychs.csv"), sep = ";") # psychologic variables
subco <- read.csv( here("helpers","subcortical.csv"), sep = ",") # subcortical structures
hippo <- read.csv( here("helpers","hippocampus.csv"), sep = ",") # hippocampal structures


# UTILS ----

rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) ) # print rounded number
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) ) # delete leading zero
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) ) # calculate and print mean and SD

# fit regressions
fit_reg <- function(d, outcomes, X = "SUBJ * AHI.F + AGE + GENDER + SBTIV", w = F) {
  
  lapply(
    
    setNames(outcomes, outcomes),
    function(y) {
      
      if (w == T) lm(formula = as.formula( paste0(y," ~ ",X) ), data = d, weights = weights)
      else lm(formula = as.formula( paste0(y," ~ ",X) ), data = d, weights = NULL)
      
    }
  )
  
}

# extract linear regressions model diagnostics
lm_dia <- function(fit, type = "frequentist") {
  
  if (type == "frequentist") {
    
    sapply(
      
      names(fit),
      function(y)
        data.frame(
          X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ),
          p_breusch_pagan = c( check_heteroscedasticity(fit[[y]]) ),
          n_cook = sum( check_outliers(fit[[y]]), na.rm = T ),
          p_shapiro_wilk = c( check_normality(fit[[y]]) )
        )
    ) %>%
      
      t() %>%
      as.data.frame() %>%
      mutate( across( everything(), ~ unlist(.x, use.names = F) ) ) %>%
      mutate( heteroscedasticity = ifelse( p_breusch_pagan < .05, "!", ""), .after = p_breusch_pagan ) %>%
      mutate( nonnormality = ifelse( p_shapiro_wilk < .05, "!", ""), .after = p_shapiro_wilk ) %>%
      rownames_to_column("y")
    
  } else if (type == "Bayesian") sapply(
    
    names(fit),
    function(y)
      data.frame(
        X = sub( ".* ~ ", "", as.character( formula(fit[[y]]) )[1] ),
        sigma = if_else( grepl( "sigma", formula(fit[[y]])[2] ), sub( ")", "", sub( ".* ~ ", "", as.character( formula(fit[[y]]) )[2] ) ), "1" ),
        p_breusch_pagan = c( check_heteroscedasticity(fit[[y]]) ),
        n_cook = sum( check_outliers(fit[[y]], method = "mahalanobis_robust"), na.rm = T ),
        n_pareto = sum( check_outliers(fit[[y]], method = "pareto"), na.rm = T )
      )
  ) %>%
    
    t() %>%
    as.data.frame() %>%
    mutate( across( everything(), ~ unlist(.x, use.names = F) ) ) %>%
    mutate( heteroscedasticity = ifelse( p_breusch_pagan < .05, "!", ""), .after = p_breusch_pagan ) %>%
    `rownames<-`( 1:nrow(.) )
  
}

# Benjamini-Hochberg adjustment for 5% FDR
bh_adjust <- function(p) {
  
  # extract threshold
  bh_thres <- data.frame(
    p = sort(p), # sort p values from smallest to largest
    thres = .05 * ( 1:length(p) ) / length(p) # prepare BH thresholds for each p value
  ) %>%
    mutate( sig = if_else( p <= thres, T, F ) ) %>%
    filter( sig == T ) %>%
    select(thres) %>%
    max()
  
  # return stars based on this threshold
  return( if_else(p < bh_thres, "*", "") )
}

# extract coefficients for interactions and with them associated p-values (frequentist)
# (equivalent to ANOVAs with type 3 sum of squares which are appropriate for interactions)
# or posterior summaries (Bayesian)
lm_coeff <- function(fit, term = "SUBJ1:AHI.F1", type = "frequentist") {
  
  if (type == "frequentist") {
    
    sapply(
      
      names(fit),
      function(y)
        
        summary(fit[[y]])$coefficients[term, ] %>%
        t() %>%
        as.data.frame() %>%
        rename("p value" = "Pr(>|t|)") %>%
        mutate( `s value` = -log(`p value`, base = 2) ) %>%
        cbind( t( confint(fit[[y]])[term, ] ) ) %>%
        relocate(`2.5 %`, .before = `t value`) %>%
        relocate(`97.5 %`, .before = `t value`) %>%
        mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1)
      
    ) %>%
      
      t() %>%
      as.data.frame() %>%
      mutate(coefficient = term, .after = X) %>%
      mutate_if(is.list, unlist) %>%
      rownames_to_column("y") %>%
      mutate(
        sig_PCER = if_else(`p value` < .05, "*", ""),
        sig_FDR = bh_adjust(`p value`),
        sig_FWER = if_else(`p value` < .05/nrow(.), "*", "")
      )
    
  } else if (type == "Bayesian") sapply(
    
    names(fit),
    function(y)
      
      fixef(fit[[y]])[term, ] %>%
      t() %>%
      as.data.frame() %>%
      mutate(
        X = sub( ".* ~ ", "", as.character( formula(fit[[y]]) )[1] ),
        sigma = if_else( grepl( "sigma", formula(fit[[y]])[2] ), sub( ")", "", sub( ".* ~ ", "", as.character( formula(fit[[y]]) )[2] ) ), "1" ),
        .before = 1            
      )
    
  ) %>%
    
    t() %>%
    as.data.frame() %>%
    mutate(coefficient = term, .after = sigma) %>%
    mutate_if(is.list, unlist) %>%
    rownames_to_column("y")
  
}

# extract average per diagnosis slopes and interaction estimates using marginaleffects
meff <- function(fit, fit0, type = "moderation") {
  
  lapply(
    
    names(fit),
    function(y) {
      
      if (type == "moderation") {
        
        full_join(
          
          avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", wts = "weights", vcov = ~subclass) %>%
            as.data.frame() %>%
            mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
          
          avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", wts = "weights", vcov = ~subclass, hypothesis = "revpairwise") %>%
            as.data.frame() %>%
            mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 )

        ) %>%
          
          as.data.frame() %>%
          select( -starts_with("predicted") ) %>%
          mutate(y = y, .before = 1)
        
      } else if (type == "full") {
        
        reduce(
          list(
            avg_comparisons(fit0[[y]], variables = "SUBJ", wts = "weights", vcov = ~subclass) %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit0[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", wts = "weights", vcov = ~subclass) %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", wts = "weights", vcov = ~subclass) %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", wts = "weights", vcov = ~subclass, hypothesis = "revpairwise") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 )
          ),
          full_join
        ) %>%
          as.data.frame() %>%
          select( -starts_with("predicted") ) %>%
          mutate(y = y, .before = 1)
        
        
      }
    }
    
  ) %>% do.call( rbind.data.frame, . )

} 


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
    across( c("EDU.Y","AGE"), ~ as.numeric( scale(.x) ) )
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
  } ) %>%

  # scale neuropsychology variables (keeping MoCA non-transformed and will use binomial likelihood instead)
  mutate(
    across( all_of( subset(psych, domain == "Memory")$variable ), ~ as.numeric( scale(.x) ) ),
    across( all_of( subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable ), ~ -as.numeric( scale( log(.x) ) ) ),
    moca_max = 30
  ) %>%
  
  # scale total score of MDS-UPDRS III
  mutate( mds_updrs_iii = as.numeric( scale(mds_updrs_iii) ) )

# check for missing values
sapply( 1:nrow(df),  function(i) if ( sum( is.na(df[i, ]) ) > 0 ) df[i, ] ) %>%
  do.call( rbind.data.frame , . ) %>%
  mutate_if( is.numeric, round, 2 )


## ---- Propensity score matching -----

# fit a logistic regression based propensity score model based on sex and age
m.out <- matchit(
  formula = PD ~ AGE + GENDER,
  data = df,
  method = "full",
  distance = "glm",
  link = "logit"
)

# matching looks to have worked well with these data and variables
plot(m.out, type = "density", interactive = F)

# prepare a data set with propensity scores weights
df.matched <- match.data(m.out)


# RQ1: DOES OSA AFFECT SUBCORTICAL STRUCTURES VOLUME DIFFERENTLY IN PD PATIENTS COMPARED TO HEALTHY CONTROLS? ----
# RQ2: DOES OSA AFFECT HIPPOCAMPAL STRUCTURES VOLUME DIFFERENTLY IN PD PATIENTS COMPARED TO HEALTHY CONTROLS? ----

## Interaction boxplots ----

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


## ---- Linear regressions ----

# loop through structures and types of analysis
lapply(
  
  c("subco","hippo"),
  function(i) {

    with(

      get(i), {
        
        ### ---- classical regressions with threshold-based decisions ----
        
        fit0 <- fit_reg(df, name, X = "SUBJ * AGE + GENDER + SBTIV", w = F) # main effect of disease
        fit1 <- fit_reg(df, name, X = "SUBJ * AHI.F * AGE + GENDER + SBTIV", w = F) # main effect of OSA & disease/OSA interaction
        
        # do full analysis for subcortical structures and interaction only for hippocampal substructures
        if (i == "subco") write.table(
          
          # extract and write main effects and interactions for subcortical structures
          x = left_join(

            rbind.data.frame( lm_coeff(fit0, term = "SUBJ1"), lm_coeff(fit1, term = "AHI.F1"), lm_coeff(fit1, term = "SUBJ1:AHI.F1") ),
            rbind.data.frame( lm_dia(fit0), lm_dia(fit1) ),
            by = c("y","X")

          ) %>% mutate( sig_FDR = bh_adjust(`p value`) ), # re-calculate Benjamini-Hochberg adjusted significance statements
  
          file = here( "tables", paste0(i,"_classical_regressions.csv") ),
          sep = ",", row.names = F, quote = F

        ) else if (i == "hippo") write.table(
          
          # extract and write 'main effects and interactions only for hippocampal substructures
          x = left_join( lm_coeff(fit1, term = "SUBJ1:AHI.F1") ,lm_dia(fit1), by = c("y","X") ),
          file = here( "tables", paste0(i,"_classical_regressions.csv") ),
          sep = ",", row.names = F, quote = F
        )
        
        
        ### ---- weighted regressions with g-computation ----
        
        fit0 <- fit_reg(df.matched, name, X = "SUBJ * AGE + GENDER + SBTIV", w = T) # main effect of disease
        fit1 <- fit_reg(df.matched, name, X = "SUBJ * AHI.F * AGE + GENDER + SBTIV", w = T) # main effect of OSA & disease/OSA interaction
        
        # weighted regressions with g-computation
        write.table(
          
          # g-computation results
          x = left_join(
            meff( fit = fit1, fit0 = fit0, type = ifelse(i == "subco", "full", "moderation") ), # marginalised effects
            rbind.data.frame( lm_dia(fit0), lm_dia(fit1) ), # diagnostics
            by = c("y","X")
          ),

          file = here( "tables", paste0(i,"_weighted_regressions.csv") ),
          sep = ",", row.names = F, quote = F
        )
        
      }
    )
    
  }
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


## MDS-UPDRS III ----

# read the data
d1 <-
  read.csv( here("_data","mds_updrs_iii.csv"), sep = "," ) %>%
  left_join( d0[ , c("Study.ID","AGE","GENDER") ], by = "Study.ID" ) %>%
  mutate(
    AHI.F = as.factor(AHI.F),
    RESP = as.ordered(response_num),
    SEX = factor( case_when(GENDER == 1 ~ "male", GENDER == 0 ~ "female") ),
    ITEM = unlist(sapply( 1:nrow(.), function(i) paste( strsplit(item[i], ".", fixed = T)[[1]][1:2], collapse = "." ) ) ),
    NAME = unlist(sapply( 1:nrow(.), function(i) paste( strsplit(item[i], ".", fixed = T)[[1]][-c(1:2)], collapse = "." ) ) )
  ) %>%
  within( . , contrasts(AHI.F) <- -contr.sum(2)/2 )
  

### ---- item response modelling ----

# prepare formulas
formula1 <- list(
  
  `1PL_fixed` = bf( RESP ~ 1 + AHI.F +  (1 + AHI.F | ITEM) + (1 | Study.ID) ),
  `1PL_flex` = bf( RESP | thres(gr = ITEM) ~ 1 + AHI.F + AHI.F:ITEM + (1 | Study.ID) )
  
)

# re-fit the models for TMT-B and GPT via brms with differing variance for the groups
prior1 <- list(
  
  `1PL_fixed` =
    prior("normal(0, 3)", class = "Intercept") +
    prior("normal(0, 3)", class = "b") +
    prior("normal(0, 3)", class = "sd", group = "ITEM") +
    prior("normal(0, 3)", class = "sd", group = "Study.ID") +
    prior("lkj(1)", class = "cor"),

  `1PL_flex` =
    prior("normal(0, 3)", class = "Intercept") +
    prior("normal(0, 3)", class = "b") +
    prior("normal(0, 3)", class = "sd", group = "Study.ID")
  
)

# fit the models
fit1 <- lapply(
  
  setNames( names(prior1), names(prior1) ),
  function(k)
    
    brm( formula = formula1[[k]],
         data = d1,
         family = cumulative(link = "logit"),
         prior = prior1[[k]],
         seed = 87542
         )
  
)

# compare models' expected predictive performance
with( fit1, loo(`1PL_fixed`, `1PL_flex`) ) # selecting the fixed model for further description of patterns in the data (based on the model)


#### ---- plots ----

# some posterior predictive checks
lapply(
  
  c("ITEM","Study.ID"),
  function(k) {
    
    pp_check(fit1$`1PL_fixed`, type = "bars_grouped", ndraws = NULL, group = k) +
      theme(legend.position = "bottom")
    
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("motor_ppc_",tolower(k),".jpg") ),
      dpi = 300,
      width = ifelse(k == "ITEM", 10, 12),
      height = ifelse(k == "ITEM", 10, 12)
    )
    
  }
)

# plot item parameters
as_draws_df(fit1$`1PL_fixed`) %>%
  select( starts_with("r_ITEM") ) %>%
  pivot_longer(everything(), names_to = "coefficient", values_to = "value") %>%
  mutate(
    coefficient = sub("r_ITEM[", "", sub("]", "", coefficient, fixed = T), fixed = T),
    item = sub(",.*", "", coefficient),
    parameter = sub(".*,", "", coefficient)
  ) %>%
  group_by(item, parameter) %>%
  summarise(
    Md = median(value),
    Q2.5 = quantile(value, probs = .025),
    Q97.5 = quantile(value, probs = .975)
  ) %>%
  ungroup() %>%
  mutate( `Parameter: ` = factor(parameter, levels = c("Intercept","AHI.F1"), ordered = T) ) %>%
  
  ggplot() +
  aes(x = Md, xmin = Q2.5, xmax = Q97.5, y = reorder(item, Md, decreasing = T), colour = `Parameter: `) +
  geom_point( size = 4, position = position_dodge(width = -.75) ) +
  geom_pointrange( linewidth = 1.25, position = position_dodge(width = -.75) ) +
  scale_colour_manual( values = c("navyblue","red") ) +
  labs(
    y = NULL,
    x = "Parameter value",
    title = "Inverse difficulty item parameters of MDS-UPDRS III",
    subtitle = "The parameters are based on Item Response Theory logistic ordered model\nwith fixed item thresholds and hierarchical priors"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5)
  )

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","motor_item_parameters.jpg"),
  dpi = 300,
  width = 6,
  height = 11
)


# Show descriptively the profile PD-OSA+ vs PD-OSA- in MDS-UPDRS III
# Some type of forest plot with beta estimates
# BMI?

# dbs_combSTIM: re-test the same patients


## MoCA ----

# do it akin to this way
glm( cbind(moca, moca_max-moca) ~ SUBJ * AGE + EDU.Y + GENDER, data = df.matched, weights = weights, family = binomial(link = "logit") ) %>%
  avg_predictions(
    variables = "SUBJ",
    type = "response",
    vcov = ~subclass,
    transform = function(x) 30 * x,
    hypothesis = "revpairwise"
  )





# table for the classical analysis
fit_reg(df, subco$name, X = "SUBJ * AHI.F + AGE + GENDER + SBTIV", w = F) %>%
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

