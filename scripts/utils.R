# This script includes utility functions used throughout analyses for this project.
# Ought to be loaded before analysis starts.


# print rounded number
rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )

# delete leading zero
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) )

# calculate and print mean and SD
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )

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

# extract contrasts for marginalised effects calculation
contr_extr <- function(vars, data) lapply( setNames(vars,vars), function(x) contrasts( data[[x]] ) )

# extract contrasts of interest
contr_comp <- function(fit, var, contr, resc = NULL, negexp = F, summarise = T) {
  
  # extract posterior draws
  drws <-
    as_draws_df(fit) %>%
    select( starts_with("b") ) # keep group-level effects only
  
  # rescale if asked for
  if ( !is.null(resc) ) drws <- drws %>%
      
      mutate(
        b_Intercept = b_Intercept * resc[var,"SD"] + resc[var,"M"],
        across( -b_Intercept, ~ .x * resc[var,"SD"] )
      )
  
  # compute marginal posteriors
  drws <- drws %>%
    
    mutate(
      
      # event
      
      # group
      
      # osa
      
      # event:group
      `enrollment HC` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$event["enrollment", ] * contr$group["HC", ] * `b_event1:group1`,
      
      `retest HC` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$event["retest", ] * contr$group["HC", ] * `b_event1:group1`,
      
      `enrollment PD` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$event["enrollment", ] * contr$group["PD", ] * `b_event1:group1`,
      
      `retest PD` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$event["retest", ] * contr$group["PD", ] * `b_event1:group1`,
      
      # event:osa
      `enrollment OSA-` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["enrollment", ] * contr$osa["-", ] * `b_event1:osa1`,
      
      `retest OSA-` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["retest", ] * contr$osa["-", ] * `b_event1:osa1`,
      
      `enrollment OSA+` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["enrollment", ] * contr$osa["+", ] * `b_event1:osa1`,
      
      `retest OSA+` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["retest", ] * contr$osa["+", ] * `b_event1:osa1`,
      
      # group:osa
      
      # event:group:osa
      
      
    ) %>%
    
    # keep only the newly computed columns
    # they are the only one with spces
    select( contains(" ") )
  
  # if reaction times need to exponentiated, do it
  if (negexp == T) drws <- drws %>%
    
    mutate( across( everything(), ~ -exp(.x) ) )
  
  # calculate differences
  drws <- drws %>%
    
    mutate(
      
      # event:group
      `retest - enrollment HC` = `retest HC` - `enrollment HC`,
      `retest - enrollment PD` = `retest PD` - `enrollment PD`,
      `decline PD - decline HC` = `retest - enrollment PD` - `retest - enrollment HC`
      
      # 
      
    )
  
  # if to be summarised, do it
  if (summarise == T) drws <- drws %>% describe_posterior(test = "p_direction")
  
}

