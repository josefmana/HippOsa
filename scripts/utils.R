# This script includes utility functions used throughout analyses for this project.
# Ought to be loaded before analysis starts.


# print rounded number ----
rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )

# delete leading zero ----
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) )

# calculate and print mean and SD ----
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )

# get frequency and proportion of binary variables ----
freqprop <- function(x, d = 0) paste0( table(x)[2], " (", rprint( 100*prop.table( table(x) )[2], d = d ), "%)" )

# extract t/z values and p values ----
statextract <- function(coeffs, y, stat = "t") coeffs %>%
  
  as.data.frame() %>%
  mutate(y = y, .before = 1) %>%
  rownames_to_column("coefficient") %>%
  
  mutate(
    out = paste0(
      stat," = ", rprint( get( paste0(stat," value") ), 2 ),
      ", p ",
      if_else(
        get( paste0("Pr(>|",stat,"|)") ) < .001,
        zerolead( get( paste0("Pr(>|",stat,"|)") ) ),
        paste0("= ", zerolead( get( paste0("Pr(>|",stat,"|)") ) ) )
      )
    )
  ) %>%
  select( -paste0(stat," value"), -paste0("Pr(>|",stat,"|)") ) %>%
  pivot_wider(
    names_from = coefficient,
    values_from = out
  ) %>%
  select(-`(Intercept)`)


# fit regressions ----
fit_reg <- function(d, outcomes, X = "SUBJ * AHI.F + AGE + GENDER + SBTIV", w = F) {
  
  lapply(
    
    setNames(outcomes, outcomes),
    function(y) {
      
      if (w == T) lm(formula = as.formula( paste0(y," ~ ",X) ), data = d, weights = weights)
      else lm(formula = as.formula( paste0(y," ~ ",X) ), data = d, weights = NULL)
      
    }
  )
  
}

# extract linear regressions model diagnostics ----
lm_dia <- function(fit) sapply(
  
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
  

# Benjamini-Hochberg adjustment for 5% FDR ----
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


# extract coefficients ----
# for interactions and with them associated p-values (frequentist) (equivalent to ANOVAs with type 3 sum of squares which are appropriate for interactions)
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

# extract average per diagnosis slopes and interaction estimates using marginaleffects ----
mass <- function(fit, type = "moderation") {
  
  lapply(
    
    names(fit),
    function(y) {
      
      if (type == "moderation") {
        
        full_join(
          
          avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ") %>%
            as.data.frame() %>%
            mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
          
          avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", hypothesis = "revpairwise") %>%
            as.data.frame() %>%
            mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 )
          
        ) %>%
          
          as.data.frame() %>%
          select( -starts_with("predicted") ) %>%
          mutate(y = y, .before = 1)
        
      } else if (type == "full") {
        
        reduce(
          list(
            avg_comparisons(fit[[y]], variables = "SUBJ") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", hypothesis = "revpairwise") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 )
          ),
          full_join
        ) %>%
          as.data.frame() %>%
          select( -starts_with("predicted") ) %>%
          mutate(y = y, .before = 1)
        
        
      } else if (type == "fullest") {
        
        reduce(
          list(
            # 'main effects'
            avg_comparisons(fit[[y]], variables = "SUBJ") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            
            # 'simple main effects'
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "SUBJ", by = "AHI.F") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            
            # interactions
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", hypothesis = "revpairwise") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "SUBJ", by = "AHI.F", hypothesis = "pairwise") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 )
            
            
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

# plot posterior predictive stats for a set of models ----
plot_ppc_stat <- function(fit, data, labs = psych, sleep = 5, stat = "sd", var = "group_osa") lapply(
  
  names(fit),
  function(y) {
    
    d <- subset( data, complete.cases( get(y) ) ) %>%
      mutate(
        x = case_when(
          var == "group_osa" ~ paste0(SUBJ,"_",AHI.F),
          var == "group" ~ SUBJ,
          var == "osa" ~ AHI.F
        )
      )
    
    print(
      pp_check(
        object = d[ , y],
        yrep = posterior_predict(fit[[y]], newdata = d),
        fun = ppc_stat_grouped,
        stat = stat,
        group = d$x
      ) +
        labs(
          title = with(labs , label[variable == y] ),
          subtitle = paste0("Observed (thick line) vs predicted (histogram) ", stat)
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
