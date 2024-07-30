# This script includes utility functions used throughout analyses for this project.
# Ought to be loaded before analysis starts.


# print rounded number ----
rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )

# delete leading zero ----
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) )

# calculate and print mean and SD ----
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )

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


# extract coefficients for interactions and with them associated p-values ----
# (equivalent to ANOVAs with type 3 sum of squares which are appropriate for interactions)
lm_coeff <- function(fit, term = "SUBJ1:AHI.F1") sapply(
  
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

# extract average per diagnosis slopes and interaction estimates using marginaleffects ----
meff <- function(fit, fit0, type = "moderation") {
  
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
            avg_comparisons(fit0[[y]], variables = "SUBJ") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit0[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 ),
            avg_comparisons(fit[[y]], variables = "AHI.F", by = "SUBJ", hypothesis = "revpairwise") %>% as.data.frame() %>% mutate( X = sub( paste0(y," ~ "), "", c( formula( fit[[y]] ) ) ), .before = 1 )
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
plot_ppc_stat <- function(fit, data, labs = psych, sleep = 5, stat = "sd", var = "event_group_osa") lapply(
  
  names(fit),
  function(i) {
    
    y <- paste0(i,"_sc")
    d <- subset( data, complete.cases( get(y) ) ) %>%
      mutate(
        x = case_when(
          var == "event_group_osa" ~ paste0(event,"_",group,"_OSA",osa),
          var == "event_group" ~ paste0(event,"_",group),
          var == "event_osa" ~ paste0(event,"_OSA",osa),
          var == "group_osa" ~ paste0(group,"_OSA",osa),
          var == "event" ~ event,
          var == "group" ~ group,
          var == "osa" ~ osa
        )
      )
    
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

# extract contrasts for marginalised effects calculation ----
contr_extr <- function(vars, data) lapply( setNames(vars,vars), function(x) contrasts( data[[x]] ) )

# extract contrasts of interest ----
contr_comp <- function(fit, y, contr, resc = NULL, negexp = F, summarise = T) {
  
  # extract posterior draws
  drws <-
    as_draws_df(fit) %>%
    select( starts_with("b") ) # keep group-level effects only
  
  # rescale if asked for
  if ( !is.null(resc) ) drws <- drws %>%
      
      mutate(
        b_Intercept = b_Intercept * resc[y,"SD"] + resc[y,"M"],
        across( -b_Intercept, ~ .x * resc[y,"SD"] )
      )
  
  ## compute marginal posteriors ----
  drws <- drws %>%
    
    mutate(
      
      ### ---- event ----
      
      Enrollment =
        b_Intercept +
        contr$event["enrollment", ] * b_event1,
      
      Retest =
        b_Intercept +
        contr$event["retest", ] * b_event1,
      
      ### ---- group ----
      
      HC =
        b_Intercept +
        contr$group["HC", ] * b_group1,
      
      PD =
        b_Intercept +
        contr$group["PD", ] * b_group1,
      
      ### ---- osa ----
      
      `OSA-` =
        b_Intercept +
        contr$osa["-", ] * b_osa1,
      
      `OSA+` =
        b_Intercept +
        contr$osa["+", ] * b_osa1,
      
      ### ---- event:group ----
      
      `Enrollment HC` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$event["enrollment", ] * contr$group["HC", ] * `b_event1:group1`,
      
      `Retest HC` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$event["retest", ] * contr$group["HC", ] * `b_event1:group1`,
      
      `Enrollment PD` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$event["enrollment", ] * contr$group["PD", ] * `b_event1:group1`,
      
      `Retest PD` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$event["retest", ] * contr$group["PD", ] * `b_event1:group1`,
      
      ### ---- event:osa ----
      
      `Enrollment OSA-` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["enrollment", ] * contr$osa["-", ] * `b_event1:osa1`,
      
      `Retest OSA-` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["retest", ] * contr$osa["-", ] * `b_event1:osa1`,
      
      `Enrollment OSA+` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["enrollment", ] * contr$osa["+", ] * `b_event1:osa1`,
      
      `Retest OSA+` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["retest", ] * contr$osa["+", ] * `b_event1:osa1`,
      
      ### ---- group:osa ----
      `HC OSA-` =
        b_Intercept +
        contr$group["HC", ] * b_group1 +
        contr$osa["-", ] * b_osa1 +
        contr$group["HC", ] * contr$osa["-", ] * `b_group1:osa1`,
      
      `HC OSA+` =
        b_Intercept +
        contr$group["HC", ] * b_group1 +
        contr$osa["+", ] * b_osa1 +
        contr$group["HC", ] * contr$osa["+", ] * `b_group1:osa1`,
      
      `PD OSA-` =
        b_Intercept +
        contr$group["PD", ] * b_group1 +
        contr$osa["-", ] * b_osa1 +
        contr$group["PD", ] * contr$osa["-", ] * `b_group1:osa1`,
      
      `PD OSA+` =
        b_Intercept +
        contr$group["PD", ] * b_group1 +
        contr$osa["+", ] * b_osa1 +
        contr$group["PD", ] * contr$osa["+", ] * `b_group1:osa1`,
      
      ### ---- event:group:osa ----
      
      #### ---- healthy controls ----
      
      `Enrollment HC OSA-` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["enrollment", ] * contr$group["HC", ] * `b_event1:group1` +
        contr$event["enrollment", ] * contr$osa["-", ] * `b_event1:osa1` +
        contr$group["HC", ] * contr$osa["-", ] * `b_group1:osa1` +
        contr$event["enrollment", ] * contr$group["HC", ] * contr$osa["-", ] * `b_event1:group1:osa1`,
      
      `Enrollment HC OSA+` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["enrollment", ] * contr$group["HC", ] * `b_event1:group1` +
        contr$event["enrollment", ] * contr$osa["+", ] * `b_event1:osa1` +
        contr$group["HC", ] * contr$osa["+", ] * `b_group1:osa1` +
        contr$event["enrollment", ] * contr$group["HC", ] * contr$osa["+", ] * `b_event1:group1:osa1`,
      
      `Retest HC OSA-` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["retest", ] * contr$group["HC", ] * `b_event1:group1` +
        contr$event["retest", ] * contr$osa["-", ] * `b_event1:osa1` +
        contr$group["HC", ] * contr$osa["-", ] * `b_group1:osa1` +
        contr$event["retest", ] * contr$group["HC", ] * contr$osa["-", ] * `b_event1:group1:osa1`,
      
      `Retest HC OSA+` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["HC", ] * b_group1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["retest", ] * contr$group["HC", ] * `b_event1:group1` +
        contr$event["retest", ] * contr$osa["+", ] * `b_event1:osa1` +
        contr$group["HC", ] * contr$osa["-", ] * `b_group1:osa1` +
        contr$event["retest", ] * contr$group["HC", ] * contr$osa["+", ] * `b_event1:group1:osa1`,
      
      #### ---- Parkinson's disease patients ----
      
      `Enrollment PD OSA-` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["enrollment", ] * contr$group["PD", ] * `b_event1:group1` +
        contr$event["enrollment", ] * contr$osa["-", ] * `b_event1:osa1` +
        contr$group["PD", ] * contr$osa["-", ] * `b_group1:osa1` +
        contr$event["enrollment", ] * contr$group["PD", ] * contr$osa["-", ] * `b_event1:group1:osa1`,
      
      `Enrollment PD OSA+` =
        b_Intercept +
        contr$event["enrollment", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["enrollment", ] * contr$group["PD", ] * `b_event1:group1` +
        contr$event["enrollment", ] * contr$osa["+", ] * `b_event1:osa1` +
        contr$group["PD", ] * contr$osa["+", ] * `b_group1:osa1` +
        contr$event["enrollment", ] * contr$group["PD", ] * contr$osa["+", ] * `b_event1:group1:osa1`,
      
      `Retest PD OSA-` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$osa["-", ] * b_osa1 +
        contr$event["retest", ] * contr$group["PD", ] * `b_event1:group1` +
        contr$event["retest", ] * contr$osa["-", ] * `b_event1:osa1` +
        contr$group["PD", ] * contr$osa["-", ] * `b_group1:osa1` +
        contr$event["retest", ] * contr$group["PD", ] * contr$osa["-", ] * `b_event1:group1:osa1`,
      
      `Retest PD OSA+` =
        b_Intercept +
        contr$event["retest", ] * b_event1 +
        contr$group["PD", ] * b_group1 +
        contr$osa["+", ] * b_osa1 +
        contr$event["retest", ] * contr$group["PD", ] * `b_event1:group1` +
        contr$event["retest", ] * contr$osa["+", ] * `b_event1:osa1` +
        contr$group["PD", ] * contr$osa["-", ] * `b_group1:osa1` +
        contr$event["retest", ] * contr$group["PD", ] * contr$osa["+", ] * `b_event1:group1:osa1`
      
    ) %>%
    
    # keep only the newly computed columns
    # all others (parameters) begin with "b_"
    select( !starts_with("b_") )
  
  # if reaction times need to exponentiated, do it
  if (negexp == T) drws <- drws %>% mutate( across( everything(), ~ exp(-.x) ) )
  
  # calculate differences
  drws <- drws %>%
    
    # event
    mutate(`Retest - Enrollment` = Retest - Enrollment, .after = Retest) %>%
    
    # group
    mutate(`PD - HC` = PD - HC, .after = PD) %>%
    
    # osa
    mutate(`OSA+ - OSA-` = `OSA+` - `OSA-`, .after = `OSA+`) %>%
      
    # event:group
    mutate(
      `Retest - Enrollment HC` = `Retest HC` - `Enrollment HC`,
      `Retest - Enrollment PD` = `Retest PD` - `Enrollment PD`,
      `Decline PD - Decline HC` = `Retest - Enrollment PD` - `Retest - Enrollment HC`,
      .after = `Retest PD`
    ) %>%
      
    # event:osa
    mutate(
      `Retest - Enrollment OSA-` = `Retest OSA-` - `Enrollment OSA-`,
      `Retest - Enrollment OSA+` = `Retest OSA+` - `Enrollment OSA+`,
      `Decline OSA+ - Decline OSA-` = `Retest - Enrollment OSA+` - `Retest - Enrollment OSA-`,
      .after = `Retest OSA+`
    ) %>%
    
    # group:osa
    mutate(
      `OSA+ - OSA- HC` = `HC OSA+` - `HC OSA-`,
      `OSA+ - OSA- PD` = `PD OSA+` - `PD OSA-`,
      `OSA difference PD - OSA difference HC` = `OSA+ - OSA- PD` - `OSA+ - OSA- HC`,
      .after = `PD OSA+`
    ) %>%
    
    # event:group:osa
    mutate(
      
      `Retest - Enrollment HC OSA-` = `Retest HC OSA-` - `Enrollment HC OSA-`,
      `Retest - Enrollment HC OSA+` = `Retest HC OSA+` - `Enrollment HC OSA+`,
      `Retest - Enrollment PD OSA-` = `Retest PD OSA-` - `Enrollment PD OSA-`,
      `Retest - Enrollment PD OSA+` = `Retest PD OSA+` - `Enrollment PD OSA+`,
      
      `Decline HC OSA+ - Decline HC OSA-` = `Retest - Enrollment HC OSA+` - `Retest - Enrollment HC OSA-`,
      `Decline PD OSA+ - Decline PD OSA-` = `Retest - Enrollment PD OSA+` - `Retest - Enrollment PD OSA-`,
      
      `Difference in differences PD - HC` = `Decline PD OSA+ - Decline PD OSA-` - `Decline HC OSA+ - Decline HC OSA-`,
      
      .after = `Retest PD OSA+`
    )
  
  # if to be summarised, do it
  if (summarise == T) drws <- drws %>%
    
    # summarise
    describe_posterior(test = "p_direction") %>%
    as.data.frame() %>%
    
    # prepare description columns
    mutate(
      Y = y,
      Term = c(
        rep("Event", 3), rep("Group", 3), rep("OSA", 3), # main effects
        rep("Event * Group",7), rep("Event * OSA",7), rep("Group * OSA",7), # two-way interactions,
        rep("Event * Group * OSA", 15) # three-way interaction
      ),
      Type = case_when(
        !grepl(" - ", Parameter) ~ "Estimate",
        grepl("Retest - Enrollment|PD - HC", Parameter) ~ "Simple main effect",
        grepl("OSA+ - OSA-", Parameter, fixed = T) ~ "Simple main effect",
        grepl("Decline|OSA difference", Parameter) ~ "Two-way interaction",
        grepl("Difference in differences", Parameter) ~ "Three-way interaction"
      ),
      .before = 1
    ) %>%
    
    # tidy it up
    mutate(
      Median = rprint(Median, 2),
      PPI = paste0("[",rprint(CI_low,2),", ",rprint(CI_high,2),"]"),
      pdir = paste0(rprint(100 * pd, 2), "%")
    ) %>%
    
    # keep only variables of interest
    select(Y, Term, Type, Parameter, Median, PPI, pdir)
  
  # return the transformed draws/their summaries
  return(drws)
  
}

