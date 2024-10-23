# Computes regressions evaluating the association between OSA and brain structures' volume/cognition:
#
# RQ1) the distribution of subcortical brain structures' volume in de novo PD conditional on OSA
# RQ2) the distribution of  hippocampal substructures' volume in de novo PD conditional on OSA
# RQ3) the distribution of cognitive performance in de novo PD conditional on OSA


# LINEAR REGRESSIONS ----
fit_regressions <- function(df, help){
  
  # extract single helpers
  for (i in names(help) ) assign(i, help[[i]])
  
  # set-up formulas
  forms <- data.frame(
    
    object = c("subco","hippo","psych"),
    y = c("subcortical", "hippocampi","cognition"),
    X = c( rep("SUBJ * AHI.F + AGE + GENDER + SBTIV", 2), "SUBJ * AHI.F + AGE + GENDER + EDU.Y")
    
  )
  
  # loop through types of regressions (base vs interaction), and structures
  lapply(
    
    set_names(x = 1:nrow(forms), nm = forms$y),
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
        if (y == "hippocampi") return(
          
          left_join( lm_coeff(fit, term = "SUBJ1:AHI.F1") ,lm_dia(fit), by = c("y","X") )
          
        ) else  return(
          
          left_join(
            
            rbind.data.frame(
              lm_coeff(fit, term = "SUBJ1"),
              lm_coeff(fit, term = "AHI.F1"),
              lm_coeff(fit, term = "SUBJ1:AHI.F1")
            ),
            lm_dia(fit),
            by = c("y","X")
            
          ) %>% mutate( sig_FDR = bh_adjust(`p value`) ) # re-calculate Benjamini-Hochberg adjusted significance statements

        )
      }
    )
  ) %>% return()

}

# SET-UP FORMULAS FOR HETEROSCEDASTIC MODELS ----
set_formulas <- function()  list(
  
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

# BAYESIAN REGRESSIONS WITH HETEROSCEDASTICITY ----
fit_bayesian <- function(df, formulas) {
  
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
        function(j) brm(
          
          formula = formulas[[i]][[j]],
          data = df,
          family = gaussian(link = "identity"),
          prior = priors[[i]],
          seed = 87542
          
        )
      )
  )
  
  # return models
  return(fit)
  
}


## ---- POSTERIOR PREDICTIVE CHECKS ----
model_check <- function(fit, help, formulas, which = "prior_sense") {
  
  if ( which == "prior_sense") return(
    
    # check prior-sensitivity
    lapply(
      set_names( names(fit) ),
      function(i)
        lapply(
          set_names( names(fit[[i]]) ),
          function(y)
            powerscale_sensitivity(fit[[i]][[y]])
        )
    )
    
  ) else if (which == "ppc") return(
    
    lapply(
      
      set_names( c("dens_overlay_grouped","stat_grouped") ),
      function(k) {
        
        ppc <- lapply(
          
          set_names( names(fit) ),
          function(i)
            
            lapply(
              
              set_names( names(fit[[i]]) ),
              function(y) {
                
                # set distinct colour palletes for different types of models
                if (i == "varequal") bayesplot::color_scheme_set( scheme = "red" )
                else if (i == "heteroscedastic") bayesplot::color_scheme_set( scheme = "blue" )
                
                set.seed(87542) # seed for reproducibility
                
                # plotting proper
                pp_check(fit[[i]][[y]], ndraws = 100, type = k, group = "SUBJ", stat = "sd", bins = 12) +
                  theme_bw(base_size = 14) +
                  labs(
                    x = paste0( with(help$psych, label[variable == y]), " (-log seconds)" ),
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
        return(
          with(
            ppc,
            ( varequal$tmt_b | heteroscedastic$tmt_b ) /
              ( varequal$gpt_phk | heteroscedastic$gpt_phk ) /
              ( varequal$gpt_lhk | heteroscedastic$gpt_lhk )
          )
          
        )
        
        # save it
        #ggsave(
        #  plot = last_plot(),
        #  filename = here( "figures", paste0("cognition_ppc_", sub("_.*","",k), ".jpg") ),
        #  dpi = 300,
        #  width = 15,
        #  height = 15
        #)
        
      }
    )
    
  )
  
}


## ---- COEFFICIENTS EXTRACTION ----
extract_coefficients <- function(fit) lapply(
  
  names(fit),
  function(i) do.call(
    
    rbind.data.frame,
    list(
      lm_coeff(fit[[i]], term = "SUBJ1", type = "Bayesian"),
      lm_coeff(fit[[i]], term = "AHI.F1", type = "Bayesian"),
      lm_coeff(fit[[i]], term = "SUBJ1:AHI.F1", type = "Bayesian")
    )
    
  )
  
)  %>% do.call( rbind.data.frame, . )
