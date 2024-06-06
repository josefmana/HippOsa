# Computes regressions evaluating the effect of hippocampal volume on cognitive performance conditional on AHI and diagnosis.

rm( list = ls() ) # clear environment

# load packages
library(here)
library(openxlsx)
library(tidyverse)
library(ggdag)
library(patchwork)
library(brms)
#library(performance)

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("mods", "figs", "tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_bw(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# list all tests in the battery
# note: need to extract only those measures that are common across PD and CON, it ain't a level-II battery
psych <- list(
  memory = paste0( "avlt_", c("1_5","6","8","r_fp","r_fn") ),
  attention = c( "tmt_a", paste0("stroop_", c("body","slova") ) ),
  executive = c( "tmt_b", "stroop_barvy" ),
  speed = paste0( "gpt_", c("phk","lhk") )
)

# read the data set
d1 <- read.xlsx( here("_raw","TabHippAm1.xlsx") )
d2 <- read.csv( here("_raw","20221120_redcap_export.csv"), sep = "," )


# PRE-PROCESSING  ----

# prepare a data frame for analyses
df <- left_join(
  
  d1 %>% mutate( TIV = as.numeric( scale(`SBTIVmm^3`) ), across( contains("_"), ~ as.numeric( scale(.x) ) ) ),
  d2 %>% mutate( Study.ID = sub("-","",study_id) ) %>% select( Study.ID, all_of( unlist(psych, use.names = F) ) ),
  by = "Study.ID"
  
  
) %>%
  
  # pre-process predictor terms
  mutate_if( is.character, as.factor ) %>%
  mutate( GENDER = as.factor(GENDER), across( c("EDU.Y","AGE"), ~ as.numeric( scale(.x) ) ) ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(SUBJ) <- contr.sum(2)/2 # CON = 0.5, PD = -0.5
    contrasts(AHI.F) <- contr.sum(2)/2 # High = 0.5, Low = -0.5
    contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  } ) %>%
  
  # scale outcome and demographic variables
  mutate(
    across( "avlt_1_5", ~ as.numeric( scale(.x) ) ),
    across( all_of( unlist(psych[2:4], use.names = F) ), ~ as.numeric( scale(.x, center = F) ) )
  )

# check for missing values
sapply( 1:nrow(df),  function(i) if ( sum( is.na(df[i,]) ) > 0 ) df[i,] ) %>% do.call( rbind.data.frame , . )


# CAUSAL ASSUMPTIONS ----

# set-up coordinates for nodes
coords <- data.frame( name = c("cog","AHI","PD","Hippo","Age","Edu","Sex","TIV"), x = c(2,0,1,1,0,2,1,2), y = c(0,0,1,2,-2,-2,-2,2) )

# DAG with AHI ~ Hippocampus being confounded
dag1 <- dagify(
  
  cog ~ AHI + PD + Hippo + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  Hippo ~ PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Sex,
  AHI ~~ Hippo,
  coords = coords

)

# DAG with AHI causing Hippocampal volume
dag2 <- dagify(
  
  cog ~ AHI + PD + Hippo + Age + Edu + Sex + TIV,
  AHI ~ PD + Age + Sex,
  PD ~ Age + Sex,
  Hippo ~ AHI + PD + Age + TIV,
  Edu ~ Sex,
  TIV ~ PD + Sex,
  coords = coords
  
)

# set-up theme
theme_set( theme_dag() )

# set-up function for adjustment sets plotting
dag_plot <- function(dag, iv, dv, eff = "direct", leg = "none") ggdag_adjustment_set(dag, exposure = iv, outcome = dv, effect = eff, type = "minimal", shadow = F) + theme(legend.position = leg)

# plot it
( ggdag(dag1) | ggdag(dag2)  ) /
( dag_plot(dag1, "PD", "cog", "total", "none") |  dag_plot(dag2, "PD", "cog", "total", "none") ) /
( dag_plot(dag1, c("PD","AHI"), "cog", "total", "none") | dag_plot(dag2, c("PD","AHI"), "cog", "total", "none") ) /
( dag_plot(dag1, c("PD","AHI","Hippo"), "cog", "direct", "none") | dag_plot(dag2, c("PD","AHI","Hippo"), "cog", "direct", "none") ) +
  plot_annotation(tag_levels = "A")

# save it
ggsave( here("figs","hippo_dags.jpg"), dpi = 300, width = 10, height = 15 )
  

# STAT MODELS ----

# set-up a table with outcome information
outs <- data.frame(
  outcome = unlist(psych, use.names = F),
  likelihood = c( "gaussian", rep("binomial",4), rep("shifted_lognormal",7) ),
  tilde = c( " ~", rep(" | trials(15) ~",2), " | trials(35) ~", " | trials(15) ~", rep(" ~",7)  )
)

# set-up prediction terms
preds <- c(
  "SUBJ + AGE + GENDER + EDU.Y",
  "SUBJ * AHI.F + AGE + GENDER + EDU.Y",
  "SUBJ * AHI.F * R_Hipp_tail + SUBJ * AHI.F * L_Hipp_tail + TIV + AGE + GENDER + EDU.Y",
  "SUBJ * AHI.F * R_Hipp_body + SUBJ * AHI.F * L_Hipp_body + TIV + AGE + GENDER + EDU.Y",
  "SUBJ * AHI.F * R_Hipp_head + SUBJ * AHI.F * L_Hipp_head + TIV + AGE + GENDER + EDU.Y",
  "SUBJ * AHI.F * R_hippocampus + SUBJ * AHI.F * L_hippocampus + TIV + AGE + GENDER + EDU.Y",
  "SUBJ * AHI.F * R_amygdala + SUBJ * AHI.F * L_amygdala + TIV + AGE + GENDER + EDU.Y"
)

# fit a series of univariate regressions
fit <- lapply(
  
  setNames( 1:nrow(outs), outs$outcome ),
  function(i)
    
    lapply(
      
      setNames(preds,preds),
      function(j)
        
        brm(
          formula = bf( as.formula( paste0( outs$outcome[i], outs$tilde[i], " ", j ) ) ),
          data = df,
          family = outs$likelihood[i],
          prior = NULL
        )
      
    )
  
)

# fit a series of univariate regressions with neuropsychology outcomes and gradually increasing complexity
# of the linear predictor
# start by listing or predictor terms to use


# alternatively use predictor terms including adjustments for demographics
#preds <- c(" ~ 1 + Gender..0.F + Education + AGE",
#           " ~ 1 + Gender..0.F + Education + AGE + GROUP",
#           " ~ 1 + Gender..0.F + Education + AGE + GROUP * AHI_LH",
#           " ~ 1 + Gender..0.F + Education + AGE + GROUP * AHI_LH * HBT.R.BODY + GROUP * AHI_LH * HBT.L.BODY" 
#           )

# lastly, test for the AHI * PD * AGE interaction adjusting estimates for demographics (similarly to biopd_psychoAHI_pd_only.R)
#preds <- c(" ~ 1 + Gender..0.F + Education",
#           " ~ 1 + Gender..0.F + Education + AGE",
#           " ~ 1 + Gender..0.F + Education + AGE * GROUP",
#           " ~ 1 + Gender..0.F + Education + AGE * GROUP * AHI_LH" 
#           )

# prepare a list for formulas
f <- list()

# set-up the linear models
for ( i in names(psych) ) {
  for ( j in psych[[i]] ) {
    for ( k in preds ) f[[i]][[j]][[k]] <- paste0( j, k ) %>% as.formula()
  }
}


# ---- implementation: learning from the data  ----

# use the default lm() function to fit the regressions via the QR-decomposition and least squares
m <- list()

# loop through all linear models in f, select appropriate likelihoods and links (via family argument)
for ( i in names(f) ) {
  for ( j in names(f[[i]]) ) {
    for ( k in names(f[[i]][[j]]) ) m[[i]][[j]][[k]] <- lm( formula = f[[i]][[j]][[k]], data = df )
  }
}


# ---- post-processing: model summaries ----

# prepare a list for ANOVA results and model checks
t <- list()

# fill-in all the values
for ( i in names(m) ) {
  for ( j in names(m[[i]]) ) {
    
    # fill-in all values we care about at this point
    for ( k in 1:length(preds) ) t[[i]][[j]][[preds[[k]]]] <- c(
      
      # number of subjects per group
      ( !is.na( df[ df$GROUP == "CON", j ] ) ) %>% sum(), # number of controls
      ( !is.na( df[ df$GROUP == "PD", j ] ) ) %>% sum(), # number of patients
      
      # ANOVA results
      if( k == 1 ) rep(NA, length(preds) ) %>% as.numeric(),
      if( k != 1 ) anova( m[[i]][[j]][[k-1]] , m[[i]][[j]][[k]] )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric(),
      
      # empty column for a statistical significance after Benjamini-Hochberg correction
      NA,
      
      # Bayes Factor via the BIC approximation (https://doi.org/10.3758/BF03194105)
      if ( k == 1 ) NA %>% as.numeric(),
      if ( k != 1 ) exp( ( BIC(m[[i]][[j]][[k-1]]) - BIC(m[[i]][[j]][[k]]) ) / 2 ),
      
      # model checks
      check_normality( m[[i]][[j]][[k]] ) %>% as.numeric(), NA, # p ≤ .05 rejects normality
      check_heteroskedasticity( m[[i]][[j]][[k]] ) %>% as.numeric(), NA, # p ≤ .05 rejects homogeneity of variances
      check_outliers( m[[i]][[j]][[k]] ) %>% sum() # number of outliers according to the default setting in the "performance" package
      
    ) %>% `names<-`( c("nCON","nPD","F","df1","df2","Pr(>F)","sig.","BF10","normality (p-value)","non-normality","homoscedasticity (p-value)","heteroscedasticity","outliers") )
    
    # collapse tables for a single outcome to a neat table
    t[[i]][[j]] <- do.call( cbind.data.frame, t[[i]][[j]] ) %>%
      t() %>% # flip dimensions
      as.data.frame() %>% rownames_to_column( var = "predictors" ) # add column denoting model type
  }
  
  # collapse tables in a single domain to an even neater table
  t[[i]] <- do.call( rbind.data.frame, t[[i]] ) %>%
    rownames_to_column( var = "outcome" ) %>% # add column denoting model type
    mutate( outcome = sub( "\\..*", "", outcome) ) # tidy the outcome variable
  
}

# pull all model summaries into a one neat table
t <- do.call( rbind.data.frame, t ) %>%
  rownames_to_column( var = "domain" ) %>%
  mutate( domain = sub( "\\..*", "", domain) )

# extract the Benjamini-Hochberg corrected threshold
bh_thres <- data.frame( p = sort( t$`Pr(>F)`), # order the p-values from lowest to largest
                        thres = .05 * (1:nrow( t[complete.cases(t$`Pr(>F)`),] ) ) / nrow( t[complete.cases(t$`Pr(>F)`),] ) # prepare BH thresholds for each p-value
                        ) %>%
  # flag BH-significant p-values and extract the largest threshold (https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
  mutate( sig = ifelse( p <= thres, T, F ) ) %>% filter( sig == T ) %>% select(thres) %>% max()

# add final touches, i.e., statistical significance decisions for ANOVAs (based on BH correction),
# non-normality (based on nominal p < .05) and heteroscedasticity (based on nominal p < .05)
t <- t %>% mutate( `non-normality` = ifelse( `normality (p-value)` <= .05, "!", NA ),
                   heteroscedasticity = ifelse( `homoscedasticity (p-value)` <= .05, "!", NA ),
                   `sig.` = ifelse( `Pr(>F)` < bh_thres, ":-)", NA )
                   )

# print as a csv (use the second line if demographics were adjusted for)
write.table( t, "tables/models_comps_&_checks.csv", sep = ",", row.names = F, na = "" )
#write.table( t, "tables/demographic_adjusted_models_comps_&_checks.csv", sep = ",", row.names = F, na = "" )
#write.table( t, "tables/age_interaction_comps_&_checks.csv", sep = ",", row.names = F, na = "" )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/psychoANOVAs.txt" )
