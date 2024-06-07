# Computes regressions evaluating the effect of hippocampal volume on cognitive performance conditional on AHI and diagnosis.

rm( list = ls() ) # clear environment

# load packages
library(here)
library(openxlsx)
library(tidyverse)
library(ggdag)
library(patchwork)
library(car) # for easy to compute Type II and Type III sum of squares
library(performance)

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("mods", "figs", "tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )

# sensible contrasts for Type III Anovas
options( contrasts = c("contr.sum","contr.poly") )

# list all tests in the battery
# note: need to extract only those measures that are common across PD and CON, it ain't a level-II battery
psych <- read.csv("psychs.csv", sep = ";")
preds <- read.csv("brains.csv", sep = ";")

# read the data set
d1 <- read.xlsx( here("_raw","TabHippAm1.xlsx") )
d2 <- read.csv( here("_raw","20221120_redcap_export.csv"), sep = "," )


# PRE-PROCESSING  ----

# prepare a data frame for analyses
d0 <- left_join( d1, d2 %>% mutate( Study.ID = sub("-","",study_id) ) %>% select( Study.ID, all_of(psych$variable) ), by = "Study.ID" )

# format it for analyses
df <- d0 %>%
  
  # pre-process predictor terms
  mutate_if( is.character, as.factor ) %>%
  mutate( TIV = as.numeric( scale(`SBTIVmm^3`) ), across( contains("L_") | contains("R_"), ~ as.numeric( scale(.x) ) ) ) %>%
  mutate( GENDER = as.factor(GENDER), across( c("EDU.Y","AGE"), ~ as.numeric( scale(.x) ) ) ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  #within( . , {
  #  contrasts(SUBJ) <- contr.sum(2)/2 # CON = 0.5, PD = -0.5
  #  contrasts(AHI.F) <- contr.sum(2)/2 # High = 0.5, Low = -0.5
  #  contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
  #} ) %>%

  # scale outcome variables
  mutate(
#    across( "avlt_1_5", ~ as.numeric( scale(.x) ) ),
    across( all_of( subset(psych, domain == "Memory")$variable ), ~ as.numeric( scale(.x) ) ),
    across( all_of( subset(psych, domain != "Memory")$variable ), ~ as.numeric( scale( log(.x) ) ) )
  )

# check for missing values
sapply( 1:nrow(df),  function(i) if ( sum( is.na(df[i, ]) ) > 0 ) df[i, ] ) %>%
  do.call( rbind.data.frame , . ) %>%
  mutate_if( is.numeric, round, 2 )


# CAUSAL ASSUMPTIONS ----

# set-up coordinates for nodes
coords <- data.frame(

  name = c("cog","AHI","PD","Hippo","Age","Edu","Sex","TIV"),
  x = c(2,0,1,1,0,2,1,2),
  y = c(0,0,1,2,-2,-2,-2,2)

)

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

# fit a series of univariate regressions
fit <- lapply(
  
  setNames( 1:nrow(psych), psych$label ),
  function(i)
    
    lapply(
      
      setNames( 1:nrow(preds),preds$predictor ),
      function(j) {
        
        form <- paste0( psych$LHS[i], preds$RHS[j] )
        #print( paste0("outcome: ",psych$label[i],", predictor: ",preds$label[j],", likelihood: ",psych$likelihood[i],", formula: ",form) )
        
        #if( psych$likelihood[i] == "gaussian" ) mod <-  lm( formula = as.formula(form), data = df )
        #else if( psych$likelihood[i] == "binomial" ) mod <- glm( formula = as.formula(form), data = df, family = binomial() )
        return( lm( formula = as.formula(form), data = df ) )
        
      }
    )
  
)


# MODEL SUMMARIES ----

# set theme for plotting
theme_set( theme_bw(base_size = 14) )

# extract statistical tests results
tab2 <- lapply(
  
  with( preds, setNames(predictor,predictor) ),
  function(i)
    
    sapply(
      
      names(fit),
      function(j) {
        
        print( paste0("outcome: ",j,", predictor: ",i) )
        
        Anova(fit[[j]][[i]], type = 3)[i, ] %>%
          as.data.frame() %>%
          #rename_with( ~ "p value", starts_with("Pr(>") ) %>%
          #rename_with( ~ "Test statistic", any_of( c("F value","LR Chisq") ) ) %>%
          #select(`Test statistic`, Df, `p value`) %>%
          select( starts_with("F value"), Df, starts_with("Pr(>") ) %>%
          return()
        
      }
    ) %>%
    
    t() %>%
    as.data.frame() %>%
    mutate( across( everything(), ~ unlist(.x, use.names = F) ) ) %>%
    rownames_to_column("Outcome") %>%
    mutate( Predictor = preds[ preds$predictor == i, "label2" ], .before = 1 ) %>%
    mutate( `Domain: ` = sapply( 1:nrow(.), function(j) psych[ psych$label == Outcome[j], "domain"] ) )
  
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  mutate( Predictor = factor(Predictor, levels = rev(preds$label2), ordered = T) ) %>%
  mutate( Outcome = factor(Outcome, levels = psych$label, ordered = T) ) %>%
  mutate( `Domain: ` = factor(`Domain: `, levels = unique(psych$domain), ordered = T) )

# write as csv
write.table( x = tab2, file = here("tabs","hippo_anovas.csv"), sep = ",", row.names = F, quote = F )

# plot p-values
tab2 %>%
  
  ggplot() +
  aes(x = `Pr(>F)`, y = Predictor, fill = `Domain: `) +
  geom_bar(stat = "identity") +
  labs( x = "p-value", y = NULL ) +
  geom_vline( xintercept = .05, linetype = "solid", linewidth = 1, colour = "black" ) +
  facet_wrap( ~ Outcome ) +
  scale_fill_manual( values = c("#F9A729","#64CDCC","#A4A4D5","#5FBB68") ) +
  theme(legend.position = "bottom")

# save it
ggsave( plot = last_plot(), file = here("figs","hippo_pvalues.jpg"), dpi = 300, width = 12.6, height = 13.3 )

# plot interactions per outcome
fig2 <-
  
  lapply(
    
    with( psych, setNames(variable,label) ),
    function(i)
      
      d0 %>%
      select( SUBJ, AHI.F, all_of(i), all_of(sub("SUBJ:AHI.F:","c",preds$predictor)[-c(1:2)]) ) %>%
      pivot_longer( cols = starts_with("c"), names_to = "Structure", values_to = "Volume" ) %>%
      
      mutate( `Group: ` = case_when(SUBJ == "CON" ~ "HC", SUBJ == "PD" ~ "PD") ) %>%
      mutate( AHI = factor( case_when(AHI.F == "L" ~ "low AHI", AHI.F == "H" ~ "high AHI"), levels = paste0( c("low","high"), " AHI"), ordered = T ) ) %>%
      
      mutate( Structure = sub("c","SUBJ:AHI.F:",Structure) ) %>%
      mutate( Structure = sapply( 1:nrow(.), function(i) preds[ preds$predictor == Structure[i], "label2" ], USE.NAMES = F ) ) %>%
      mutate( Structure = factor( Structure, levels = preds$label2[-c(1:2)], ordered = T ) ) %>%
      
      ggplot() +
      aes(x = Volume, y = get(i), colour = `Group: `, fill = `Group: `) +
      geom_point(size = 3) +
      labs( x = bquote("Standardized volume"~("mm"^3) ), y = with( psych, label2[variable == i] ) ) +
      geom_smooth(method = "lm", linewidth = 1.5, alpha = .25) +
      ggh4x::facet_grid2( Structure ~ AHI, scales = "free_x", independent = "x" ) +
      theme_bw( base_size = 16 ) +
      theme(legend.position = "bottom")
    
  )

# save it
for( i in names(fig2) ) ggsave( plot = fig2[[i]], file = here( "figs", paste0(i,"_intplot.jpg") ), dpi = 300, width = 12.6, height = 26.6 )




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
