# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "tidyverse", "dplyr", # data wrangling
           "ggplot2", "patchwork", # plotting
           "performance" # model quality checking
           )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_minimal(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# list all tests in the battery
# note: need to extract only those measures that are common across PD and CON, it ain't a level-II battery
psych <- list( memory = paste0( "avlt_", c("1_5","6","8","r_fp","r_fn") ),
               attention = c( "tmt_a", paste0("stroop_", c("body","slova") ) ),
               executive = c( "tmt_b", "stroop_barvy" ),
               speed = paste0( "gpt_", c("phk","lhk") )
               )

# read the data set
d1 <- read.csv( "data/pd_con.csv" , sep = "," )
d2 <- read.csv( "data/20221120_redcap_export.csv", sep = "," )


# ---- pre-processing  ----

# keep only included subjects in d2
d2 <- d2 %>%
  mutate( Study.ID = sub( "-", "", study_id ) ) %>%
  slice( which( Study.ID %in% unique(d1$Study.ID) ) )

# prepare a data frame for analyses
df <- d1 %>% left_join( d2, by = "Study.ID" ) %>% # glue d2 (redcap export of neuropsychology) to d1 (Kris's data file)
  
  # select only those columns (variables) that we will need
  select( Study.ID, # subject's id
          which( names(.) %in% unlist(psych) ), # neuropsychological tests
          GROUP, AHI_LH, Gender..0.F, AGE, Education, # diagnosis (GROUP) and demographics
          contains( "HBT" ) & contains( "BODY" ) # hippocampus body volumes
          ) %>%
  
  # pre-process predictor terms
  mutate_if( is.character, as.factor ) %>%
  mutate_if( grepl( "HBT", names(.) ), function(x) scale(x) ) %>%
  mutate( Gender..0.F = as.factor( Gender..0.F), Education = scale( Education ), AGE = scale(AGE) ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(GROUP) <- contr.treatment(2) # CON = 0, PD = 1
    contrasts(AHI_LH) <- contr.treatment(2) # High = 0, Low = 1
    contrasts(Gender..0.F) <- contr.treatment(2) # female = 0, male = 1
  } ) %>%
  
  # log-transform reaction times and false positives/negatives
  mutate_if( names(.) %in% unlist( psych[2:4] ), function(x) log(x) ) %>%
  mutate_if( grepl( "avlt_r", names(.) ), function(x) log(x+1) )

# check for missing values
sapply( 1:nrow(df), # loop through all the rows (subjects)
        function(i) if ( sum( is.na(df[i,]) ) > 0 ) df[i,] # print the full entry of each subject with missing value
        ) %>%
  do.call( rbind.data.frame , . ) # put all rows with missing values together


# ---- representation: statistical model  ----

# fit a series of univariate regressions with neuropsychology outcomes and gradually increasing complexity
# of the linear predictor
f <- list()

# set-up the linear models
for ( i in names(psych) ) {
  for ( j in psych[[i]] ) {
    for ( k in c(" ~ 1", # 1) intercept only
                 " ~ 1 + GROUP", # 2) effect of diagnosis
                 " ~ 1 + GROUP * AHI_LH", # 3) added effect of AHI
                 " ~ 1 + GROUP * AHI_LH * HBT.R.BODY + GROUP * AHI_LH * HBT.L.BODY" ) # 4) mediation via hippocampal volume
    ) f[[i]][[j]][[k]] <- paste0( j, k ) %>% as.formula()
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
    
    # first fill-in all diagnostics we care about at this point
    for ( k in names(m[[i]][[j]]) ) t[[i]][[j]][[k]] <- c(
      
      # number of subjects per group
      ( !is.na( df[ df$GROUP == "CON", j ] ) ) %>% sum(), # number of controls
      ( !is.na( df[ df$GROUP == "PD", j ] ) ) %>% sum(), # number of patients
      
      # ANOVA results
      case_when(
        k == " ~ 1" ~ rep(NA, 4) %>% as.numeric(),
        k == " ~ 1 + GROUP" ~ with( m[[i]][[j]], anova( ` ~ 1`, ` ~ 1 + GROUP` ) )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric(),
        k == " ~ 1 + GROUP * AHI_LH" ~ with( m[[i]][[j]], anova( ` ~ 1 + GROUP`, ` ~ 1 + GROUP * AHI_LH` ) )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric(),
        k == " ~ 1 + GROUP * AHI_LH * HBT.R.BODY + GROUP * AHI_LH * HBT.L.BODY" ~ with( m[[i]][[j]], anova( ` ~ 1 + GROUP * AHI_LH`, ` ~ 1 + GROUP * AHI_LH * HBT.R.BODY + GROUP * AHI_LH * HBT.L.BODY` ) )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric()
      ),
      
      # empty column for a statistical significance after Benjamini-Hochberg correction
      NA,
      
      # model checks
      check_normality( m[[i]][[j]][[k]] ) %>% as.numeric(), NA, # p ≤ .05 rejects normality
      check_heteroskedasticity( m[[i]][[j]][[k]] ) %>% as.numeric(), NA, # p ≤ .05 rejects homogeneity of variances
      check_outliers( m[[i]][[j]][[k]] ) %>% sum() # number of outliers according to the default setting in the "performance" package
      
    ) %>% `names<-`( c("nCON","nPD","F","df1","df2","Pr(>F)","sig.","normality (p-value)","non-normality","homoscedasticity (p-value)","heteroscedasticity","outliers") )
    
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
  # flag BH-significant p-values and extract the largest threshold as per https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
  mutate( sig = ifelse( p <= thres, T, F ) ) %>% filter( sig == T ) %>% select(thres) %>% max()

# add final touches, i.e., statistical significance decisions for ANOVAs (based on BH correction),
# non-normality (based on nominal p < .05) and heteroscedasticity (based on nominal p < .05)
t <- t %>% mutate( `non-normality` = ifelse( `normality (p-value)` <= .05, "!", NA ),
                   heteroscedasticity = ifelse( `homoscedasticity (p-value)` <= .05, "!", NA ),
                   `sig.` = ifelse( `Pr(>F)` < bh_thres, ":-)", NA )
                   )

# print as a csv
write.table( t, "tables/models_comps_&_checks.csv", sep = ",", row.names = F, na = "" )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/psychoANOVAs.txt" )
