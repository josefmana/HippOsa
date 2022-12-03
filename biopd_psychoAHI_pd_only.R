# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "tidyverse", "dplyr", # data wrangling
           "ggplot2", "patchwork", # plotting
           "ggdag", "dagitty", # causal models
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
# note: extracting only those measures that are present in most patients (some older don't have the same battery),
# it's the set of tests common across Bio-PD and CON, it ain't a level-II battery
psych <- list( memory = paste0( "avlt_", c("1_5","6","8","r_fp","r_fn") ),
               attention = c( "tmt_a", paste0("stroop_", c("body","slova") ) ),
               executive = c( "tmt_b", "stroop_barvy" ),
               speed = paste0( "gpt_", c("phk","lhk") )
               )

# read the data set
d1 <- read.csv( "data//pd_datscan.csv" , sep = "," )
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
          GrAHi, Gender.0.F, AGE, Education, # AHI grouup and demographics
          starts_with( "DaTSC_") # datscan
          ) %>%
  
  # pre-process predictor terms
  mutate_if( is.character, as.factor ) %>%
  mutate_if( grepl( "DaTSC", names(.) ), function(x) scale(x) ) %>%
  mutate( Gender.0.F = as.factor( Gender.0.F), Education = scale( Education ), AGE = scale(AGE) ) %>%
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(GrAHi) <- 1-contr.treatment(2) # Low = 0, High = 1
    contrasts(Gender.0.F) <- contr.treatment(2) # female = 0, male = 1
    } ) %>%

  # log-transform reaction times, false positives/negatives
  mutate_if( grepl( "tmt|stroop|gpt", names(.) ), function(x) log(x) ) %>%
  mutate_if( grepl( "avlt_r", names(.) ), function(x) log(x+1) )

# check for missing values
sapply( 1:nrow(df), # loop through all the rows (subjects)
        function(i) if ( sum( is.na(df[i,]) ) > 0 ) df[i,] # print the full entry of each subject with missing value
        ) %>%
  do.call( rbind.data.frame , . ) # put all rows with missing values together


# ---- computation: causal assumptions ----

# prepare a list for DAG figures
f.dag <- list()

# three plausible DAGs ()
dag <- list( dagitty( "dag { datsc <- ahi <- age -> cog <- edu <- sex -> ahi -> cog <- datsc <- age -> cog <- sex }" ),
             dagitty( "dag { datsc -> ahi <- age -> cog <- edu <- sex -> ahi -> cog <- datsc <- age -> cog <- sex }" ),
             dagitty( "dag { datsc <-> ahi <- age -> cog <- edu <- sex -> ahi -> cog <- datsc <- age -> cog <- sex }" )
             )

# prepare DAG figures including the DAG (first row) and adjustment sets for the total (second row) and
# direct (third row) effect of the regression of outcome y = {cog} on the set of predictors x = {ahi,age}
for ( i in 1:length(dag) ) {
  
  # add coordinates
  dag[[i]] <- dag[[i]] %>% `coordinates<-`( list( x = c(sex = 0, ahi = 0, datsc = 1, age = 1, edu = 2, cog = 2),
                                                  y = c(sex = 4, ahi = 3, datsc = 2, age = 1, edu = 4, cog = 3) )
                                            )
  
  # prepare a list for each DAG separately
  f.dag[[i]] <- list()
  
  # plot the full DAG
  f.dag[[i]][["full"]] <- dag[[i]] %>% ggdag() + theme_dag()
  
  # plot adjustment sets
  for( j in c("total","direct") ) f.dag[[i]][[j]] <- dag[[i]] %>%
    ggdag_adjustment_set( exposure = c("ahi","age"), outcome = "cog", effect = j, type = "minimal", shadow = F ) +
    theme_dag() + theme( legend.position = "none" )
  
}

# arrange the DAGs
( f.dag[[1]]$full | f.dag[[2]]$full | f.dag[[3]]$full ) /
( f.dag[[1]]$total | f.dag[[2]]$total | f.dag[[3]]$total ) /
( f.dag[[1]]$direct | f.dag[[2]]$direct | f.dag[[3]]$direct ) +
  plot_annotation( tag_levels = "A" )

# save it
ggsave( "figures/some_dags.jpg", dpi = 300, width = 1.1*9.64, height = 1.1*11.7 )


# ---- representation: statistical model  ----

# based on the DAGs above, adjusting for sex and datsc will lead to unbiased estimate of the joint direct effect
# of age and ahi on cognition, the total effect of the same set of predictors can be estimated only if the DAG #1
# is the best via adjusting for sex only (and possibly education, see below),
# adjusting for education won't hurt causal identification according to any of the DAGs but it can increase
# precision of ACE estimate (see Model 8 in https://www.doi.org/10.1177/00491241221099552 )

# fit a series of univariate regressions with neuropsychology outcomes
# prepare a list for formulas
f <- list()

# set-up the linear models
for ( i in names(psych) ) {
  for ( j in psych[[i]] ) f[[i]][[j]] <- paste0( j, " ~ AGE * GrAHi + DaTSC_striatum_worse + Gender.0.F + Education" ) %>% as.formula()
}


# ---- implementation: learning from the data  ----

# use the default lm() function to fit the regressions via the QR-decomposition and least squares
m <- list()

# loop through all linear models in f, select appropriate likelihoods and links (via family argument)
for ( i in names(f) ) {
  for ( j in names(f[[i]]) ) m[[i]][[j]] <- lm( formula = f[[i]][[j]], data = df )
}


# ---- post-processing: model summaries ----

# prepare a list for model parameters and model checks
t <- list()

# fill-in all the values
for ( i in names(m) ) {
  for ( j in names(m[[i]]) ) t[[i]][[j]] <- cbind.data.frame(
    
    # number of patients
    rep( ( !is.na( df[ , j ] ) ) %>% sum(), 3 ),
    
    # model parameters and associated statistics
    ( summary( m[[i]][[j]] ) )$coefficients[ c("AGE","GrAHi2","AGE:GrAHi2"), ],
   
    # empty column for a statistical significance (raw, i.e. adjusted for per comparison error rate, PCER,
    # and after Benjamini-Hochberg correction for false discovery rate, FDR) and model checks
    matrix(NA, nrow = 3, ncol = 9 )
    
  ) %>%
      
      # give proper names to each column
      `names<-`( c("N","Effect size","Std. Error","t value","Pr(>|t|)","sig. (PCER = 5%)","sig. (FDR = 5%)",
                   "normality (p-value)","non-normality","homoscedasticity (p-value)","heteroscedasticity",
                   "outliers", "VIF (max)","multicollinearity") ) %>%
      
      # add a column including relevant prediction terms
      rownames_to_column( var = "Predictor term" )
  # collapse tables for a single outcome to a neat table
  t[[i]] <- do.call( rbind.data.frame, t[[i]] ) %>%
    as.data.frame() %>% rownames_to_column( var = "Outcome" ) %>% # add column denoting the outcome
    mutate( Outcome = sub( "\\..*", "", Outcome) ) # tidy the outcome variable
  
}

# pull all model summaries into a one neat table
t <- do.call( rbind.data.frame, t ) %>%
  rownames_to_column( var = "Domain" ) %>%
  mutate( Domain = sub( "\\..*", "", Domain) )

# add model checks
for ( i in 1:nrow(t) ) t[ i , c("normality (p-value)","homoscedasticity (p-value)","outliers","VIF (max)") ] <- c(
  with( t, check_normality( m[[Domain[i]]][[Outcome[i]]] ) ) %>% as.numeric(), # p ≤ .05 rejects normality
  with( t, check_heteroskedasticity( m[[Domain[i]]][[Outcome[i]]] ) ) %>% as.numeric(), # p ≤ .05 rejects homogeneity of variances
  with( t, check_outliers( m[[Domain[i]]][[Outcome[i]]] ) ) %>% sum(), # number of outliers according to the default setting in the "performance" package
  with( t, check_collinearity( m[[Domain[i]]][[Outcome[i]]] ) )[ , "VIF" ] %>% max() # maximal VIF
)

# extract the Benjamini-Hochberg corrected threshold
bh_thres <- data.frame( p = sort( t$`Pr(>|t|)`), # order the p-values from lowest to largest
                        thres = .05 * (1:nrow( t[complete.cases(t$`Pr(>|t|)`),] ) ) / nrow( t[complete.cases(t$`Pr(>|t|)`),] ) # prepare BH thresholds for each p-value
                        ) %>%
  # flag BH-significant p-values and extract the largest threshold (https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
  mutate( sig = ifelse( p <= thres, T, F ) ) %>% filter( sig == T ) %>% select(thres) %>% max()

# add final touches, i.e., statistical significance decisions for ANOVAs (based on BH correction),
# non-normality (based on nominal p < .05) and heteroscedasticity (based on nominal p < .05)
t <- t %>% mutate( `non-normality` = ifelse( `normality (p-value)` < .05, "!", NA ),
                   heteroscedasticity = ifelse( `homoscedasticity (p-value)` < .05, "!", NA ),
                   multicollinearity = ifelse( `VIF (max)` > 5, "!", NA ),
                   `sig. (PCER = 5%)` = ifelse( `Pr(>|t|)` < .05, ":-)", NA ),
                   `sig. (FDR = 5%)` = ifelse( `Pr(>|t|)` < bh_thres, ":-)", NA )
                   )

# print as a csv (use the second line if demographics were adjusted for)
write.table( t, "tables/pd_only_direct_effects.csv", sep = ",", row.names = F, na = "" )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/pd_only_regrs.txt" )