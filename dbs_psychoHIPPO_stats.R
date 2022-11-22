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
  
  # set-up contrasts to avoid multicollinearity in interaction terms
  within( . , {
    contrasts(GROUP) <- contr.treatment(2) # CON = 0, PD = 1
    contrasts(AHI_LH) <- contr.treatment(2) # High = 0, Low = 1
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
  for ( j in psych[[i]] ) f[[i]][[j]] <- list(
    intercept = paste0( j, " ~ 1") %>% as.formula(), # 1) intercept only
    diagnosis = paste0( j, " ~ 1 + GROUP" ) %>% as.formula(), # 2) effect of diagnosis
    ahi = paste0( j, "~ 1 + GROUP * AHI_LH" ) %>% as.formula(), # 3) added effect of AHI
    hippo =  paste0( j, " ~ 1 + GROUP * AHI_LH * HBT.R.BODY + GROUP * AHI_LH * HBT.L.BODY" ) # 4) mediation via hippocampal volume
  )
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


# ---- post-processing: model checks ----

# prepare a list for model checks
check <- list()

# fill-in the checks
for ( i in names(m) ) {
  for ( j in names(m[[i]]) ) {
    
    # first fill-in all diagnostics we care about at this point
    for ( k in names(m[[i]][[j]]) ) check[[i]][[j]][[k]] <- c(
      normality = check_normality( m[[i]][[j]][[k]] ) %>% as.numeric(), # p ≤ .05 rejects normality
      homoscedasticity = check_heteroskedasticity( m[[i]][[j]][[k]] ) %>% as.numeric(), # p ≤ .05 rejects homogeneity of variances
      outliers = check_outliers( m[[i]][[j]][[k]] ) %>% sum() # number of outliers according to the default setting in the "performance" package
    )
    
    # collapse tables for a single outcome to a neat table
    check[[i]][[j]] <- do.call( cbind.data.frame, check[[i]][[j]] ) %>%
      t() %>% # flip dimensions
      as.data.frame() %>% rownames_to_column( var = "predictors" ) # add column denoting model type
  }
  
  # collapse tables in a single domain to an even neater table
  check[[i]] <- do.call( rbind.data.frame, check[[i]] ) %>%
    rownames_to_column( var = "outcome" ) %>% # add column denoting model type
    mutate( outcome = sub( "\\..*", "", outcome) ) # tidy the outcome variable
  
}

# pull all models' diagnostics into one and flag potentially problematic models
check <- do.call( rbind.data.frame, check ) %>%
  rownames_to_column( var = "domain" ) %>% # add domain names
  mutate( domain = sub( "\\..*", "", domain), # tidy the domain column
          `non-normality` = ifelse( normality <= .05, 1, 0 ),
          heteroscedasticity = ifelse( homoscedasticity <= .05, 1, 0 ),
          )

# print as a csv
write.table( check, "tables/model_checks.csv", sep = ",", row.names = F )


# ---- post-processing: models' comparisons (via ANOVA) ---

# prepare a list for "stepwise" ANOVA results
aov <- list()

# fill-in the stats
for ( i in names(m) ) {
  for ( j in names(m[[i]]) ) {
    
    # first fill-in all comparisons for a single outcome
    for ( k in names(m[[i]][[j]]) ) aov[[i]][[j]][[k]] <- c(
      ( !is.na( df[ df$GROUP == "CON", j ] ) ) %>% sum(), # number of controls
      ( !is.na( df[ df$GROUP == "PD", j ] ) ) %>% sum(), # number of patients
      case_when(
        k == "intercept" ~ rep(NA, 4) %>% as.numeric(),
        k == "diagnosis" ~ with( m[[i]][[j]], anova( intercept, diagnosis) )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric(),
        k == "ahi" ~ with( m[[i]][[j]], anova( diagnosis, ahi ) )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric(),
        k == "hippo" ~ with( m[[i]][[j]], anova( ahi, hippo ) )[ 2, c("F","Df","Res.Df","Pr(>F)") ] %>% as.numeric()
      )
    ) %>% `names<-`( c("nCON", "nPD", "F", "df1", "df2", "Pr(>F)") )
    
    # collapse tables for a single outcome to a neat table
    aov[[i]][[j]] <- do.call( cbind.data.frame, aov[[i]][[j]] ) %>%
      t() %>% # flip dimensions
      as.data.frame() %>% rownames_to_column( var = "predictors" ) # add column denoting model type
  }
  
  # collapse tables in a single domain to an even neater table
  aov[[i]] <- do.call( rbind.data.frame, aov[[i]] ) %>%
    rownames_to_column( var = "outcome" ) %>% # add column denoting model type
    mutate( outcome = sub( "\\..*", "", outcome) ) # tidy the outcome variable
  
}

# pull all ANOVAs into one table and flag significant comparisons after BH correction
aov <- do.call( rbind.data.frame, aov ) %>%
  rownames_to_column( var = "domain" ) %>% # add domain names
  mutate( domain = sub( "\\..*", "", domain) # tidy the domain column
  )

# extract the Benjamini-Hochberg corrected threshold
bh_thres <- data.frame( p = sort( aov$`Pr(>F)`), # order the p-values from lowest to largest
                        thres = .05 * ( (1:nrow(na.omit(aov)))/nrow(na.omit(aov)) ) # prepare BH thresholds for each p-value
                        ) %>%
  # flag BH-significant p-values and extract the largest threshold as per https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
  mutate( sig = ifelse( p <= thres, T, F ) ) %>% filter( sig == T ) %>% select(thres) %>% max()

# flag BH-corrected significant comparisons in the aov table
aov <- aov %>% mutate( sig = ifelse( `Pr(>F)` <= bh_thres, "*", "" ) )
