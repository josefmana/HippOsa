# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "dplyr", "tidyverse", # for data wrangling
           "ggdag", "dagitty", # for DAGs
           "brms", "tidybayes", # for Bayesian analyses 
           "ggplot2", "patchwork" # for plotting
           )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set some values for later
s = 87542 # seed for reproducibility

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .95 # adapt_delta parameter

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_minimal(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# list all tests in the battery
bat <- list( att = c("tmt_a","wais_iii_ciselny_pozadu"), # attention and working memory
             ef = c("tol","vf_zvi"), # executive functions
             lan = c("wais_iii_podob","bnt60"), # language
             mem = c("avlt_8","bvmt_r_delay"), # memory
             vs = c("jolo","clox_i") # visuospatial
             )

# read the data set
d1 <- read.csv( "data/pd_con.csv" , sep = "," )
d2 <- read.csv( "data/pd_datscan.csv" , sep = "," )
d3 <- read.csv( "data/20221109_redcap_export.csv", sep = "," )


# ---- pre-processing  ----

# check there are the same IDs in d1 and d2
isTRUE( all.equal( with( d1, Study.ID[ grepl( "BIOPD", Study.ID) ] ), unique(d2$Study.ID) ) ) # A-OK

# keep only included patients in d3
d3 <- d3 %>%
  mutate( Study.ID = sub( "-", "", study_id ) ) %>%
  slice( which( Study.ID %in% unique(d2$Study.ID) ) )

# prepare a data frame for analyses
df <- d2 %>% left_join( d3, by = "Study.ID" ) %>%
  rename( id = Study.ID, age_y = AGE, ahi = AHI.value, edu_y = Education ) %>%
  mutate( ahi_gr = ifelse( AHI == 0, "low", "high" ), sex = ifelse( Gender.0.F == 0, "f", "m") ) %>%
  mutate_if( is.character, as.factor )

# check for missing values
sapply( names(df), function(i) sum( is.na(df[[i]]) ) ) # two in avlt (one ain't speaking Czech)


# ---- computation: causal model  ----

# let`s draw some
dag <- dagitty( "dag {
                hippo -> cogn <- hippo <- ahi -> cogn <- datsc
                hippo <-  ahi -> cogn
                datsc -> ahi
                cogn [outcome]
                }" )

# print it
ggdag( tidy_dagitty(dag) ) + theme_dag()


# ---- representation: statistical model  ----

# outcome: psycho
# exposure: pd * ahi * hippo_body (L vs R)


# ---- implementation: learning from the data  ----