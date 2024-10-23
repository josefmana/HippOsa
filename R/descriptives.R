# Describes the sample of de novo PD patients and HC for the study of association between hippocampus and apnea

describe_data <- function(data) {
  
  # list variables to be described
  cont <- c("AGE","EDU.Y","BMI","age_first_symptom","disease_duration","moca","mds_updrs_i","mds_updrs_ii","mds_updrs_iii_total","mds_updrs_iii_axial","mds_updrs_iii_rigidityakineasia","mds_updrs_iii_tremor")
  nomin <- c("GENDER","RBD")
  
  # pre-process data
  d0 <- data %>%
    
    # prepare GENDER and iRBD as factors
    mutate( across(all_of( c(nomin,"SUBJ","AHI.F") ), as.factor) ) %>%
    within( . , {
      contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
      contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
    } )
  
  # prepare table with descriptions
  tab1 <- d0 %>%
    
    group_by(GROUP) %>%
    summarise(
      across( all_of(nomin), freqprop ),
      across( all_of(cont), msd )
    ) %>%
    ungroup() %>%
    column_to_rownames("GROUP") %>%
    t() %>%
    as.data.frame() %>%
    mutate( across( everything(), ~ if_else(grepl("NaN",.x), "-", .x) ) ) %>%
    rownames_to_column("y")
  
  # add statistical analysis results of variables present in both PD and HC
  tab1 <- left_join(
    
    tab1, lapply(
      
      cont[c(1:3,6)],
      function(y)
        
        summary( lm(as.formula( paste0(y," ~ SUBJ * AHI.F") ), data = d0) )$coefficients[ , c("t value","Pr(>|t|)") ] %>%
        statextract(y = y, stat = "t")
      
    ) %>% do.call( rbind.data.frame, . ),
    
    by = "y"
    
  )
  
  # add logistic regression for gender
  tab1[tab1$y == "GENDER", c("SUBJ1","AHI.F1","SUBJ1:AHI.F1")] <-
    
    summary( glm(GENDER ~ SUBJ * AHI.F, family = binomial(), data = d0) )$coefficients[ , c("z value","Pr(>|z|)") ] %>%
    statextract(y = "GENDER", stat = "z") %>%
    select(-y)
  
  # add logistic regression of iRBD
  tab1[tab1$y == "RBD", "AHI.F1"] <-
    
    summary( glm( RBD ~ AHI.F, family = binomial(), data = d0 ) )$coefficients[ , c("z value","Pr(>|z|)") ] %>%
    statextract(y = "RBD", stat = "z") %>%
    select(-y)
  
  # add results of continuous variables for PD only
  for( i in with( tab1, y[CON_H == "-"] )[-1] ) tab1[ tab1$y == i, "AHI.F1"] <-
    
    summary( lm( as.formula( paste0(i," ~ AHI.F") ), data = d0 ) )$coefficients[ , c("t value","Pr(>|t|)") ] %>%
    statextract(y = i, stat = "t") %>%
    select(-y)
  
  # return the table
  return(tab1)
  
}
