# Imports and pre-process data to analysis-ready format

rm( list = ls() ) # clear environment

# load packages
library(here)
library(openxlsx)
library(tidyverse)

if( !dir.exists("_data") ) dir.create("_data")

# all tests in the battery
psych <- read.csv( here("helpers","psychs.csv"), sep = ";")


# the primary data

d1 <- left_join(
  
  # read hippocampus data first
  read.xlsx( here("_raw","Tabhipposubfields_rhx.xlsx") ) %>% rename_with( ~ gsub("-", "_", .x) ),
  read.xlsx( here("_raw","Tabhipposubfields_lhx.xlsx") ) %>% rename_with( ~ gsub("-", "_", .x) ),
  by = "Study.ID",
  suffix = c("_rhx","_lhx") # add suffixes to differentiate betwen right and left sides
  
) %>%
  
  # add all other subscortical structures and cognitive data
  left_join(
    
    # start with the subcortical data as they include demography and such
    read.xlsx( here("_raw","Tabsubcortex1_RBD.xlsx") ) %>% rename_with( ~ gsub("-", "_", .x) ),
    . ,
    by = "Study.ID"
    
  ) %>%
  left_join(
    
    read.csv( here("_raw","RBDBIOPDCON_DATA_2024-07-17_1146.csv"), sep = "," ) %>%
      mutate( Study.ID = sub("-","",study_id) ) %>%
      mutate( event = sub("_arm_2|_arm_3", "", redcap_event_name) ) %>%
      select( Study.ID, event, all_of( c(psych$variable,"moca") ) ),
    by = "Study.ID"
    
  ) %>%
  
  # correct formatting errors
  mutate( across( ends_with("_rhx"), as.numeric ) ) %>%
  
  # add MDS-UPDRS III
  left_join(
    
    read.xlsx( here("_raw","BIOPD_MDSUPDRSIII.xlsx"), startRow = 2, check.names = T ) %>%
      mutate(
        age_first_symptom = time_length(difftime(První.motorický.příznak.PN..MM.RRRR., Datum.narození), "years"),
        disease_duration = time_length(difftime(Datum.vyšetření, První.motorický.příznak.PN..MM.RRRR.), "years"),
        event = "enrollment"
      ) %>%
      select(Study.ID, event, age_first_symptom, disease_duration, MDS.UPDRS.Část.I.sumární.skóre, MDS.UPDRS.Část.II.sumární.skóre..počítaná.hodnota., MDS.UPDRS.Část.III..celkové.score., Axiasubscore, Rigidityakinesia, Tremorsubscoret) %>%
      rename(
        "mds_updrs_i" = "MDS.UPDRS.Část.I.sumární.skóre",
        "mds_updrs_ii" = "MDS.UPDRS.Část.II.sumární.skóre..počítaná.hodnota.",
        "mds_updrs_iii_total" = "MDS.UPDRS.Část.III..celkové.score.",
        "mds_updrs_iii_axial" = "Axiasubscore",
        "mds_updrs_iii_rigidityakineasia" = "Rigidityakinesia",
        "mds_updrs_iii_tremor" = "Tremorsubscoret"
      ),
    
    by = c("Study.ID","event")

  )


# save the data
write.table(d1, file = here("_data","primary_dataset.csv"), sep = ",", row.names = F, quote = F)
