# Imports and pre-process data to analysis-ready format

import_data <- function(files, helpers) {
  
  # helpers
  with(
    helpers, {
      psychohelp <<- read.csv(psychvar, sep = ";")
      calculator <<- read.xlsx(calculat, sheet = calc_sheet, startRow = 2)
    }
  )
  
  # read raw data
  data <- with(
    
    files, list(
      
      # MRI volummetry estimates
      rhippo =
        read.xlsx(rhippo) %>%
        rename_with( ~ gsub("-", "_", .x) ) %>%
        mutate( across(!ends_with("ID"), as.numeric) ),
      lhippo = read.xlsx(lhippo) %>% rename_with( ~ gsub("-", "_", .x) ),
      subcor = read.xlsx(subcor) %>% rename_with( ~ gsub("-", "_", .x) ),
      
      # psychology data
      psych =
        read.csv(psych, sep = ",") %>%
        mutate(
          Study.ID = sub("-","",study_id),
          event = sub("_arm_2|_arm_3", "", redcap_event_name),
          fluency = vf_zvi + vf_oble + vf_obch,
        ) %>%
        select( Study.ID, event, all_of( c(psychohelp$variable,"moca") ) ) %>%
        filter(event == "enrollment"),
      
      # MDS-UPDRS data
      motor =
        read.xlsx(motor, startRow = 2, check.names = T) %>%
        mutate(
          age_first_symptom = time_length(difftime(První.motorický.příznak.PN..MM.RRRR., Datum.narození), "years"),
          disease_duration = time_length(difftime(Datum.vyšetření, První.motorický.příznak.PN..MM.RRRR.), "years"),
          event = "enrollment"
        ) %>%
        select(
          Study.ID, event, age_first_symptom, disease_duration,
          MDS.UPDRS.Část.I.sumární.skóre,
          MDS.UPDRS.Část.II.sumární.skóre..počítaná.hodnota.,
          MDS.UPDRS.Část.III..celkové.score., Axiasubscore, Rigidityakinesia, Tremorsubscoret
        ) %>%
        rename(
          "mds_updrs_i" = "MDS.UPDRS.Část.I.sumární.skóre",
          "mds_updrs_ii" = "MDS.UPDRS.Část.II.sumární.skóre..počítaná.hodnota.",
          "mds_updrs_iii_total" = "MDS.UPDRS.Část.III..celkové.score.",
          "mds_updrs_iii_axial" = "Axiasubscore",
          "mds_updrs_iii_rigidityakineasia" = "Rigidityakinesia",
          "mds_updrs_iii_tremor" = "Tremorsubscoret"
        )
      
    )
  ) %>%
    
    # pull the data to a single file
    reduce(left_join)
  
  # prepare objects for PD-MCI labelling
  calculator <-
    calculator %>%
    filter( X1 %in% with( psychohelp, label[!is.na(pairs)] ) )
  
  # prepare threholds for BNT-60
  bnt_thresh <- data.frame(
    age_bottom = c(0, 0, 60, 60),
    age_top = c(60, 60, Inf, Inf),
    education = c("low", "high", "low", "high"),
    threshold = c(49, 52, 50, 53)
  )
  
  # 

}
