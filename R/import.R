# z-scores calculation
zscore <- function(calc, x, nam, AGE, GEN, EDU) with(
  
  calc, {
    
    # prepare a matrix of data and vector of parameters
    pars <- as.numeric( c( Constant[var == nam], age[var == nam], gender[var == nam], education[var == nam] ) )
    data <- as.matrix( cbind( rep( 1, length(x) ), AGE, GEN, EDU ) )
    
    # compute the quantities
    x_bar <- data %*% pars
    z <- c( sign[var == nam] * (x - x_bar) / RMSE[var == nam] )
    return(z)
    
  }
)


# List paths to data ----
data_paths <- function() data.frame(
  
  lhippo = here("_raw", "Tabhipposubfields_lhx.xlsx"),
  rhippo = here("_raw", "Tabhipposubfields_rhx.xlsx"),
  subcor = here("_raw", "Tabsubcortex1_RBD.xlsx"),
  psych = here("_raw", "RBDBIOPDCON_DATA_2024-07-17_1146.csv"),
  motor = here("_raw", "BIOPD_MDSUPDRSIII.xlsx")
  
)


# List paths to helper files ----
extract_helpers <- function() with(
  
  data.frame(
    psychvar = here("helpers","psychs.csv"),
    calculat = here("_raw","calculator_final_v7_c_301116.xlsx"),
    calc_sheet = "equations",
    subcortex = here("helpers","subcortical.csv"),
    hippocampi = here("helpers","hippocampus.csv")
  ),
  
  list(
    psychohelp = read.csv(psychvar, sep = ";"), # psychological variables
    calculator = read.xlsx(calculat, sheet = calc_sheet, startRow = 2),
    psych = read.csv(psychvar, sep = ";") %>% filter( !is.na(domain) ), # psychological variables for PD vs CON comparisons
    subco = read.csv(subcortex, sep = ","), # subcortical structures
    hippo = read.csv(hippocampi, sep = ",") %>% filter( complete.cases(name) ) # hippocampal structures
  )
  
)


# Extract response time variable names to be log-transformed ----
extract_rt_variables <- function(helpers) with(
  
  helpers,
  subset(psych, domain %in% c("Attention","Executive function","Processing speed") )$variable
  
)


# Imports and pre-process data to analysis-ready format ----
import_data <- function(files, helpers) {
  
  # read raw data
  data <- with(
    
    files, list(
      
      # MRI volummetry estimates
      hippo = left_join(
        
        read.xlsx(rhippo) %>% rename_with( ~ gsub("-", "_", .x) ) %>% mutate( across(!ends_with("ID"), as.numeric) ),
        read.xlsx(lhippo) %>% rename_with( ~ gsub("-", "_", .x) ),
        by = "Study.ID",
        suffix = c("_rhx", "_lhx")
        
      ),
      subcor = read.xlsx(subcor) %>% rename_with( ~ gsub("-", "_", .x) ),
      
      # psychology data
      psych =
        read.csv(psych, sep = ",") %>%
        mutate(
          Study.ID = sub("-","",study_id),
          event = sub("_arm_2|_arm_3", "", redcap_event_name)
        ) %>%
        select( Study.ID, event, all_of( c(helpers$psychohelp$variable,"moca") ) ) %>%
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
    helpers$calculator %>%
    filter( X1 %in% with( helpers$psychohelp, label[!is.na(pairs)] ) ) %>%
    mutate(
      sign = if_else(X1 %in% c("TMT-A","PST-C"), -1, 1),
      var = unlist(
        sapply( 1:nrow(.), function(i) helpers$psychohelp[helpers$psychohelp$label == X1[i], "variable"] ),
        use.names = T
      )
    )
  
  # prepare thresholds for BNT-60
  bnt_thresh <- data.frame(
    age_bottom = c(0, 0, 60, 60),
    age_top = c(60, 60, Inf, Inf),
    edu_bottom = c(0, 12, 0, 12),
    edu_top = c(12, Inf, 12, Inf),
    threshold = c(49, 52, 50, 53)
  )

  # add MCI flag for each test of each patient
  mci <-
    
    # prepare data
    data %>%
    filter(SUBJ == "PD") %>% # keep patients only
    
    # add MCI flags (i.e., z ≤ -1.5)
    mutate(
      across(
        .cols = calculator$var,
        .fns = ~ zscore(calculator, .x, cur_column(), AGE, GENDER, EDU.Y),
        .names = "z_{.col}"
      ),
      across(
        .cols = starts_with("z_"),
        .fns = ~ if_else(.x <= -1.5, 1, 0),
        .names = "mci_{.col}"
      ),
      mci_bnt = sapply(
        1:nrow(.),
        function(i)
          ifelse(
            test = bnt60[i] <= with(bnt_thresh, threshold[AGE[i] > age_bottom & AGE[i] <= age_top & EDU.Y[i] > edu_bottom & EDU.Y[i] <= edu_top] ),
            yes = 1,
            no = 0
          )
      )
    ) %>%
    
    # evaluate PD-MCI
    select(Study.ID, starts_with("mci_") ) %>%
    mutate(
      flags = rowSums( across( starts_with("mci_") ), na.rm = T),
      nas = rowSums( is.na( across( starts_with("mci_") ) ) ),
      pd_mci = case_when(
        flags >= 2 ~ 1,
        flags == 0 & nas < 2 ~ 0,
        flags == 1 & nas == 0 ~ 0,
        .default = NA
      )
    ) %>%
    select(Study.ID, pd_mci)
  
  # add PD-MCI to the data and return
  return( left_join(data, mci) )

}


# Pre-processes the data ----
preprocess_data <- function(data, help, rt_vars, return = "df") {
  
  # read the data
  d0 <- data %>%
    
    # re-calculate hippocampal fields according to the legend in hippo
    mutate(
      !!!setNames( rep(NA, length( unique(help$hippo$name) ) ), unique(help$hippo$name) ),
      across(
        .cols = unique(help$hippo$name),
        .fns = ~ rowSums( across( all_of( with( help$hippo, var[name == cur_column()] ) ) ) )
      )
    )
  
  # format it for analyses
  df <- d0 %>%
    
    # pre-process brain and demography variables
    mutate_if(is.character, as.factor) %>%
    mutate(
      PD = if_else(SUBJ == "PD", 1, 0),
      GENDER = as.factor(GENDER),
      across( all_of(rt_vars), ~ -log(.x) ) # cognition
    ) %>%
    
    # set-up contrasts to avoid multicollinearity in interaction terms
    within( . , {
      contrasts(SUBJ) <- -contr.sum(2)/2 # CON = -0.5, PD = 0.5
      contrasts(AHI.F) <- -contr.sum(2)/2 # High = -0.5, Low = 0.5
      contrasts(GENDER) <- contr.sum(2)/2 # female = 0.5, male = -0.5
    } )
  
  # extract scaling values, i.e.,
  # enrollment full sample means and SDs
  scl <- sapply(
    
    c("AGE", "EDU.Y", "BMI", "SBTIV", help$subco$name, unique(help$hippo$name), help$psych$variable),
    function(i)
      with(df, c(M = mean( get(i), na.rm = T ), SD = sd( get(i), na.rm = T ) ) )
    
  ) %>% t()
  
  # scale the variables
  df <- df %>%
    
    mutate(
      across(
        all_of( rownames(scl) ),
        ~ (.x - scl[cur_column(), "M"]) / scl[cur_column(), "SD"]
      )
    )
  
  # return the data
  return( get(return) )
  
}