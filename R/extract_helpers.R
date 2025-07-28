#' Prepare 'helper' files
#' 
#' This function takes no parameter but needs to be adjusted manually if desired.
#' It is used in the \pkg{targets} pipeline to check for changes in files that
#' do not contain data proper but 'helper' data, i.e., data mapping raw data
#' structures to functions used throughout the pipeline.
#' 
#' @exports A list with data frames containg different kinds of helper files.
extract_helpers <- function() {
  with(data.frame(
    psychvar = here("helpers", "psychs.csv"),
    calculat = here("_raw", "calculator_final_v7_c_301116.xlsx"),
    calc_sheet = "equations",
    subcortex = here("helpers", "subcortical.csv"),
    hippocampi = here("helpers", "hippocampus.csv")
  ),
  list(
    psychohelp = read.csv(psychvar, sep = ";"), # psychological variables
    calculator = read.xlsx(calculat, sheet = calc_sheet, startRow = 2),
    psych = read.csv(psychvar, sep = ";") |> filter(!is.na(domain)), # psychological variables for PD vs CON comparisons
    subco = read.csv(subcortex, sep = ","), # subcortical structures
    hippo = read.csv(hippocampi, sep = ",") |> filter(complete.cases(name)) # hippocampal structures
  ))
}
