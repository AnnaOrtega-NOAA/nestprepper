#' Run Full NestPrepper Workflow
#'
#' @description Chaining the data prep, abundance calculation, and model run into one command.
#' @param df A dataframe with Year, Site, and Count columns.
#' @param clutch_freq Average nests per female per season.
#' @param remig_int Average years between nesting seasons.
#' @param iter Number of MCMC iterations.
#' @param quiet Logical. If TRUE, skips interactive console prompts.
#' @return A list containing the model fit and the processed abundance data.
#' @export
run_nestprepper_workflow <- function(df, clutch_freq = 5.5, remig_int = 3.0, iter = 2000, quiet = TRUE) {

  # 1. QAQC and Cleaning
  df_clean <- prep_nesting_data(df, quiet = quiet)

  # 2. Abundance Calculation
  df_abund <- calculate_abundance(
    df_clean,
    clutch_freq = clutch_freq,
    remig_int = remig_int,
    quiet = quiet
  )

  # 3. Bayesian Model Execution
  model_results <- run_turtle_model(df_abund, iter = iter)

  return(list(
    res = model_results,
    abund = df_abund
  ))
}
