#' Convert Nest Counts to Abundance
#' @param df Cleaned dataframe from prep_nesting_data
#' @param clutch_freq Average nests per female per season
#' @param remig_int Average years between nesting seasons
#' @export
calculate_abundance <- function(df, clutch_freq = NULL, remig_int = NULL) {
  
  cat("\n--- Abundance Conversion ---\n")
  
  if(is.null(clutch_freq)) {
    clutch_freq <- as.numeric(readline("Nests per female per season (e.g. 5.5): "))
  }
  
  if(is.null(remig_int)) {
    remig_int <- as.numeric(readline("Remigration interval in years (e.g. 3.0): "))
  }
  
  res <- df %>%
    dplyr::mutate(
      Annual_Nesters = Count / clutch_freq,
      Total_Adult_Females = (Count / clutch_freq) * remig_int
    )
  
  message("Success: Abundance metrics calculated.")
  return(res)
}