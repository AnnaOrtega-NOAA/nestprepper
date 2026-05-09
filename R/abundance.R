#' Convert Nest Counts to Abundance
#' @param df Cleaned dataframe from prep_nesting_data (Year, Site, Count)
#' @param clutch_freq Average nests per female per season
#' @param remig_int Average years between nesting seasons
#' @param quiet Logical. If TRUE, skips interactive console prompts.
#' @export
calculate_abundance <- function(df, clutch_freq = NULL, remig_int = NULL, quiet = FALSE) {

  if(!quiet) cat("\n--- Abundance Conversion ---\n")

  if(is.null(clutch_freq)) {
    if(quiet) stop("clutch_freq must be provided in non-interactive mode.")
    clutch_freq <- as.numeric(readline("Nests per female per season (e.g. 5.5): "))
  }

  if(is.null(remig_int)) {
    if(quiet) stop("remig_int must be provided in non-interactive mode.")
    remig_int <- as.numeric(readline("Remigration interval in years (e.g. 3.0): "))
  }

  res <- df %>%
    dplyr::mutate(
      Annual_Nesters = Count / clutch_freq,
      Total_Adult_Females = Annual_Nesters * remig_int
    )

  if(!quiet) message("Success: Abundance metrics calculated.")
  return(res)
}
