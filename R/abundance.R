#' Convert Nest Counts to Abundance (Nesters & Females)
#' @param df Dataframe from prep_nesting_data (Year, Site, Count)
#' @param clutch_freq Nests per female per season
#' @param remig_int Average years between nesting seasons
#' @export
calculate_abundance <- function(df, clutch_freq = 5.5, remig_int = 3.06, quiet = FALSE) {

  if(!quiet) message("--- Standard Abundance Conversion ---")

  res <- df %>%
    dplyr::mutate(
      Annual_Nesters = Count / clutch_freq,
      Total_Adult_Females = Annual_Nesters * remig_int
    )

  return(res)
}
