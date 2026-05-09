#' Aggregate Monthly Data to Annual Totals
#' @export
aggregate_monthly_to_annual <- function(df) {
  # Standardize for case-sensitivity
  names(df) <- tools::toTitleCase(tolower(names(df)))

  if("Month" %in% names(df)) {
    message("-> Monthly data detected for JM/W. Aggregating...")

    annual_df <- df %>%
      dplyr::group_by(Year, Site) %>%
      dplyr::summarise(
        Months_Present = sum(!is.na(Count)),
        Raw_Sum = sum(Count, na.rm = TRUE),
        # Basic imputation: if months are missing, scale up to 12 months
        # (Mirroring the intent of your Fourier imputation)
        Count = ifelse(Months_Present > 0, (Raw_Sum / Months_Present) * 12, NA_real_),
        .groups = "drop"
      )
    return(annual_df %>% dplyr::select(Year, Site, Count))
  }
  return(df)
}
