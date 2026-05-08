#' QAQC and Prep Nesting Data
#' @param df A dataframe with Year, Site, and Count columns
#' @export
prep_nesting_data <- function(df) {
  cat("\n--- Starting NestPrepper QAQC ---\n")
  
  # 1. Outlier Detection (3x IQR Rule)
  check_df <- df %>%
    dplyr::group_by(Site) %>%
    dplyr::mutate(
      med = median(Count, na.rm = TRUE),
      iqr = IQR(Count, na.rm = TRUE),
      limit = med + (3 * iqr)
    )
  
  outliers <- check_df %>% dplyr::filter(Count > limit)
  
  if(nrow(outliers) > 0) {
    message("! WARNING: Potential Outliers Detected (Magnitude Check)")
    print(outliers %>% dplyr::select(Year, Site, Count))
    ans_outlier <- readline("Is this a data entry mistake? (y = stop and fix / n = biological spike): ")
    if(tolower(ans_outlier) == 'y') stop("Process halted by user to fix data.")
  }
  
  # 2. Missing Data Questions
  if(any(is.na(df$Count))) {
    message("\n! Found missing values (NAs).")
    ans_na <- readline("Should NAs be zeros (z) or treated as unmonitored/missing (m)?: ")
    if(tolower(ans_na) == 'z') {
      df$Count[is.na(df$Count)] <- 0
      message("-> NAs converted to 0.")
    } else {
      message("-> NAs left as missing. Model will interpolate.")
    }
  }
  
  # 3. Effort Check
  effort <- readline("\nWas monitoring effort consistent across all years? (y/n): ")
  if(tolower(effort) == 'n') {
    warning("Uneven effort detected. Interpret abundance trends with caution!")
  }
  
  return(df %>% dplyr::select(Year, Site, Count))
}