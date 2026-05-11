[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The goal of `nestprepper` is to streamline the transition from messy, real-world sea turtle nesting counts to robust biological abundance estimates and Bayesian trend analyses. This package supports NOAA Open Science initiatives by providing a transparent, reproducible pipeline for population status assessments.

## Key Features

* **Smart Mapping**: Automatically handles both "Wide" format (multiple beach columns) and "Long" format (single count column) data.
* **Built-in QAQC**: Interactively flags biological spikes (outliers) and handles missing monitoring years.
* **Bayesian Engine**: Fits State-Space models to estimate population trends ($U$) and process variance ($Q$).
* **Visual Status**: Generates "Ghost Plots" that communicate population status and uncertainty in a way that is accessible to managers.

## Data Input Requirements

`nestprepper` is designed for high-level population assessments. To ensure the accuracy of Bayesian trend estimates, the package currently requires **Annual Nest Counts**. 

### **Standardized Data Structure**
For the most reliable results, your input CSV should follow a "one-row-per-year" structure. The app is optimized to recognize:

* **Year Mapping**: Automatically identifies columns named `Year`, `Season`, or `Year_begin`.
* **Count Mapping**: Users can select the specific column containing total annual nest counts (e.g., `JM_Nests` or `Nests`).
* **Site Identification**: Supports data with multiple sites/beaches by allowing the user to select the relevant site-specific count column, or multiple sites/beaches/columns.

### **Current Constraints**
* **Annual Totals Only**: The app is currently optimized for aggregated annual totals. 
* **Structure**: If your data is in a monthly format, it must be aggregated into annual sums before being uploaded into the `nestprepper` pipeline. We are working on another R package to do this!
* **Missing Values**: The Bayesian engine handles missing monitoring years (NAs), but requires a designated `Year` and `Count` column to build the time series.

```r
# Ideal structure for nestprepper:
# Year | Nests
# 2003 | 379
# 2004 | 124
# 2005 | NA (Missing monitoring year)
```

## NOAA Disclaimer
This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce. All NOAA GitHub project code is provided "as is" and the user assumes the entire risk as to its quality and performance. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the Department of Commerce.

## Installation
You can install the development version of `nestprepper` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("AnnaOrtega-NOAA/nestprepper")
```

## Quick Start (Shiny Dashboard)
For most users, the easiest way to use the package is through the built-in interactive dashboard:

```r
library(nestprepper)
launch_app()
```

## Command Line Workflow
For users looking to integrate `nestprepper` into a scripted pipeline:

```r
library(nestprepper)

# 1. Prepare/Clean your data
clean_data <- prep_nesting_data(my_raw_data)

# 2. Convert Nests to Females
abundance <- calculate_abundance(clean_data, clutch_freq = 5.5, remig_int = 3.0)

# 3. Model Trends (Bayesian State-Space)
results <- run_turtle_model(abundance, iter = 10000, parallel = FALSE)

# 4. Visualize Status
plot_turtle_status(results, abundance, species_name = "Leatherback")
```

## Understanding the Output

| Technical Term | Plain English Translation |
| :--- | :--- |
| **Annual Nesters** | The estimated number of individual females that nested this season. |
| **95% Confidence Range** | The "Window of Truth." We are 95% certain the real population size falls within this area. |
| **Growth Rate (%)** | The population "Speedometer." A positive percentage suggests recovery. |
| **Model Precision (MCMC)** | The number of internal "checks" the model runs to ensure the population estimate is stable and reliable. |

## Citations
If you use `nestprepper` for a DPS assessment or publication, please cite:

> Ortega, A. (2026). nestprepper: An R package for sea turtle nesting data integration and trend modeling. Available at: https://github.com/AnnaOrtega-NOAA/nestprepper

The underlying Bayesian state-space methodology for this package is based on:

> Martin SL, Siders Z, Eguchi T, Langseth B, Yau A, Baker J, Ahrens R, Jones TT. 2020. Assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions in the Hawaii-based shallow-set longline fishery. U.S. Dept. of Commerce, NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-95, 183 p.doi:10.25923/ydp1-f891

> Martin SL, Siders Z, Eguchi T, Langseth B, Yau A, Baker J, Ahrens R, Jones TT. 2020. Update to assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions: inclusion of the Hawaii-based deep-set and American Samoa-based longline fisheries. U.S. Dept. of Commerce, NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-101, 67 p. doi:10.25923/pnf2-2q77

## Contributing
We welcome contributions! To ensure scientific integrity, major changes to the modeling logic may undergo internal peer review by internal NOAA staff before being merged. Please open an Issue to discuss proposed changes.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
