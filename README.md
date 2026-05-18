# nestprepper (`marss-integration` branch)

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)

> **Branch Note (`marss-integration`):** This active development branch introduces a robust **Multi-Engine Comparative Suite**. It pairs the default Bayesian assessment framework with a frequentist companion engine (MARSS) to validate population trend assumptions and diagnose localized tracking anomalies.

The goal of `nestprepper` is to streamline the transition from messy, real-world sea turtle nesting counts to robust biological abundance estimates and verified population trends. This package supports NOAA Open Science initiatives by providing a transparent, reproducible pipeline for formal population status assessments.

## Key Features

* **Smart Mapping & Aggregation**: Automatically maps "Wide" or "Long" data layouts and dynamically aggregates monthly nesting records into standard annual sums.
* **Built-in QAQC**: Interactively evaluates biological anomalies and provides clear options for handling missing monitoring seasons.
* **Multi-Engine Framework**: Estimates population trends using Bayesian State-Space models (via JAGS) and cross-validates results using frequentist maximum-likelihood frameworks (via MARSS).
* **Structural Diagnostics Sandbox**: Allows managers to unconstrain site-specific parameters to test if individual beaches are moving in harmony with the region or dynamically decoupling.
* **Interactive Tooltip Architecture**: Embedded contextual tooltips translate complex variance parameters into concrete fieldwork and environmental feedback.

## Data Input Requirements

`nestprepper` is optimized for regional and localized population status reviews. 

### Standardized Data Structure
Your input CSV can contain either monthly tracking logs or pre-aggregated annual sums. The pipeline is optimized to identify:

* **Temporal Mapping**: Automatically flags columns named `Year`, `Season`, `Year_begin`, or `Month`.
* **Count Mapping**: Users can isolate specific columns containing survey counts (e.g., `JM_Nests` or `Nests`).
* **Multi-Site Arrays**: Supports multiple beach columns simultaneously to analyze complex regional complexes.

```r
# Ideal structure for nestprepper:
# Year  | Month    | Beach_A | Beach_B
# 2024  | May      | 45      | 12
# 2024  | June     | 112     | 38
# 2025  | May      | NA      | 14  (Handles missing monitoring slots)
```

### Constraints & Assumptions
* **Time-Series Integrity**: While the state-space engines handle missing survey years (`NA` values), a consistent time-series sequence is required to properly isolate trends from background noise.

## NOAA Disclaimer
This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce. All NOAA GitHub project code is provided "as is" and the user assumes the entire risk as to its quality and performance. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the Department of Commerce.

## Installation
You can install this specific integration version of `nestprepper` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("AnnaOrtega-NOAA/nestprepper", ref = "marss-integration")
```

## Quick Start (Shiny Dashboard)
Launch the multi-engine comparison sandbox interface directly from your console:

```r
library(nestprepper)
launch_app()
```

## Command Line Workflow
For scripted, reproducible assessment pipelines:

```r
library(nestprepper)

# 1. Prepare and filter raw multi-beach data
clean_data <- prep_nesting_data(my_raw_data)

# 2. Convert raw counts to estimated adult female abundance
abundance <- calculate_abundance(clean_data, clutch_freq = 5.5, remig_int = 3.06)

# 3. Model regional trends and cross-validate across engines
jags_results  <- run_turtle_model(abundance, iter = 50000)
marss_results <- run_marss_models(abundance)

# 4. Generate management-ready status trajectories
plot_turtle_status(jags_results, abundance)
```

## Understanding the Management Diagnostics

| Technical Term | Management Context & Interpretation |
| :--- | :--- |
| **Annual Population Trend ($U$)** | The long-term trajectory of the population. A positive percentage indicates a recovering population trend. |
| **Environmental Fluctuations ($Q$)** | True year-to-year shifts in population size driven by changing environmental conditions (e.g., climate cycles, foraging availability). |
| **Monitoring & Survey Noise ($R$)** | Data variation introduced by human tracking limits (e.g., missed survey days, weather dropouts, or shifting field effort). |
| **Shared Regional Assumption Mismatch** | Occurs when individual beaches decouple from the regional average. This flags that the sub-populations are behaving independently and forcing a single regional trend line may mask localized risks. |

## Citations
If you utilize `nestprepper` for a Distinct Population Segment (DPS) assessment or public report, please cite:

> Ortega, A. (2026). nestprepper: An R package for sea turtle nesting data integration and trend modeling. Available at: https://github.com/AnnaOrtega-NOAA/nestprepper

The underlying Bayesian and frequentist state-space methodology for this package is based on:

> Martin SL, Siders Z, Eguchi T, Langseth B, Yau A, Baker J, Ahrens R, Jones TT. 2020. Assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions in the Hawaii-based shallow-set longline fishery. U.S. Dept. of Commerce, NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-95, 183 p. doi:10.25923/ydp1-f891

> Martin SL, Siders Z, Eguchi T, Langseth B, Yau A, Baker J, Ahrens R, Jones TT. 2020. Update to assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions: inclusion of the Hawaii-based deep-set and American Samoa-based longline fisheries. U.S. Dept. of Commerce, NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-101, 67 p. doi:10.25923/pnf2-2q77

## Contributing
We welcome contributions! To ensure scientific integrity, major changes to the modeling logic may undergo internal peer review by NOAA senior assessment staff before being merged into the production branch. Please open an Issue to discuss structural proposals.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
