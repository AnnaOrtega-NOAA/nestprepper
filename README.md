[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
The goal of `nestprepper` is to streamline the transition from messy, real-world sea turtle nesting counts to robust biological abundance estimates and Bayesian trend analyses. This package supports NOAA Open Science initiatives by providing a transparent, reproducible pipeline for population status assessments.

## Key Features

* **Smart Mapping & Multi-Site Scaling**: Automatically handles both "Wide" format (multiple beach columns) and "Long" format data. It dynamically accommodates an infinite number of sub-beaches/time-series per island without hardcoding.
* **Fourier Monthly Imputation**: Houses an integrated Bayesian AR1 Fourier series engine (6-month or 12-month periods) to dynamically fill monthly observation gaps and aggregate them to the annual level.
* **Built-in QA/QC Controls**: Interactively flags extreme biological spikes (>3 SD outliers) and automatically scrubs backward-casted model hallucinations prior to the first true year of empirical beach monitoring.
* **Unified Take Integration**: Replicates historical fishery interaction adjustments by dynamically allocating Adult Nester Equivalents (ANE) to sub-sites proportional to their annual nesting density.
* **Bayesian State-Space Engine**: Fits exact Boyd et al. (2016) state-space models (`singleUQ` architecture) via JAGS to separate true population trends ($U$) from observational noise ($R$) and environmental process variance ($Q$).

## Data Input Requirements

`nestprepper` is built to handle the biological noise of standard population assessments. The app natively ingests two data configurations:

### 1. Raw Monthly Nest Records (Highly Recommended)
Your CSV can contain raw, monthly observations spanning multiple sub-beaches. The integrated pipeline automatically identifies missing survey months and handles them mathematically.
* **Required Mappings**: A `Year` column, a `Month` column, and one or more distinct `Site_Count` columns (e.g., `A_Nests`, `B_Nests`, `Nests_C`).

### 2. Standardized Annual Totals
If your data have already been aggregated, your input CSV can follow a clean "one-row-per-year" structure.
* **Required Mappings**: A `Year` column and selected `Site_Count` columns. 

```r
# Ideal multi-site monthly structure for nestprepper:
# Year | Month | A_Nests | B_Nests
# 2005 | 5     | 184      | 12
# 2005 | 6     | 379      | NA  (Internal gap handled by Fourier imputation)
# 2005 | 7     | 210      | 45
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
For most users, the easiest way to interact with the package is through the built-in, reactive Shiny user interface:

```r
library(nestprepper)
# Run the local development app layout
runApp('inst/shiny-examples/nestprepper-app')
```

## Command Line Workflow
For power users looking to integrate `nestprepper` into a scripted, reproducible markdown pipeline:

```r
library(nestprepper)

# 1. Prepare, clean, and impute raw monthly data (Autodetects 6 vs 12 month peaks)
annual_data <- aggregate_monthly_to_annual(my_raw_monthly_data)

# 2. Convert raw nest counts to Total Adult Female Abundance
abundance <- calculate_abundance(annual_data, clutch_freq = 5.5, remig_int = 3.06)

# 3. Optional: Inject historical bycatch distributions (ANE)
# intrinsic_abundance <- inject_historical_take(abundance, take_data)

# 4. Model Trends using the Tech Memo State-Space configuration
results <- run_turtle_model(abundance, iter = 150000, burnin = 50000, thin = 10)

# 5. Visualize Status and Posterior Densities
plot_turtle_status(results, abundance)
```

## Understanding the Output

| Technical Term | Plain English Translation |
| :--- | :--- |
| **Annual Nesters** | The estimated number of individual females that nested this season across all selected beaches combined. |
| **Total Females** | The regional population size estimate, scaling annual nesters out by the species' remigration interval. |
| **Growth Rate (U)** | A negative percentage indicates a long-term decline ($U < 0$); a positive percentage indicates active recovery. |
| **Model Precision (MCMC)** | The number of internal Markov Chain Monte Carlo checks the JAGS engine runs to guarantee the posterior distributions are perfectly converged and stable. |

## Citations
If you use `nestprepper` for a Distinct Population Segment (DPS) assessment or public report, please cite:

> Ortega, A. (2026). nestprepper: An R package for sea turtle nesting data integration and trend modeling. Available at: https://github.com/AnnaOrtega-NOAA/nestprepper

The underlying Bayesian state-space and imputation methodologies for this package are based on:

> Martin SL, Siders Z, Eguchi T, Langseth B, Yau A, Baker J, Ahrens R, Jones TT. 2020. Assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions in the Hawaii-based shallow-set longline fishery. U.S. Dept. of Commerce, NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-95, 183 p. doi:10.25923/ydp1-f891

> Martin SL, Siders Z, Eguchi T, Langseth B, Yau A, Baker J, Ahrens R, Jones TT. 2020. Update to assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions: inclusion of the Hawaii-based deep-set and American Samoa-based longline fisheries. U.S. Dept. of Commerce, NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-101, 67 p. doi:10.25923/pnf2-2q77

## Contributing
We welcome contributions! To ensure scientific integrity, major adjustments to the underlying modeling logic or priors may undergo internal statistical peer review by NOAA Science Center staff before being merged. Please open an Issue to discuss proposed structural updates.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
```
