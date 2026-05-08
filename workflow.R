# --- THE NESTPREPPER WORKFLOW ---

# 1. Setup
library(dplyr)
library(ggplot2)
library(jagsUI)
library(tidyr)

# 2. Generate Play Data
raw_data <- simulate_turtle_data(true_trend = -0.04) # Simulate 4% decline

# 3. QAQC (The Gatekeeper)
clean_data <- prep_nesting_data(raw_data)

# 4. Biological Conversion
abundance_data <- calculate_abundance(clean_data, clutch_freq = 5.5, remig_int = 3.2)

# 5. Run the Math (Bayesian Model)
model_results <- run_turtle_model(abundance_data, iter = 5000)

# 6. Generate the Output
final_plot <- plot_turtle_status(model_results, abundance_data, "Leatherback")
