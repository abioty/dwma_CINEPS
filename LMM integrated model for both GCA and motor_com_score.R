# Load required libraries
library(lme4)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(grid)
library(MuMIn)
library(scales)
library(plyr)
library(broom.mixed)
library(lmerTest)

# Read the data
data <- read.csv('merged_data_with_globalbrainscore.csv')

# Convert categorical variables to factors
categorical_vars <- c("twins", "hdp", "ante_steroids", "sex", "dc_mommilk", "hriskses", "globalcatmod")
data[categorical_vars] <- lapply(data[categorical_vars], as.factor)

# Scale continuous variables
continuous_vars <- c("DWMA", "ga", "pmamri")
data[paste0(continuous_vars, "_scaled")] <- lapply(data[continuous_vars], scale)

# Center variables used in interaction (important for interpretation of main effects)
# Convert factors to numeric for centering
data$hriskses_c <- as.numeric(data$hriskses) - mean(as.numeric(data$hriskses), na.rm=TRUE)
data$globalcatmod_c <- as.numeric(data$globalcatmod) - mean(as.numeric(data$globalcatmod), na.rm=TRUE)

# Update the model formula in the fit_lmm_model function
fit_lmm_model <- function(outcome_var) {
  na_count <- sum(is.na(data[[outcome_var]]))
  print(paste("Number of NA values in", outcome_var, ":", na_count))
  data_clean <- data[!is.na(data[[outcome_var]]), ]
  
  formula <- as.formula(paste(outcome_var, "~ DWMA_scaled + ga_scaled + pmamri_scaled +
                     hdp + ante_steroids + sex + bpdgrade +
                     dc_mommilk + globalcatmod_c + hriskses_c + globalcatmod_c:hriskses_c + (1|twins)"))
  lmm_model <- lmer(formula, data = data_clean)
  
  # Create summary dataframe
  coef_summary <- summary(lmm_model)$coefficients
  coef_summary <- data.frame(
    term = rownames(coef_summary),
    estimate = coef_summary[, "Estimate"],
    std.error = coef_summary[, "Std. Error"],
    statistic = coef_summary[, "t value"],
    p.value = coef_summary[, "Pr(>|t|)"]
  )
  
  # Calculate confidence intervals for fixed effects only
  ci <- confint(lmm_model, parm = "beta_")
  
  # Add confidence intervals to coef_summary
  coef_summary$conf.low <- ci[,1]
  coef_summary$conf.high <- ci[,2]
  
  print(coef_summary)
  write.csv(coef_summary, paste0("lmm_coefficients_", outcome_var, "_with_CI.csv"), row.names = FALSE)
  print(summary(lmm_model))
  
  r.squaredGLMM_values <- r.squaredGLMM(lmm_model)
  cat("\nMarginal R-squared:", r.squaredGLMM_values[1], "\n")
  cat("Conditional R-squared:", r.squaredGLMM_values[2], "\n")
  
  data_clean$predicted_values <- fitted(lmm_model)
  
  return(list(model = lmm_model, coef_summary = coef_summary, 
              r_squared = r.squaredGLMM_values, data = data_clean))
}

# Fit models for both outcomes
gca_model <- fit_lmm_model("a_das_gca_stnd")
motor_model <- fit_lmm_model("motor_comp_score")

# Function to create a nicely formatted summary table
create_summary_table <- function(model, outcome_var) {
  coef_summary <- model$coef_summary
  
  # Round numeric columns to 3 decimal places
  coef_summary <- coef_summary %>%
    mutate_if(is.numeric, round, 3)
  
  # Create a column with estimate and CI
  coef_summary$estimate_ci <- sprintf("%.3f (%.3f, %.3f)", 
                                      coef_summary$estimate, 
                                      coef_summary$conf.low, 
                                      coef_summary$conf.high)
  
  # Select columns and rename them manually
  summary_table <- data.frame(
    Term = coef_summary$term,
    `Estimate (95% CI)` = coef_summary$estimate_ci,
    `P-value` = coef_summary$p.value
  )
  
  # Write to CSV
  write.csv(summary_table, paste0("summary_table_", outcome_var, ".csv"), row.names = FALSE)
  
  return(summary_table)
}

# Custom theme for plots
theme_clean <- function(base_size = 24) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 5),
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 40, hjust = -0.05, vjust = 2),
      axis.text = element_text(size = 36, face = "bold"),
      axis.title = element_text(size = 40, face = "bold")
    )
}

# Function to create scatter plot
create_scatter_plot <- function(data, x, y, r_squared, title) {
  lm_fit <- lm(data[[y]] ~ data[[x]])
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  # Determine the range for both axes
  min_val <- min(min(data[[x]]), min(data[[y]]))
  max_val <- max(max(data[[x]]), max(data[[y]]))
  
  # Note: Swapped x and y in the aes() mapping here
  ggplot(data, aes(x = .data[[y]], y = .data[[x]])) +
    geom_point(alpha = 0.7, size = 8) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "aquamarine2", size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +  # Add perfect calibration line
    # Swapped the labels as well
    labs(x = "Actual Values", y = "Predicted Values") +
    theme_clean() +
    theme(
      axis.line = element_line(size = 4),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      aspect.ratio = 1  # Force square aspect ratio
    ) +
    # Updated annotation to reflect correct relationship
    annotate("text", x = min_val, y = max_val, 
             label = sprintf("R²m = %.2f\nR²c = %.2f\nSlope = %.2f\nIntercept = %.2f", 
                             r_squared[1], r_squared[2], slope, intercept), 
             hjust = 0, vjust = 1, size = 12, fontface = "bold") +
    coord_fixed(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +  # Set same limits for both axes
    ggtitle(title)
}

# Function to create fixed effects plot
create_fixed_effects_plot <- function(coef_summary, title) {
  # Map variable names to their abbreviated versions
  predictor_names <- c(
    "pmamri_scaled" = "PMA", 
    "DWMA_scaled" = "DWMA",
    "ga_scaled" = "GA",
    "hdp1" = "HDP", 
    "ante_steroids1" = "ACS", 
    "sex1" = "SEX",
    "dc_mommilk1" = "MMDD", 
    "hriskses_c" = "HRSS", 
    "globalcatmod_c" = "msBA",
    "globalcatmod_c:hriskses_c" = "msBA*HRSS", 
    "interaction_term" = "msBA*HRSS",
    "bpdgrade" = "BPD"
  )
  
  fe_estimates_plot <- coef_summary %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term = plyr::revalue(term, replace = predictor_names)
    ) %>%
    arrange(estimate) %>%
    mutate(term = factor(term, levels = term))
  
  ggplot(fe_estimates_plot, aes(x = term, y = estimate)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  width = 0.5, size = 3, color = "black") +  # Thick error bars
    geom_point(size = 12, color = "black") +  # Large points
    coord_flip() +
    theme_clean() +
    labs(x = "Fixed Effects", y = "Estimate") +
    theme(
      axis.line = element_line(size = 4),
      axis.text.y = element_text(size = 34),
      axis.text.x = element_text(size = 30),
      axis.title = element_text(size = 36, face = "bold"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.title = element_text(size = 40, face = "bold")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 2) +
    ggtitle(title)
}

# Create plots
panel_a <- create_scatter_plot(gca_model$data, "predicted_values", "a_das_gca_stnd", gca_model$r_squared, "A")
panel_b <- create_fixed_effects_plot(gca_model$coef_summary, "B")
panel_c <- create_scatter_plot(motor_model$data, "predicted_values", "motor_comp_score", motor_model$r_squared, "C")
panel_d <- create_fixed_effects_plot(motor_model$coef_summary, "D")

# Combine plots
combined_plot <- grid.arrange(panel_a, panel_b, panel_c, panel_d, 
                              ncol = 2, nrow = 2, 
                              widths = c(1, 1.2), heights = c(1, 1))

# Save combined plot
ggsave("GCA_and_Motor_models2.tiff", 
       plot = combined_plot, 
       device = "tiff", 
       width = 24, height = 24,
       dpi = 300, 
       compression = "lzw")

##Additional codes to obtain uni-variate results:
# Function to fit LMM model with only DWMA and extract results
fit_dwma_only_model <- function(outcome_var, data) {
  data_clean <- data[!is.na(data[[outcome_var]]), ]
  
  formula <- as.formula(paste(outcome_var, "~ DWMA + (1|twins)"))
  lmm_model <- lmerTest::lmer(formula, data = data_clean)
  
  dwma_result <- summary(lmm_model)$coefficients["DWMA", ]
  ci <- confint(lmm_model)["DWMA", ]
  
  estimate <- dwma_result["Estimate"]
  p_value <- dwma_result["Pr(>|t|)"]
  
  result_string <- sprintf("β = %.2f; 95%% CI: %.2f to %.2f; p %s",
                           estimate,
                           ci[1],
                           ci[2],
                           ifelse(is.na(p_value) || p_value < 0.001, 
                                  "< 0.001", 
                                  sprintf("= %.3f", p_value)))
  
  return(result_string)
}

# Fit models and print results
outcomes <- c("a_das_gca_stnd", "motor_comp_score")
names(outcomes) <- c("GCA", "Motor")

for (name in names(outcomes)) {
  result <- fit_dwma_only_model(outcomes[name], data)
  print(paste(name, ":", result))
}

# Function to extract formatted results for all predictors in a model
extract_all_predictors_formatted <- function(model, outcome_name) {
  # Get model summary
  model_summary <- summary(model)$coefficients
  
  # Get confidence intervals
  ci <- confint(model)
  
  # Create empty data frame for results
  results <- data.frame(
    Predictor = rownames(model_summary),
    Beta = numeric(nrow(model_summary)),
    CI_Lower = numeric(nrow(model_summary)),
    CI_Upper = numeric(nrow(model_summary)),
    P_Value = numeric(nrow(model_summary))
  )
  
  # Fill in values
  for (i in 1:nrow(model_summary)) {
    var_name <- rownames(model_summary)[i]
    results$Predictor[i] <- var_name
    results$Beta[i] <- model_summary[var_name, "Estimate"]
    results$CI_Lower[i] <- ci[var_name, 1]
    results$CI_Upper[i] <- ci[var_name, 2]
    results$P_Value[i] <- model_summary[var_name, "Pr(>|t|)"]
  }
  
  # Replace predictor names with more readable names
  name_mapping <- c(
    "(Intercept)" = "Intercept",
    "DWMA_scaled" = "DWMA",
    "ga_scaled" = "GA",
    "pmamri_scaled" = "PMA",
    "hdp1" = "HDP",
    "ante_steroids1" = "ACS",
    "sex1" = "Sex",
    "dc_mommilk1" = "MMDD", 
    "hriskses_c" = "HRSS",
    "globalcatmod_c" = "msBA",
    "globalcatmod_c:hriskses_c" = "msBA*HRSS",
    "interaction_term" = "msBA*HRSS",
    "bpdgrade" = "BPD"
  )
  
  # Apply name mapping
  results$Predictor <- sapply(results$Predictor, function(x) {
    if (x %in% names(name_mapping)) name_mapping[x] else x
  })
  
  # Format p-values
  results$P_Value_Formatted <- ifelse(
    results$P_Value < 0.001, 
    "<0.001", 
    sprintf("%.3f", results$P_Value)
  )
  
  # Round numeric values to 1 decimal place
  results$Beta <- round(results$Beta, 1)
  results$CI_Lower <- round(results$CI_Lower, 1)
  results$CI_Upper <- round(results$CI_Upper, 1)
  
  # Format exactly as in the table
  formatted_table <- data.frame(
    Predictors = results$Predictor,
    Beta = results$Beta,
    CI = sprintf("%.1f, %.1f", results$CI_Lower, results$CI_Upper),
    P = results$P_Value_Formatted
  )
  
  # Set row names to NULL
  rownames(formatted_table) <- NULL
  
  # Remove intercept from results
  formatted_table <- formatted_table[formatted_table$Predictors != "Intercept", ]
  
  # Write to CSV
  write.csv(formatted_table, paste0("formatted_table_", outcome_name, ".csv"), row.names = FALSE)
  
  return(formatted_table)
}

# Extract results for both models
gca_table <- extract_all_predictors_formatted(gca_model$model, "cognitive")
motor_table <- extract_all_predictors_formatted(motor_model$model, "motor")

# Print tables
print(gca_table)
print(motor_table)


#Sensitivity analysis (after removing subjects with severe brain abnormality score)

# Create filtered dataset excluding subjects with moderate-severe brain abnormality
data_filtered <- data[data$globalbrainscore2 <= 12, ]

# Print the number of subjects included/excluded
cat("Original dataset:", nrow(data), "subjects\n")
cat("Filtered dataset:", nrow(data_filtered), "subjects\n")
cat("Excluded", nrow(data) - nrow(data_filtered), "subjects with moderate-severe brain abnormality (score > 7)\n")

# Update the model formula in the fit_lmm_model function
fit_lmm_model <- function(outcome_var, dataset) {
  na_count <- sum(is.na(dataset[[outcome_var]]))
  print(paste("Number of NA values in", outcome_var, ":", na_count))
  data_clean <- dataset[!is.na(dataset[[outcome_var]]), ]
  
  formula <- as.formula(paste(outcome_var, "~ DWMA_scaled + ga_scaled + pmamri_scaled +
                     hdp + ante_steroids + sex + bpdgrade +
                     dc_mommilk + globalcatmod_c + hriskses_c + globalcatmod_c:hriskses_c + (1|twins)"))
  lmm_model <- lmer(formula, data = data_clean)
  
  # Create summary dataframe
  coef_summary <- summary(lmm_model)$coefficients
  coef_summary <- data.frame(
    term = rownames(coef_summary),
    estimate = coef_summary[, "Estimate"],
    std.error = coef_summary[, "Std. Error"],
    statistic = coef_summary[, "t value"],
    p.value = coef_summary[, "Pr(>|t|)"]
  )
  
  # Calculate confidence intervals for fixed effects only
  ci <- confint(lmm_model, parm = "beta_")
  
  # Add confidence intervals to coef_summary
  coef_summary$conf.low <- ci[,1]
  coef_summary$conf.high <- ci[,2]
  
  return(list(model = lmm_model, coef_summary = coef_summary))
}

# Fit models for both outcomes with original dataset
gca_model <- fit_lmm_model("a_das_gca_stnd", data)
motor_model <- fit_lmm_model("motor_comp_score", data)

# Fit models with filtered dataset (excluding moderate-severe brain abnormality)
gca_model_filtered <- fit_lmm_model("a_das_gca_stnd", data_filtered)
motor_model_filtered <- fit_lmm_model("motor_comp_score", data_filtered)

# Extract DWMA results from original models
extract_dwma_results <- function(model_result) {
  coef_row <- model_result$coef_summary[model_result$coef_summary$term == "DWMA_scaled", ]
  
  estimate <- round(coef_row$estimate, 1)
  ci_low <- round(coef_row$conf.low, 1)
  ci_high <- round(coef_row$conf.high, 1)
  p_value <- coef_row$p.value
  
  # Format p-value
  if (p_value < 0.001) {
    p_formatted <- "<0.001"
  } else {
    p_formatted <- sprintf("%.3f", p_value)
  }
  
  result <- sprintf("β=%.1f; 95%% CI, %.1f to %.1f; p=%s", 
                    estimate, ci_low, ci_high, p_formatted)
  
  return(result)
}

# Extract DWMA results
gca_result_original <- extract_dwma_results(gca_model)
motor_result_original <- extract_dwma_results(motor_model)
gca_result_filtered <- extract_dwma_results(gca_model_filtered)
motor_result_filtered <- extract_dwma_results(motor_model_filtered)

# Print results for comparison
cat("DWMA effects on Cognitive outcome (original dataset):\n", gca_result_original, "\n\n")
cat("DWMA effects on Cognitive outcome (excluding moderate-severe brain abnormality):\n", gca_result_filtered, "\n\n")
cat("DWMA effects on Motor outcome (original dataset):\n", motor_result_original, "\n\n")
cat("DWMA effects on Motor outcome (excluding moderate-severe brain abnormality):\n", motor_result_filtered, "\n\n")

# For CP model (if you have this variable)
# Function to fit GLMM model (for CP outcome)
fit_glmm_model <- function(outcome_var, dataset) {
  dataset_clean <- dataset[complete.cases(dataset[, c(outcome_var, "DWMA_scaled", "ga_scaled", "pmamri_scaled", 
                                                      "hdp", "ante_steroids", "sex", "dc_mommilk", 
                                                      "globalcatmod_c", "hriskses_c", "bpdgrade", "twins")]), ]
  
  formula <- as.formula(paste(outcome_var, "~ DWMA_scaled + ga_scaled + pmamri_scaled +
                     hdp + ante_steroids + sex + dc_mommilk + 
                     globalcatmod_c * hriskses_c + bpdgrade + (1|twins)"))
  
  model <- glmer(formula, data = dataset_clean, family = binomial, 
                 control = glmerControl(optimizer = "Nelder_Mead"))
  
  # Extract summary
  model_summary <- summary(model)$coefficients
  
  # Extract DWMA coefficient
  dwma_coef <- model_summary["DWMA_scaled", ]
  
  # Calculate CI
  ci <- confint(model)["DWMA_scaled", ]
  
  # Convert to OR
  or <- exp(dwma_coef["Estimate"])
  ci_or <- exp(ci)
  
  # Format result
  p_value <- dwma_coef["Pr(>|z|)"]
  p_formatted <- ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  
  result <- sprintf("aOR=%.1f; 95%% CI, %.1f to %.1f; p=%s", 
                    round(or, 1), round(ci_or[1], 1), round(ci_or[2], 1), p_formatted)
  
  return(list(model = model, result = result))
}
