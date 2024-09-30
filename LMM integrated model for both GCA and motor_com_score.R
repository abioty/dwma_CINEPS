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
data <- read.csv('dwma_gca_cp_f.csv')

# Convert categorical variables to factors
categorical_vars <- c("twins", "hdp", "ante_steroids", "sex", "dc_mommilk", "hriskses", "globalcatmod")
data[categorical_vars] <- lapply(data[categorical_vars], as.factor)

# Scale continuous variables
continuous_vars <- c("DWMA", "ga", "pmamri")
data[paste0(continuous_vars, "_scaled")] <- lapply(data[continuous_vars], scale)

# Function to fit LMM model and extract results
fit_lmm_model <- function(outcome_var) {
  na_count <- sum(is.na(data[[outcome_var]]))
  print(paste("Number of NA values in", outcome_var, ":", na_count))
  data_clean <- data[!is.na(data[[outcome_var]]), ]
  
  formula <- as.formula(paste(outcome_var, "~ DWMA_scaled + ga_scaled + pmamri_scaled +
                     hdp + ante_steroids + sex + 
                     dc_mommilk + globalcatmod + hriskses + bpdgrade + (1|twins)"))
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
  
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.7, size = 8) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "aquamarine2", size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +  # Add perfect calibration line
    labs(x = "Predicted Values", y = "Actual Values") +
    theme_clean() +
    theme(
      axis.line = element_line(size = 4),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      aspect.ratio = 1  # Force square aspect ratio
    ) +
    annotate("text", x = min_val, y = max_val, 
             label = sprintf("R²m = %.2f\nR²c = %.2f\nSlope = %.2f\nIntercept = %.2f", 
                             r_squared[1], r_squared[2], slope, intercept), 
             hjust = 0, vjust = 1, size = 12, fontface = "bold") +
    coord_fixed(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +  # Set same limits for both axes
    ggtitle(title)
}

# Function to create fixed effects plot
create_fixed_effects_plot <- function(coef_summary, title) {
  predictor_names <- c(
    "pmamri_scaled" = "PMA (scaled)", "DWMA_scaled" = "DWMA (scaled)",
    "hdp1" = "HDP", "ante_steroids1" = "ACS", "sex1" = "SEX",
    "dc_mommilk1" = "MMDD", "hriskses1" = "HRSS", "globalcatmod1" = "msBAS",
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

 