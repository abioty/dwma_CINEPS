# Load required libraries
library(lme4)
library(ggplot2)
library(dplyr)
library(pROC)
library(cowplot)

# Read and prepare data
data <- read.csv('dwma_gca_cp_f.csv')
categorical_vars <- c("twins", "hdp", "ante_steroids", "sex", "dc_mommilk", "hriskses", "globalcatmod", "motor_cp2")
data[categorical_vars] <- lapply(data[categorical_vars], as.factor)
continuous_vars <- c("DWMA", "ga", "pmamri")
data[paste0(continuous_vars, "_scaled")] <- lapply(data[continuous_vars], scale)

# Center variables used in interaction (important for interpretation of main effects)
data$hriskses_c <- as.numeric(data$hriskses) - mean(as.numeric(data$hriskses), na.rm=TRUE)
data$globalcatmod_c <- as.numeric(data$globalcatmod) - mean(as.numeric(data$globalcatmod), na.rm=TRUE)

# Fit GLMM model with centered variables
fit_glmm_model <- function(outcome_var, data) {
  data_clean <- data[complete.cases(data[, c(outcome_var, "DWMA_scaled", "ga_scaled", "pmamri_scaled", 
                                             "hdp", "ante_steroids", "sex", "dc_mommilk", "globalcatmod_c", 
                                             "hriskses_c", "bpdgrade", "twins")]), ]
  formula <- as.formula(paste(outcome_var, "~ DWMA_scaled + ga_scaled + pmamri_scaled +
                     hdp + ante_steroids + sex + dc_mommilk + globalcatmod_c * hriskses_c + bpdgrade + (1|twins)"))
  model <- glmer(formula, data = data_clean, family = binomial, control = glmerControl(optimizer = "Nelder_Mead"))
  list(model = model, data = data_clean)
}

model_results <- fit_glmm_model("motor_cp2", data)
motor_cp2_model <- model_results$model
data_clean <- model_results$data

# Extract fixed effects summary
fixed_effects_summary <- function(model) {
  fe <- fixef(model)
  se <- sqrt(diag(vcov(model)))
  ci <- fe + outer(se, c(-1.96, 1.96))
  data.frame(
    Estimate = fe,
    OddsRatio = exp(fe),
    LowerCI_OR = exp(ci[,1]),
    UpperCI_OR = exp(ci[,2]),
    P_value = 2 * (1 - pnorm(abs(fe / se)))
  )
}

fe_summary <- fixed_effects_summary(motor_cp2_model)
write.csv(fe_summary, "logistic_fixed_effects_motor_cp2_summary.csv", row.names = TRUE)

# Define predictor names
predictor_names <- c(
  "globalcatmod_c:hriskses_c" = "msBA*HRSS",
  "bpdgrade" = "BPD",
  "globalcatmod_c" = "msBA",
  "hriskses_c" = "HRSS",
  "dc_mommilk1" = "MMDD",
  "sex1" = "SEX",
  "ante_steroids1" = "ACS",
  "hdp1" = "HDP",
  "DWMA_scaled" = "DWMA",
  "pmamri_scaled" = "PMA",
  "ga_scaled" = "GA"
)

# ROC Curve Plot
full_model_pred <- predict(motor_cp2_model, type = "response")
roc_curve <- roc(data_clean$motor_cp2, full_model_pred)
coords_result <- coords(roc_curve, "best")

roc_plot <- ggplot(data.frame(specificity = 1 - roc_curve$specificities, sensitivity = roc_curve$sensitivities), 
                   aes(x = specificity, y = sensitivity)) +
  geom_line(color = "darkred", size = 1.5) +
  geom_abline(linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(x = 1 - coords_result[["specificity"]], y = coords_result[["sensitivity"]]),
             aes(x = x, y = y), color = "deepskyblue", size = 3) +
  geom_text(data = data.frame(x = 0.75, y = 0.25),
            aes(x = x, y = y, label = sprintf("AUC: %.2f\nSens: %.2f\nSpec: %.2f", 
                                              auc(roc_curve), coords_result[["sensitivity"]], coords_result[["specificity"]])),
            color = "black", size = 4.5, fontface = "bold", hjust = 1) +
  theme_classic() +
  theme(axis.line = element_line(size = 1.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 24, face = "bold")) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "a") +
  coord_cartesian(expand = FALSE)

# Odds Ratio Plot
or_plot_data <- fe_summary[-1, ]  # Remove intercept
or_plot_data$Variable <- factor(rownames(or_plot_data), levels = names(predictor_names))
levels(or_plot_data$Variable) <- predictor_names[levels(or_plot_data$Variable)]
or_plot_data <- or_plot_data[order(abs(or_plot_data$Estimate), decreasing = TRUE), ]

odds_ratio_plot <- ggplot(or_plot_data, aes(x = reorder(Variable, OddsRatio), y = OddsRatio)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI_OR, ymax = UpperCI_OR), width = 0.2, size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  coord_flip() +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
                labels = c(".1", ".2", ".5", "1", "2", "5", "10")) +
  theme_classic() +
  theme(axis.line = element_line(size = 1.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 24, face = "bold")) +
  labs(x = "Predictors", y = "Odds Ratio", title = "b")

# Calibration Plot
cal_data <- data.frame(
  predicted = full_model_pred,
  observed = as.numeric(as.character(data_clean$motor_cp2))
)

cal_plot <- ggplot(cal_data, aes(x = predicted, y = observed)) +
  geom_smooth(method = "loess", se = TRUE, color = "purple4", fill = "lightblue", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red3", size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid", size = 1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_classic() +
  theme(axis.line = element_line(size = 1.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none",
        plot.title = element_text(size = 24, face = "bold")) +
  labs(x = "Predicted Probability", y = "Observed Probability", title = "c") +
  coord_cartesian(expand = FALSE)

# Add slope and intercept to calibration plot
cal_model <- lm(observed ~ predicted, data = cal_data)
cal_plot <- cal_plot +
  annotate("text", x = 0.25, y = 0.75, 
           label = sprintf("Slope: %.2f\nIntercept: %.2f", 
                           coef(cal_model)[2], coef(cal_model)[1]),
           size = 4, fontface = "bold")

# Combine the three plots
combined_plot <- plot_grid(roc_plot, odds_ratio_plot, cal_plot,
                           ncol = 3, align = "h", axis = "bt",
                           rel_widths = c(1, 1.2, 1))

# Save combined plot
ggsave("Combined_CP_Model_Plots2.tiff", combined_plot, width = 10, height = 4, dpi = 300, compression = "lzw")

# Print summary of fixed effects for verification
print(fe_summary)


 
