# ==========================================
# 1. LOAD PACKAGES
# ==========================================
library(lme4)
library(lmerTest)
library(tidyverse)
library(emmeans)     # Replaces the deprecated 'lsmeans' package
library(performance)
library(report)
library(ggplot2)
library(kableExtra)
library(sjPlot)
library(jtools)
library(insight)

# ==========================================
# 2. DATA LOADING & PREP
# ==========================================
mydata <- read.csv("regionBestfitdata.csv") %>% 
  mutate(
    stage = as.factor(stage),
    region = as.factor(region),
    genotype = as.factor(genotype)
  )

# ==========================================
# 3. MODELING
# ==========================================
# Primary model: Gamma predicted by the interaction of eta and genotype
mod <- lmer(gamma ~ eta * genotype + (1 | fish_id), data = mydata)

# Model Summary & ANOVA
summary(mod)
anova(mod)

# Optional: If you want to compare models, ensure both are defined first:
# mod1 <- glmer(log(gamma) ~ eta * stage * region * genotype + (1|fish_id), data = mydata)
# sidecomp <- compare_performance(mod1, mod)
# export_table(sidecomp, format = "html")

# ==========================================
# 4. MODEL DIAGNOSTICS
# ==========================================
# Create a dataframe for diagnostics
tdat <- data.frame(
  predicted = predict(mod), 
  residual = residuals(mod)
)

# Plot 1: Residuals vs Predicted
ggplot(tdat, aes(x = predicted, y = residual)) + 
  geom_point(alpha = 0.6) + 
  geom_hline(yintercept = 0, lty = 3, color = "red") +
  theme_minimal() +
  labs(title = "Residuals vs Predicted")

# Plot 2: Histogram of Residuals
ggplot(tdat, aes(x = residual)) + 
  geom_histogram(bins = 20, color = "black", fill = "lightblue") +
  theme_minimal() +
  labs(title = "Histogram of Residuals")

# Plot 3: QQ Plot
ggplot(tdat, aes(sample = residual)) + 
  stat_qq() + 
  stat_qq_line(color = "red") +
  theme_minimal() +
  labs(title = "QQ Plot of Residuals")

# ==========================================
# 5. POST-HOC ANALYSIS (emmeans / emtrends)
# ==========================================
# Analyze trends (slopes) of gamma over eta, by genotype
m.lst <- emtrends(mod, "genotype", var = "eta", pbkrtest.limit = 11140)
pairs(m.lst)

# (Note: To use emmeans with 'stage' or 'time', those variables 
# MUST be included in your 'mod' formula above.)
# Example if stage is added to the model:
# mod.em <- emmeans(mod, ~ genotype | stage, adjust = "bonferroni", pbkrtest.limit = 11140)
# pairs(mod.em)
# mod.ci <- confint(mod.em)

# ==========================================
# 6. REPORTING
# ==========================================
# Generate English text report of the model
report(mod)

# Export joint tests to HTML
df_tests <- joint_tests(mod, pbkrtest.limit = 11140)
export_table(df_tests, format = "html")


# ==========================================
# 7. ORPHANED / REFERENCE CODE 
# (From previous projects - commented out to prevent errors)
# ==========================================

# --- Tab Models (Undefined variables: modacfw, modpili, modle) ---
# function_test <- function(arg_1) { exp(arg_1) }
# tab_model(modacfw, modpili, modle)
# tab_model(modacfw, modpili, modle, transform = "function_test", auto.label = FALSE)
# plot_summs(modacfw, modpili, modle, robust = TRUE, plot.distributions = TRUE, model.names = c("ACFW", "PILI", "LE"))

# --- Table Generation (Requires 'time' variable, which is not in 'mod') ---
# mod.em %>% merge(mod.ci) %>% 
#   select(time, estimate, SE, df, lower.CL, upper.CL, p.value) %>% 
#   mutate(
#     estimate = exp(estimate),
#     SE = exp(SE),
#     lower.CL = exp(lower.CL),
#     upper.CL = exp(upper.CL)
#   ) %>% 
#   kbl(digits=c(0,3,3,3,3,3,3,3,3)) %>% 
#   kable_styling()