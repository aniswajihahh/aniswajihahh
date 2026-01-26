#### LOGISTIC REGRESSION IN MEDICAL STUDIES ####

# Load dataset
heart <- read.csv("HEART_FAILURE.csv", header = TRUE)

# Explore structure and summary of the dataset
str(heart)       # Check variable types and structure
summary(heart)   # Summary statistics (mean, median, min, max, etc.)

# -------------------------------
# 1. Set categorical variables as factor
# -------------------------------
# In medical logistic regression, categorical variables must be factors.
# This ensures R treats them correctly for model fitting.
heart$anaemia <- factor(heart$anaemia)
heart$diabetes <- factor(heart$diabetes)
heart$high_blood_pressure <- factor(heart$high_blood_pressure)
heart$sex <- factor(heart$sex)
heart$smoking <- factor(heart$smoking)
heart$DEATH_EVENT <- factor(heart$DEATH_EVENT)

# -------------------------------
# 2. Fit Logistic Regression Models
# -------------------------------

# Null model: only intercept, no predictors
# Useful as a baseline to compare more complex models
model_null <- glm(DEATH_EVENT ~ 1, data = heart, family = binomial)
summary(model_null)

# Full model: includes all available predictors
# Helps identify which variables may contribute to DEATH_EVENT
model_full <- glm(DEATH_EVENT ~ ., data=heart, family=binomial)
summary(model_full)

# Reduced model: only selected predictors based on clinical relevance
# Focus on age, ejection fraction, serum creatinine, and time
model_reduced <- glm(DEATH_EVENT ~ age + ejection_fraction + serum_creatinine + time, 
                     data = heart, family = binomial)
summary(model_reduced)

# -------------------------------
# 3. Model Comparison and Fit
# -------------------------------

# a) Likelihood Ratio Test (LRT) for comparing models
# Tests if adding predictors significantly improves the model fit
anova(model_null, model_full, model_reduced, test = "Chisq")

# b) Akaike Information Criterion (AIC)
# Lower AIC indicates a better-fitting model (penalizes overfitting)
AIC(model_null, model_full, model_reduced)

# c) Compare reduced model to null model specifically
anova(model_reduced, test = "Chisq")

# -------------------------------
# 4. Multicollinearity Check
# -------------------------------
library(car)
# Variance Inflation Factor (VIF) detects multicollinearity
# VIF > 5–10 indicates potential issues
vif(model_reduced)

# -------------------------------
# 5. Goodness-of-fit test
# -------------------------------
library(ResourceSelection)
# Hosmer-Lemeshow test checks how well predicted probabilities match observed outcomes
hoslem.test(heart$DEATH_EVENT, fitted(model_reduced))

# -------------------------------
# 6. Influence Diagnostics
# -------------------------------
# Identify outliers or influential observations that may affect model
plot(cooks.distance(model_reduced), type="h")  # Cook's distance
plot(hatvalues(model_reduced))                # Leverage points
plot(residuals(model_reduced, type="deviance"))  # Deviance residuals

# -------------------------------
# 7. Model Diagnostics
# -------------------------------
library(pscl)
# Pseudo R-squared measures model explanatory power (not exactly R^2)
pR2(model_reduced)

# ROC curve for model discrimination
library(pROC)
roc_curve <- roc(heart$DEATH_EVENT, fitted(model_reduced))
plot(roc_curve)       # Visualize sensitivity vs 1-specificity
auc(roc_curve)        # Area Under the Curve (AUC) measures predictive accuracy

# -------------------------------
# 8. Odds Ratios and Confidence Intervals
# -------------------------------
# OR provides clinical interpretation of predictor effect
OR <- exp(coef(model_reduced))           # Exponentiate coefficients
CI <- confint(model_reduced)             # 95% Confidence Interval for OR


