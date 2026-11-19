# Diabetes Prediction Using Health Indicators

This project uses 2021 BRFSS health survey data to study whether common health indicators can predict diabetes status.

## Goals

- Compare several classification models.
- Understand how class imbalance affects results.
- Identify health indicators that are consistently useful for prediction.

## Data

The data include 21 health indicators, such as BMI, age, general health, high blood pressure, high cholesterol, physical activity, income and education.

We used three versions of the data:

- Three-class: no diabetes, prediabetes, diabetes
- Imbalanced binary dataset
- Balanced binary dataset

## Models

We tested:

- Logistic regression
- One-versus-rest logistic regression
- KNN
- Classification tree
- Ridge and lasso regularization
- Random forest
- Boosting

## Main Findings

- Accuracy alone was misleading because the data were imbalanced.
- We focused on recall, F1 score, ROC-AUC, and confusion matrices.
- Lowering the prediction threshold to `0.20` improved recall to about `0.68` and F1 score to about `0.48`.
- Ridge performed similarly to logistic regression.
- Lasso retained 14 of 21 predictors while keeping a similar ROC-AUC.
- Important predictors included BMI, age, general health, high blood pressure, high cholesterol, and income.

These results show predictive associations, not causal relationships.

## Files

- `diabetes_eda.ipynb` - exploratory data analysis
- `diabetes_ovr_logistic.ipynb` - one-versus-rest logistic regression
- `prediction_analysis(logistic and KNN).Rmd` - logistic regression and KNN
- `tree_classification_analysis.R` - classification tree
- `random_forest_analysis.Rmd` - random forest
- `ridge_lasso_regularization.R` - ridge and lasso analysis
- `tree_results/` - tree figures and tables
- `diabetes_Filtered_CV.ipynb` - XGboost and RF for filtered data

## Team

Lasya Banka, Bailey Allison, Jiajun Guo, Luis Fuentes, and Maggio Laquidara
