# Classification Tree: Presentation Notes

## Method
- CART classification trees fit to the same 21 health indicators used by the logistic/KNN analysis.
- BMI remains continuous; documented BRFSS binary/ordinal numeric codes are retained to match the logistic/KNN baseline.
- Each dataset uses a stratified 80/20 split with seed 42.
- The complexity parameter is chosen by minimum 10-fold cross-validation error and the tree is pruned before test evaluation.

## Test-set results

### Binary tasks
           dataset accuracy precision recall specificity    f1   auc
   Balanced binary    0.719     0.706  0.752       0.686 0.728 0.770
 Imbalanced binary    0.840     0.575  0.128       0.981 0.209 0.718

### Three-class task
                dataset accuracy macro_precision macro_recall macro_specificity
 Three-class imbalanced    0.839           0.463        0.375             0.704
 macro_f1 macro_auc_ovr
    0.379         0.673

## Interpretation
- Report recall, precision, F1, and AUC alongside accuracy; accuracy is inflated when the no-diabetes class dominates.
- The three-class per-class table is the key evidence for whether the model recognizes prediabetes rather than only the majority class.
- Variable importance describes predictive splitting value, not causal effect.

## Top five splitting variables: imbalanced binary tree
- HighBP: 3672.7
- GenHlth: 2492.1
- HighChol: 880.6
- Age: 647.5
- DiffWalk: 553.2
