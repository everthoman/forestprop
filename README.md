# forestprop

A command-line QSAR tool for building conformal prediction models from SMILES data using Random Forest. Supports both regression (pIC50, Ki, etc.) and binary classification with valid coverage guarantees via conformal prediction, and optional Bayesian hyperparameter optimisation.

---

## Features

- **Conformal prediction** — generates prediction intervals (regression) or prediction sets (classification) with user-specified coverage guarantees (e.g. 90%)
- **Two conformal modes** — cross-conformal (jackknife+, default) or prefit with a dedicated calibration set
- **Molecular featurisation** — Mordred 2D descriptors (~1600 features) and/or RDKit fingerprints (Morgan, RDKit, MACCS, AtomPair)
- **Bayesian HPO** — optional Optuna search over Random Forest hyperparameters using proper nested cross-validation (test set is never touched during HPO)
- **Train/test diagnostics** — reports in-sample and out-of-sample metrics side by side to quantify overfitting
- **Inference mode** — save trained models and apply them to new compound sets
- **Outputs** — scatter plot, per-compound predictions, feature importances, metrics CSV, and (optionally) a serialised model checkpoint

---

## Installation

```bash
conda env create -f environment.yml
conda activate ml
```

The environment uses Python 3.11 with conda-forge packages. RDKit is always installed via conda (not pip) for binary compatibility. Mordred, MAPIE, and Optuna are installed via pip.

**Manual install (pip only):**
```bash
pip install scikit-learn mordred rdkit-pypi mapie optuna tqdm pandas numpy
# If mordred gives issues on Python 3.11+:
pip install mordredcommunity
```

---

## Usage

### Train

```bash
python forestprop.py train -i data.csv --activity_col pIC50
```

The `--confidence_level` is prompted interactively if not supplied.

**Key options:**

| Option | Default | Description |
|---|---|---|
| `--activity_col` | `activity` | Column name for the target variable |
| `--smiles_col` | `SMILES` | Column name for SMILES strings |
| `--id_col` | *(auto)* | Column name for compound IDs (optional) |
| `--task` | `regression` | `regression` or `classification` |
| `--threshold` | — | Binarise activity ≥ threshold (classification only) |
| `--confidence_level` | *(prompted)* | Conformal coverage, e.g. `0.9` for 90% |
| `--features` | `mordred rdkit` | Feature types: `mordred`, `rdkit`, or both |
| `--fp_types` | `morgan` | Fingerprint type(s): `morgan`, `rdkit`, `maccs`, `atompair` |
| `--fp_bits` | `2048` | Fingerprint bit length |
| `--cv_folds` | `5` | Cross-conformal folds (0 = prefit mode) |
| `--test_size` | `0.15` | Fraction held out as final test set |
| `--hpo_trials` | `0` | Optuna trials (0 = disabled) |
| `--hpo_cv_folds` | `3` | Inner CV folds for HPO |
| `--output` / `-o` | `results/` | Output directory |
| `--save_model` | off | Serialise model for later inference |
| `--seed` | `42` | Random seed |

**Examples:**

```bash
# Regression with default settings
python forestprop.py train -i compounds.csv --activity_col pIC50 --confidence_level 0.9

# With Bayesian HPO (100 trials)
python forestprop.py train -i compounds.csv --activity_col pIC50 \
    --confidence_level 0.9 --hpo_trials 100 -o results_hpo/

# Binary classification with threshold binarisation
python forestprop.py train -i compounds.csv --task classification \
    --activity_col pIC50 --threshold 7.0 --confidence_level 0.9

# Prefit mode (faster, for larger datasets)
python forestprop.py train -i compounds.csv --cv_folds 0 --confidence_level 0.9

# Save model for inference
python forestprop.py train -i compounds.csv --confidence_level 0.9 \
    --save_model --model_name SMUG1i_pIC50 -o model_output/
```

### Predict

Apply a saved model to new compounds:

```bash
python forestprop.py predict -i new_compounds.csv \
    --load_model model_output/SMUG1i_pIC50_<date>.pkl \
    -o predictions/
```

---

## Outputs

All files are written to the output directory specified with `-o`.

| File | Description |
|---|---|
| `predictions_test.csv` | Per-compound predictions for the held-out test set, including prediction intervals (regression) or prediction sets (classification) |
| `predictions_all.csv` | Predictions for all compounds (train + test), with a `split` column indicating membership |
| `metrics.csv` | Train and test metrics in a single row (R², RMSE, MAE for regression; Accuracy, MCC, ROC-AUC for classification; plus conformal coverage) |
| `feature_importance.csv` | Top-50 features by Random Forest importance |
| `scatter_pred_vs_exp.png` | Predicted vs experimental scatter plot with 90% conformal intervals on test points and train R²/test R² annotated |
| `hpo_trials.csv` | Per-trial HPO results (only with `--hpo_trials > 0`) |
| `hpo_param_importance.csv` | Optuna parameter importances (only with `--hpo_trials > 0`) |
| `*.pkl` | Serialised model checkpoint (only with `--save_model`) |

---

## Architecture

### Conformal prediction modes

**Cross-conformal (default, `--cv_folds k`):**
MAPIE trains *k* Random Forest folds; out-of-fold residuals calibrate the conformal scores. Every compound contributes to both training and calibration. Uses jackknife+ (regression) or score (classification).

**Prefit (`--cv_folds 0`):**
A dedicated holdout calibration set (controlled by `--cal_size`) is used. Simpler and faster for large datasets.

### Nested HPO

When `--hpo_trials N` is set, HPO runs inside a proper nested cross-validation:

```
┌─ Outer split: trainval vs. test (isolated throughout) ───────────────┐
│  ┌─ HPO inner CV (--hpo_cv_folds, default 3) ──────────────────────┐ │
│  │  Optuna minimises CV-RMSE (regression) / CV-log-loss (classif.) │ │
│  │  Searches: n_estimators, max_depth, min_samples_split,          │ │
│  │            min_samples_leaf, max_features, max_samples          │ │
│  └──────────────────────────────────────────────────────────────── ┘ │
│  Best params → final cross-conformal model on full trainval           │
│  Conformal coverage evaluated on outer test set only                  │
└───────────────────────────────────────────────────────────────────────┘
```

This ensures the conformal coverage guarantee is not inflated by hyperparameter selection.

### Train vs. test metrics

Both in-sample (train) and out-of-sample (test) metrics are reported. The train metrics use bare RF predictions (no conformal wrapping — in-sample conformal coverage is not meaningful). The gap between train and test performance is a direct measure of overfitting.

---

## Input format

A CSV file with at minimum a SMILES column and an activity column. An ID column is optional.

```
ID,SMILES,pIC50
CPD001,Cc1ccccc1,6.3
CPD002,c1ccccc1,5.1
...
```

Column names are matched case-insensitively.

---

## Dependencies

| Package | Purpose |
|---|---|
| `scikit-learn` | Random Forest, imputation, cross-validation |
| `mapie` | Conformal prediction wrapper |
| `rdkit` | SMILES parsing and fingerprints |
| `mordred` / `mordredcommunity` | 2D molecular descriptors |
| `optuna` | Bayesian hyperparameter optimisation |
| `pandas`, `numpy` | Data handling |
| `tqdm` | Progress bars |
| `matplotlib` | Scatter plot output |
