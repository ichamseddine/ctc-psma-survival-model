# Cox Survival Model for PSMA-PET and CTC Biomarkers

This repository contains the full analysis pipeline for building and validating a Cox proportional hazards model to stratify prostate cancer patients treated with LuPSMA, based on PSMA-PET imaging and CTC transcriptomic features.

## Overview

The analysis integrates multiple types of biomarkers to predict radiographic progression-free survival (rPFS). It includes:

1. Loading the dataset
2. Removing multicollinear features using variance inflation factor (VIF)
3. Fitting a Cox model using all remaining features
4. Generating a volcano plot of coefficients vs. -log10(p-values)
5. Selecting top features based on significance and importance
6. Re-fitting the model using selected features
7. Computing risk scores and stratifying patients
8. Generating Kaplan–Meier plots
9. Performing 5-fold cross-validation for model performance

## Input Data Format

Place the dataset in:

```
input/dataset.xlsx
```

### Required columns:

- `patient` – unique patient ID
- `rPFS_days` – time to progression or censoring
- `rPFS` – progression event indicator (1 = event, 0 = censored)
- Biomarker features (numeric):

```
FAT1, KLK2, STEAP2, TMPRSS2, AGR2, FOLH1, HOXB13, KLK3, CTCM, MYC,
DLL3, SYP, ARV7, TACSTD2, EZH2, CHGA, NEUROD1, PDX1, E2F1,
SUV_mean, tumor_volume, Age, PSA0
```

## Running the Script

First, install the required dependencies:

```
pip install -r requirements.txt
```

Then run the analysis:

```
python survival_analysis.py
```

All results will be saved in the `results/` directory.

## Output Files

| File                    | Description                                      |
|-------------------------|--------------------------------------------------|
| `cox_summary.csv`       | Summary of full model coefficients and p-values |
| `volcano_plot.png`      | Volcano plot of all features                    |
| `selected_features.txt` | Selected features used in the final model       |
| `risk_scores.xlsx`      | Patient risk scores and group assignments       |
| `km_plot.png`           | Kaplan–Meier plot with HR and p-value           |

## Model Example Summary

- Final features: SUV_mean, ARV7, E2F1, STEAP2
- Hazard Ratio (high vs. low risk): 10.43
- Log-rank p-value: 0.0003
- Cross-validated C-index: 0.87 ± 0.07


## Contact

For questions or data access, contact @ichamseddine or open a GitHub issue.
