import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from lifelines.utils import k_fold_cross_validation
from statsmodels.stats.outliers_influence import variance_inflation_factor

# ==== CONFIGURATION ====
DATA_FILE = 'input/dataset.xlsx'
OUTPUT_DIR = 'results'
TIME_COLUMN = 'rPFS_days'
EVENT_COLUMN = 'rPFS'
ALL_FEATURES = [
    'FAT1', 'KLK2', 'STEAP2', 'TMPRSS2', 'AGR2', 'FOLH1', 'HOXB13', 'KLK3', 'CTCM', 'MYC',
    'DLL3', 'SYP', 'ARV7', 'TACSTD2', 'EZH2', 'CHGA', 'NEUROD1', 'PDX1', 'E2F1',
    'SUV_mean', 'tumor_volume', 'Age', 'PSA0'
]
P_VALUE_THRESHOLD = 0.2
TOP_N_FEATURES = 7
PENALIZER = 0.5
N_FOLDS = 5

# ==== UTILS ====

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def calculate_vif(df):
    vif_data = pd.DataFrame()
    vif_data["feature"] = df.columns
    vif_data["VIF"] = [variance_inflation_factor(df.values, i) for i in range(df.shape[1])]
    return vif_data

def load_and_filter_data(file_path, time_col, event_col, feature_cols):
    df = pd.read_excel(file_path)
    df = df[feature_cols + [time_col, event_col]].dropna()

    dropped = []
    vif_df = calculate_vif(df[feature_cols])
    while vif_df['VIF'].max() > 10:
        feature_to_drop = vif_df.sort_values('VIF', ascending=False).iloc[0]['feature']
        dropped.append(feature_to_drop)
        feature_cols.remove(feature_to_drop)
        vif_df = calculate_vif(df[feature_cols])

    print("Dropped due to VIF:", dropped)
    print("Remaining features:", feature_cols)
    df = df[feature_cols + [time_col, event_col]]
    return df, feature_cols

def fit_cox_model(df, time_col, event_col, penalizer=0.0):
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(df, duration_col=time_col, event_col=event_col)
    return cph

def extract_feature_stats(cph):
    summary = cph.summary.copy()
    summary['Feature'] = summary.index
    summary['Coefficient'] = summary['coef']
    summary['p_value'] = summary['p']
    return summary[['Feature', 'Coefficient', 'p_value']]

def plot_volcano(feature_stats, threshold, save_path):
    feature_stats['-log10(p_value)'] = -np.log10(feature_stats['p_value'])
    plt.figure(figsize=(5, 5))
    plt.scatter(feature_stats['Coefficient'], feature_stats['-log10(p_value)'], c='blue', alpha=0.7)
    sig = feature_stats['p_value'] < threshold
    plt.scatter(feature_stats['Coefficient'][sig], feature_stats['-log10(p_value)'][sig], c='red', alpha=0.7)
    for _, row in feature_stats.iterrows():
        plt.text(row['Coefficient'], row['-log10(p_value)'] + 0.05, row['Feature'], fontsize=8, ha='center')
    plt.axhline(-np.log10(threshold), color='gray', linestyle='--')
    plt.axvline(0, color='black', linestyle='--')
    plt.xlabel('Coefficient')
    plt.ylabel('-log10(p-value)')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=600)
    plt.show()

def select_features(feature_stats, pval_thresh, top_n):
    filtered = feature_stats[feature_stats['p_value'] < pval_thresh]
    sorted_features = filtered.reindex(filtered['Coefficient'].abs().sort_values(ascending=False).index)
    return sorted_features['Feature'].head(top_n).tolist()

def compute_risk_scores(df, features, time_col, event_col, penalizer=0.5):
    df_model = df[features + [time_col, event_col]].copy()
    cph = fit_cox_model(df_model, time_col, event_col, penalizer=penalizer)
    df_model['risk_score'] = cph.predict_partial_hazard(df_model)
    return df_model, cph

def stratify_patients(df, score_col='risk_score'):
    median_score = df[score_col].median()
    df['risk_group'] = df[score_col].apply(lambda x: 'high' if x > median_score else 'low')
    return df

def plot_km(df, time_col, event_col, output_path):
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(6, 5))
    colors = {'high': 'red', 'low': 'blue'}
    for group in ['high', 'low']:
        ix = df['risk_group'] == group
        kmf.fit(df.loc[ix, time_col], event_observed=df.loc[ix, event_col], label=group)
        kmf.plot_survival_function(ci_show=True, color=colors[group])
    add_at_risk_counts(kmf, ax=plt.gca())

    high = df[df['risk_group'] == 'high']
    low = df[df['risk_group'] == 'low']
    results = logrank_test(high[time_col], low[time_col], event_observed_A=high[event_col], event_observed_B=low[event_col])
    pval = results.p_value

    strat_df = df[['risk_group', time_col, event_col]].copy()
    strat_df['risk_group'] = strat_df['risk_group'].map({'low': 0, 'high': 1})
    hr_model = fit_cox_model(strat_df, time_col, event_col)
    hr = hr_model.hazard_ratios_['risk_group']

    plt.text(plt.xlim()[1] * 0.6, 0.2, f"p = {pval:.4f}\nHR = {hr:.2f}", fontsize=11, bbox=dict(facecolor='white', alpha=0.5))
    plt.title('Survival Curves by Risk Group')
    plt.xlabel('Time (days)')
    plt.ylabel('Survival Probability')
    plt.ylim(0, 1)
    plt.grid(True)
    plt.savefig(output_path, dpi=600)
    plt.show()

def cross_validate_model(df, features, time_col, event_col, penalizer, n_folds):
    df_cv = df[features + [time_col, event_col]].copy()
    cph = CoxPHFitter(penalizer=penalizer)
    scores = k_fold_cross_validation(cph, df_cv, duration_col=time_col, event_col=event_col, k=n_folds)
    return np.mean(scores), np.std(scores)

# ==== MAIN ====
def main():
    ensure_dir(OUTPUT_DIR)

    df, filtered_features = load_and_filter_data(DATA_FILE, TIME_COLUMN, EVENT_COLUMN, ALL_FEATURES.copy())

    cph_all = fit_cox_model(df, TIME_COLUMN, EVENT_COLUMN)
    stats = extract_feature_stats(cph_all)
    stats.to_csv(f'{OUTPUT_DIR}/cox_summary.csv', index=False)

    plot_volcano(stats, P_VALUE_THRESHOLD, save_path=f'{OUTPUT_DIR}/volcano_plot.png')

    selected = select_features(stats, P_VALUE_THRESHOLD, TOP_N_FEATURES)
    print("Selected features:", selected)
    with open(f'{OUTPUT_DIR}/selected_features.txt', 'w') as f:
        f.write('\n'.join(selected))

    df_model, cph = compute_risk_scores(df, selected, TIME_COLUMN, EVENT_COLUMN, penalizer=PENALIZER)
    df_model = stratify_patients(df_model)
    df_model[['risk_score', 'risk_group', TIME_COLUMN, EVENT_COLUMN]].to_excel(f'{OUTPUT_DIR}/risk_scores.xlsx', index=False)

    plot_km(df_model, TIME_COLUMN, EVENT_COLUMN, output_path=f'{OUTPUT_DIR}/km_plot.png')

    mean_c, std_c = cross_validate_model(df_model, selected, TIME_COLUMN, EVENT_COLUMN, PENALIZER, N_FOLDS)
    print(f"C-index: {mean_c:.4f} ± {std_c:.4f}")

if __name__ == '__main__':
    main()
