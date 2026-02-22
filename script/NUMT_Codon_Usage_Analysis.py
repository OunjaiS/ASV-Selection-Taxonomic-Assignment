"""
NUMT Validation: Codon Usage Bias Analysis
============================================
Addresses Reviewer 1, Major Comment #1.

Compares codon usage patterns between Authenticated ASVs and putative NUMTs
to provide molecular evidence supporting NUMT classification.

4 analyses:
  1. Nucleotide composition (GC/AT content, position-specific GC)
  2. Effective Number of Codons (ENC)
  3. Per-ASV RSCU + PCA
  4. Amino acid composition + case study selection
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CONFIG
# ============================================================
DATA_PATH = '/Users/sarawut/Desktop/ASV_selection/data/ASV_Complete_Analysis_MLClassified.csv'
HC_NUMTS_PATH = '/Users/sarawut/Desktop/ASV_selection/data/ASV_Complete_Analysis_NUMPs_Abundance_Classified_HighConfidence_NUMTs.csv'
OUTPUT_DIR = '/Users/sarawut/Desktop/Manuscript_ASV_selection/data_analysis/numt_validation/'
FIG_DIR = '/Users/sarawut/Desktop/Manuscript_ASV_selection/figures/'

import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Invertebrate Mitochondrial Codon Table (NCBI Table 5)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'S', 'AGG': 'S',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Group codons by amino acid
AA_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa == '*':
        continue
    AA_CODONS.setdefault(aa, []).append(codon)

ALL_CODONS = sorted([c for c, aa in CODON_TABLE.items() if aa != '*'])

# ============================================================
# FUNCTIONS
# ============================================================

def calculate_rscu(sequence):
    """Calculate RSCU for a single DNA sequence."""
    sequence = sequence.upper()
    codons = [sequence[i:i+3] for i in range(0, (len(sequence) // 3) * 3, 3)]
    codon_counts = Counter(codons)

    rscu = {}
    for aa, syn_codons in AA_CODONS.items():
        n_syn = len(syn_codons)
        total = sum(codon_counts.get(c, 0) for c in syn_codons)
        for c in syn_codons:
            if total > 0:
                rscu[c] = (codon_counts.get(c, 0) * n_syn) / total
            else:
                rscu[c] = np.nan
    return rscu


def effect_size_rank_biserial(u_stat, n1, n2):
    """Calculate rank-biserial correlation as effect size for Mann-Whitney U."""
    return 1 - (2 * u_stat) / (n1 * n2)


# ============================================================
# 1. LOAD DATA
# ============================================================
print("=" * 70)
print("NUMT VALIDATION: CODON USAGE BIAS ANALYSIS")
print("=" * 70)

print("\n[1/6] Loading data...")
df = pd.read_csv(DATA_PATH, low_memory=False)
print(f"  Total records: {len(df):,}")

# Get unique ASVs (deduplicate by asv_id, keep first occurrence)
df_unique = df.drop_duplicates(subset='asv_id', keep='first').copy()
print(f"  Unique ASVs: {len(df_unique):,}")

# Filter groups
auth = df_unique[df_unique['ml_final_classification'] == 'Authenticated'].copy()
numts = df_unique[df_unique['ml_final_classification'] == 'Nuclear_Pseudogenes'].copy()
tech = df_unique[df_unique['ml_final_classification'] == 'Technical_Artifacts'].copy()

print(f"  Authenticated (unique): {len(auth):,}")
print(f"  Nuclear_Pseudogenes (unique): {len(numts):,}")
print(f"  Technical_Artifacts (unique): {len(tech):,}")

# Filter sequences that are valid
auth = auth[auth['corrected_sequence'].notna() & (auth['corrected_sequence'].str.len() >= 300)].copy()
numts = numts[numts['corrected_sequence'].notna() & (numts['corrected_sequence'].str.len() >= 300)].copy()
tech = tech[tech['corrected_sequence'].notna() & (tech['corrected_sequence'].str.len() >= 300)].copy()

print(f"  After sequence filter - Auth: {len(auth):,}, NUMTs: {len(numts):,}, Tech: {len(tech):,}")

# ============================================================
# 2. NUCLEOTIDE COMPOSITION ANALYSIS
# ============================================================
print("\n[2/6] Nucleotide Composition Analysis...")

metrics = {
    'GC_content': 'GC Content (%)',
    'AT_content': 'AT Content (%)',
    'gc_pos3': 'GC at 3rd Codon Position (%)',
    'gc_pos1': 'GC at 1st Codon Position (%)',
    'gc_pos2': 'GC at 2nd Codon Position (%)',
}

comp_results = []
for col, label in metrics.items():
    a_vals = auth[col].dropna()
    n_vals = numts[col].dropna()

    if len(a_vals) == 0 or len(n_vals) == 0:
        print(f"  SKIP {col}: no data")
        continue

    u_stat, p_val = stats.mannwhitneyu(a_vals, n_vals, alternative='two-sided')
    r_bs = effect_size_rank_biserial(u_stat, len(a_vals), len(n_vals))

    comp_results.append({
        'Metric': label,
        'Column': col,
        'Auth_mean': a_vals.mean(),
        'Auth_sd': a_vals.std(),
        'Auth_median': a_vals.median(),
        'NUMT_mean': n_vals.mean(),
        'NUMT_sd': n_vals.std(),
        'NUMT_median': n_vals.median(),
        'U_statistic': u_stat,
        'p_value': p_val,
        'rank_biserial_r': r_bs,
        'n_Auth': len(a_vals),
        'n_NUMT': len(n_vals),
    })

    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    print(f"  {label}:")
    print(f"    Auth: {a_vals.mean():.2f} ± {a_vals.std():.2f} (median: {a_vals.median():.2f})")
    print(f"    NUMT: {n_vals.mean():.2f} ± {n_vals.std():.2f} (median: {n_vals.median():.2f})")
    print(f"    Mann-Whitney U = {u_stat:.0f}, p = {p_val:.2e}, r = {r_bs:.3f} {sig}")

comp_df = pd.DataFrame(comp_results)
comp_df.to_csv(os.path.join(OUTPUT_DIR, 'nucleotide_composition_comparison.csv'), index=False)

# --- Figure: Violin plots ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
plot_metrics = ['GC_content', 'gc_pos3', 'AT_content']
plot_labels = ['GC Content (%)', 'GC at 3rd Codon Position (%)', 'AT Content (%)']

for ax, col, label in zip(axes, plot_metrics, plot_labels):
    plot_data = pd.DataFrame({
        'Value': pd.concat([auth[col].dropna(), numts[col].dropna()], ignore_index=True),
        'Group': ['Authenticated'] * len(auth[col].dropna()) + ['NUMT'] * len(numts[col].dropna())
    })
    sns.violinplot(data=plot_data, x='Group', y='Value', ax=ax,
                   palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
                   inner='box', cut=0)
    ax.set_title(label, fontsize=11, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel(label)

    # Add p-value annotation
    a_vals = auth[col].dropna()
    n_vals = numts[col].dropna()
    _, p = stats.mannwhitneyu(a_vals, n_vals, alternative='two-sided')
    sig_text = f"p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    ymax = max(a_vals.max(), n_vals.max())
    ax.text(0.5, ymax * 1.02, sig_text, ha='center', fontsize=9, fontstyle='italic')

plt.suptitle('Nucleotide Composition: Authenticated ASVs vs Putative NUMTs',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_NucleotideComposition.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_NucleotideComposition.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: NUMT_Validation_NucleotideComposition.png/pdf")

# ============================================================
# 3. EFFECTIVE NUMBER OF CODONS (ENC) ANALYSIS
# ============================================================
print("\n[3/6] Effective Number of Codons (ENC) Analysis...")

enc_auth = auth['effective_codon_number'].dropna()
enc_numt = numts['effective_codon_number'].dropna()

u_enc, p_enc = stats.mannwhitneyu(enc_auth, enc_numt, alternative='two-sided')
r_enc = effect_size_rank_biserial(u_enc, len(enc_auth), len(enc_numt))

print(f"  Auth ENC: {enc_auth.mean():.2f} ± {enc_auth.std():.2f} (median: {enc_auth.median():.2f})")
print(f"  NUMT ENC: {enc_numt.mean():.2f} ± {enc_numt.std():.2f} (median: {enc_numt.median():.2f})")
print(f"  Mann-Whitney U = {u_enc:.0f}, p = {p_enc:.2e}, r = {r_enc:.3f}")

# Also compare codon_diversity and codon_bias_score
for col_name in ['codon_diversity', 'codon_bias_score']:
    a_v = auth[col_name].dropna()
    n_v = numts[col_name].dropna()
    if len(a_v) > 0 and len(n_v) > 0:
        u, p = stats.mannwhitneyu(a_v, n_v, alternative='two-sided')
        r = effect_size_rank_biserial(u, len(a_v), len(n_v))
        print(f"  {col_name}: Auth={a_v.mean():.2f}±{a_v.std():.2f}, NUMT={n_v.mean():.2f}±{n_v.std():.2f}, p={p:.2e}, r={r:.3f}")

# --- Figure: ENC boxplot ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: ENC
enc_data = pd.DataFrame({
    'Value': pd.concat([enc_auth, enc_numt], ignore_index=True),
    'Group': ['Authenticated'] * len(enc_auth) + ['NUMT'] * len(enc_numt)
})
sns.boxplot(data=enc_data, x='Group', y='Value', ax=axes[0],
            palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
            showfliers=False)
sns.stripplot(data=enc_data, x='Group', y='Value', ax=axes[0],
              palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
              alpha=0.05, size=1)
axes[0].set_title('Effective Number of Codons (ENC)', fontsize=11, fontweight='bold')
axes[0].set_ylabel('ENC')
axes[0].set_xlabel('')

# Panel B: Codon Diversity
cd_auth = auth['codon_diversity'].dropna()
cd_numt = numts['codon_diversity'].dropna()
cd_data = pd.DataFrame({
    'Value': pd.concat([cd_auth, cd_numt], ignore_index=True),
    'Group': ['Authenticated'] * len(cd_auth) + ['NUMT'] * len(cd_numt)
})
sns.boxplot(data=cd_data, x='Group', y='Value', ax=axes[1],
            palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
            showfliers=False)
sns.stripplot(data=cd_data, x='Group', y='Value', ax=axes[1],
              palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
              alpha=0.05, size=1)
axes[1].set_title('Codon Diversity', fontsize=11, fontweight='bold')
axes[1].set_ylabel('Codon Diversity Score')
axes[1].set_xlabel('')

# Panel C: Codon Bias Score
cb_auth = auth['codon_bias_score'].dropna()
cb_numt = numts['codon_bias_score'].dropna()
cb_data = pd.DataFrame({
    'Value': pd.concat([cb_auth, cb_numt], ignore_index=True),
    'Group': ['Authenticated'] * len(cb_auth) + ['NUMT'] * len(cb_numt)
})
sns.boxplot(data=cb_data, x='Group', y='Value', ax=axes[2],
            palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
            showfliers=False)
sns.stripplot(data=cb_data, x='Group', y='Value', ax=axes[2],
              palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
              alpha=0.05, size=1)
axes[2].set_title('Codon Bias Score', fontsize=11, fontweight='bold')
axes[2].set_ylabel('Codon Bias Score')
axes[2].set_xlabel('')

plt.suptitle('Codon Usage Bias: Authenticated ASVs vs Putative NUMTs',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_CodonBias.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_CodonBias.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: NUMT_Validation_CodonBias.png/pdf")

# ============================================================
# 4. PER-ASV RSCU CALCULATION + PCA
# ============================================================
print("\n[4/6] Calculating per-ASV RSCU (this may take a moment)...")

# Use a sample if datasets are very large
MAX_PER_GROUP = 5000
auth_sample = auth.sample(n=min(MAX_PER_GROUP, len(auth)), random_state=42)
numt_sample = numts.sample(n=min(MAX_PER_GROUP, len(numts)), random_state=42)
tech_sample = tech.sample(n=min(2000, len(tech)), random_state=42)

combined = pd.concat([auth_sample, numt_sample, tech_sample], ignore_index=True)
print(f"  Computing RSCU for {len(combined):,} ASVs...")

rscu_records = []
for idx, row in combined.iterrows():
    seq = row['corrected_sequence']
    if not isinstance(seq, str) or len(seq) < 300:
        continue
    rscu = calculate_rscu(seq)
    rscu['group'] = row['ml_final_classification']
    rscu['asv_id'] = row['asv_id']
    rscu_records.append(rscu)

rscu_df = pd.DataFrame(rscu_records)
print(f"  RSCU computed for {len(rscu_df):,} ASVs")

# Save RSCU matrix
rscu_df.to_csv(os.path.join(OUTPUT_DIR, 'per_asv_rscu_matrix.csv'), index=False)

# --- PCA ---
print("  Running PCA on RSCU matrix...")
codon_cols = [c for c in ALL_CODONS if c in rscu_df.columns]
X = rscu_df[codon_cols].fillna(rscu_df[codon_cols].mean()).values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=5)
pcs = pca.fit_transform(X_scaled)

print(f"  PC1: {pca.explained_variance_ratio_[0]:.1%}")
print(f"  PC2: {pca.explained_variance_ratio_[1]:.1%}")
print(f"  PC3: {pca.explained_variance_ratio_[2]:.1%}")
print(f"  Total (PC1-2): {sum(pca.explained_variance_ratio_[:2]):.1%}")

# --- Figure: PCA plot ---
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Panel A: Auth vs NUMTs only
group_colors = {
    'Authenticated': '#2196F3',
    'Nuclear_Pseudogenes': '#F44336',
    'Technical_Artifacts': '#9E9E9E',
}
group_labels = {
    'Authenticated': 'Authenticated',
    'Nuclear_Pseudogenes': 'NUMT',
    'Technical_Artifacts': 'Technical Artifact',
}

for g in ['Technical_Artifacts', 'Nuclear_Pseudogenes', 'Authenticated']:
    mask = rscu_df['group'] == g
    if mask.sum() == 0:
        continue
    axes[0].scatter(pcs[mask, 0], pcs[mask, 1],
                    alpha=0.25, s=8, c=group_colors[g],
                    label=group_labels[g], edgecolors='none')

axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=11)
axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=11)
axes[0].set_title('PCA of Codon Usage (RSCU)', fontsize=12, fontweight='bold')
axes[0].legend(fontsize=9, markerscale=3)

# Panel B: Density contours for Auth vs NUMTs
for g, color in [('Authenticated', '#2196F3'), ('Nuclear_Pseudogenes', '#F44336')]:
    mask = rscu_df['group'] == g
    if mask.sum() < 10:
        continue
    x, y = pcs[mask, 0], pcs[mask, 1]
    try:
        sns.kdeplot(x=x, y=y, ax=axes[1], color=color,
                    levels=5, linewidths=1.5, label=group_labels[g])
    except Exception:
        axes[1].scatter(x, y, alpha=0.3, s=8, c=color, label=group_labels[g])

axes[1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=11)
axes[1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=11)
axes[1].set_title('Density Contours: Authenticated vs NUMT', fontsize=12, fontweight='bold')
axes[1].legend(fontsize=9)

plt.suptitle('RSCU-based PCA: Codon Usage Patterns across ASV Categories',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_RSCU_PCA.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_RSCU_PCA.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: NUMT_Validation_RSCU_PCA.png/pdf")

# --- RSCU Difference Bar Chart ---
print("  Generating RSCU difference chart...")
auth_mask = rscu_df['group'] == 'Authenticated'
numt_mask = rscu_df['group'] == 'Nuclear_Pseudogenes'

auth_mean_rscu = rscu_df.loc[auth_mask, codon_cols].mean()
numt_mean_rscu = rscu_df.loc[numt_mask, codon_cols].mean()
delta_rscu = numt_mean_rscu - auth_mean_rscu
delta_sorted = delta_rscu.sort_values()

# Statistical test per codon
codon_pvals = {}
for c in codon_cols:
    a_v = rscu_df.loc[auth_mask, c].dropna()
    n_v = rscu_df.loc[numt_mask, c].dropna()
    if len(a_v) > 10 and len(n_v) > 10:
        _, p = stats.mannwhitneyu(a_v, n_v, alternative='two-sided')
        codon_pvals[c] = p
    else:
        codon_pvals[c] = 1.0

fig, ax = plt.subplots(figsize=(16, 6))
bar_colors = ['#F44336' if v > 0 else '#2196F3' for v in delta_sorted.values]

# Mark significant codons
edge_colors = []
for c in delta_sorted.index:
    if codon_pvals.get(c, 1.0) < 0.001:
        edge_colors.append('black')
    else:
        edge_colors.append('none')

bars = ax.bar(range(len(delta_sorted)), delta_sorted.values, color=bar_colors,
              edgecolor=edge_colors, linewidth=0.8)
ax.set_xticks(range(len(delta_sorted)))
ax.set_xticklabels([f"{c}\n({CODON_TABLE.get(c, '?')})" for c in delta_sorted.index],
                   rotation=90, fontsize=6.5)
ax.set_ylabel('ΔRSCU (NUMT − Authenticated)', fontsize=11)
ax.set_title('Difference in Relative Synonymous Codon Usage between NUMTs and Authenticated ASVs',
             fontsize=12, fontweight='bold')
ax.axhline(y=0, color='black', linewidth=0.5)
ax.text(0.02, 0.95, 'Black border = p < 0.001', transform=ax.transAxes,
        fontsize=8, fontstyle='italic', va='top')
ax.text(0.02, 0.89, 'Red = NUMT uses more | Blue = Authenticated uses more',
        transform=ax.transAxes, fontsize=8, fontstyle='italic', va='top')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_RSCU_Difference.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'NUMT_Validation_RSCU_Difference.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: NUMT_Validation_RSCU_Difference.png/pdf")

# Save RSCU comparison table
rscu_comp = pd.DataFrame({
    'codon': codon_cols,
    'amino_acid': [CODON_TABLE.get(c, '?') for c in codon_cols],
    'Auth_mean_RSCU': [auth_mean_rscu[c] for c in codon_cols],
    'NUMT_mean_RSCU': [numt_mean_rscu[c] for c in codon_cols],
    'delta_RSCU': [delta_rscu[c] for c in codon_cols],
    'p_value': [codon_pvals.get(c, np.nan) for c in codon_cols],
})
rscu_comp = rscu_comp.sort_values('delta_RSCU')

# Apply Benjamini-Hochberg FDR correction
valid_pvals = rscu_comp['p_value'].dropna()
reject, pvals_corrected, _, _ = multipletests(valid_pvals, alpha=0.05, method='fdr_bh')
rscu_comp.loc[valid_pvals.index, 'p_adjusted_BH'] = pvals_corrected
rscu_comp['significant'] = rscu_comp['p_adjusted_BH'] < 0.05
rscu_comp.to_csv(os.path.join(OUTPUT_DIR, 'rscu_comparison_per_codon.csv'), index=False)

n_sig = rscu_comp['significant'].sum()
print(f"  Significant codons (BH-adjusted p < 0.05): {n_sig}/{len(codon_cols)}")

# ============================================================
# 5. AMINO ACID COMPOSITION
# ============================================================
print("\n[5/6] Amino Acid Composition Analysis...")

aa_metrics = {
    'hydrophobic_percent': 'Hydrophobic AA (%)',
    'leucine_percent': 'Leucine (%)',
    'internal_stops': 'Internal Stop Codons',
}

aa_results = []
for col, label in aa_metrics.items():
    a_v = auth[col].dropna()
    n_v = numts[col].dropna()
    if len(a_v) == 0 or len(n_v) == 0:
        continue
    u, p = stats.mannwhitneyu(a_v, n_v, alternative='two-sided')
    r = effect_size_rank_biserial(u, len(a_v), len(n_v))
    aa_results.append({
        'Metric': label, 'Column': col,
        'Auth_mean': a_v.mean(), 'Auth_sd': a_v.std(),
        'NUMT_mean': n_v.mean(), 'NUMT_sd': n_v.std(),
        'U_stat': u, 'p_value': p, 'rank_biserial_r': r,
    })
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {label}:")
    print(f"    Auth: {a_v.mean():.3f} ± {a_v.std():.3f}")
    print(f"    NUMT: {n_v.mean():.3f} ± {n_v.std():.3f}")
    print(f"    p = {p:.2e}, r = {r:.3f} {sig}")

aa_df = pd.DataFrame(aa_results)
aa_df.to_csv(os.path.join(OUTPUT_DIR, 'amino_acid_composition_comparison.csv'), index=False)

# ============================================================
# 6. CASE STUDY SELECTION
# ============================================================
print("\n[6/6] Selecting NUMT Case Studies...")

hc = pd.read_csv(HC_NUMTS_PATH, low_memory=False)
print(f"  High Confidence NUMTs: {len(hc):,} records")

# Find concordant NUMTs (both abundance-based and ML classify as NUMT)
concordant = hc[hc['ml_final_classification'] == 'Nuclear_Pseudogenes'].copy()
print(f"  Concordant (ML=Nuclear_Pseudogenes): {len(concordant):,}")

# Filter biologically plausible NUMTs for case studies:
# - Low abundance (< 10% reads) to match NUMT expectation
# - NUMT GC should differ from Main ASV GC
# - Main ASV should have typical insect mtDNA GC range (25-40%)
def get_main_gc(row):
    main_id = row.get('Main_ASV', '')
    main_row = df_unique[df_unique['asv_id'] == main_id]
    if len(main_row) > 0:
        return main_row.iloc[0].get('GC_content', np.nan)
    return np.nan

concordant['Main_GC_check'] = concordant.apply(get_main_gc, axis=1)
plausible = concordant[
    (concordant['percentage_reads'] < 10) &  # Low abundance
    (concordant['Main_GC_check'].between(25, 42))  # Main ASV has typical insect mtDNA GC
].copy()

print(f"  Biologically plausible cases (low abundance + normal Main ASV GC): {len(plausible):,}")

# Score for case study selection
plausible['case_score'] = (
    plausible['Phylogenetic_distance'].fillna(0) * 10 +  # Higher distance = better
    (100 - plausible['percentage_reads'].fillna(100)) * 0.5 +  # Lower abundance = better
    plausible['GC_content'].fillna(0) * 0.3  # Higher GC = more NUMT-like
)

# Select top candidates from different families
concordant_sorted = plausible.sort_values('case_score', ascending=False)
selected_families = set()
case_studies = []

for _, row in concordant_sorted.iterrows():
    fam = row.get('family', 'Unknown')
    if fam in selected_families:
        continue
    selected_families.add(fam)
    case_studies.append(row)
    if len(case_studies) >= 3:
        break

print(f"\n  Selected {len(case_studies)} case studies:")
case_study_data = []
for i, cs in enumerate(case_studies, 1):
    print(f"\n  --- Case Study {i} ---")
    print(f"  Specimen: {cs.get('project_sample_id', 'N/A')}")
    print(f"  Family: {cs.get('family', 'N/A')}")
    print(f"  Country: {cs.get('country', 'N/A')}")
    print(f"  ASV ID: {cs.get('asv_id', 'N/A')}")
    print(f"  Main ASV: {cs.get('Main_ASV', 'N/A')}")
    print(f"  Reads: {cs.get('reads', 'N/A')}")
    print(f"  % Reads: {cs.get('percentage_reads', 'N/A'):.2f}%")
    print(f"  Phylogenetic Distance: {cs.get('Phylogenetic_distance', 'N/A')}")
    print(f"  GC Content: {cs.get('GC_content', 'N/A'):.2f}%")
    print(f"  ENC: {cs.get('effective_codon_number', 'N/A')}")
    print(f"  Internal Stops: {cs.get('internal_stops', 'N/A')}")
    print(f"  ML Classification: {cs.get('ml_final_classification', 'N/A')}")
    print(f"  ML Confidence: {cs.get('ml_confidence', 'N/A')}")
    print(f"  NUMPs Label: {cs.get('NUMPs_Label', 'N/A')}")

    # Get Main ASV data for comparison
    main_asv_id = cs.get('Main_ASV', '')
    main_row = df_unique[df_unique['asv_id'] == main_asv_id]
    if len(main_row) > 0:
        mr = main_row.iloc[0]
        print(f"  [Main ASV] GC: {mr.get('GC_content', 'N/A')}, ENC: {mr.get('effective_codon_number', 'N/A')}")
        case_study_data.append({
            'Case': i,
            'Specimen': cs.get('project_sample_id', ''),
            'Family': cs.get('family', ''),
            'Country': cs.get('country', ''),
            'NUMT_ASV': cs.get('asv_id', ''),
            'Main_ASV': main_asv_id,
            'NUMT_reads': cs.get('reads', ''),
            'NUMT_pct_reads': cs.get('percentage_reads', ''),
            'Phylo_distance': cs.get('Phylogenetic_distance', ''),
            'NUMT_GC': cs.get('GC_content', ''),
            'Main_GC': mr.get('GC_content', ''),
            'NUMT_ENC': cs.get('effective_codon_number', ''),
            'Main_ENC': mr.get('effective_codon_number', ''),
            'NUMT_internal_stops': cs.get('internal_stops', ''),
            'Main_internal_stops': mr.get('internal_stops', ''),
            'NUMT_hydrophobic': cs.get('hydrophobic_percent', ''),
            'Main_hydrophobic': mr.get('hydrophobic_percent', ''),
            'ML_confidence': cs.get('ml_confidence', ''),
            'NUMPs_Label': cs.get('NUMPs_Label', ''),
        })

if case_study_data:
    cs_df = pd.DataFrame(case_study_data)
    cs_df.to_csv(os.path.join(OUTPUT_DIR, 'numt_case_studies.csv'), index=False)
    print("\n  Saved: numt_case_studies.csv")

# ============================================================
# COMBINED SUPPLEMENTARY FIGURE
# ============================================================
print("\n[BONUS] Generating combined supplementary figure...")

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Panel A: GC Content violin
plot_data_gc = pd.DataFrame({
    'Value': pd.concat([auth['GC_content'].dropna(), numts['GC_content'].dropna()], ignore_index=True),
    'Group': ['Authenticated'] * len(auth['GC_content'].dropna()) + ['NUMT'] * len(numts['GC_content'].dropna())
})
sns.violinplot(data=plot_data_gc, x='Group', y='Value', ax=axes[0, 0],
               palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
               inner='box', cut=0)
axes[0, 0].set_title('(A) GC Content', fontsize=12, fontweight='bold')
axes[0, 0].set_ylabel('GC Content (%)')
axes[0, 0].set_xlabel('')

# Panel B: ENC violin
plot_data_enc = pd.DataFrame({
    'Value': pd.concat([enc_auth, enc_numt], ignore_index=True),
    'Group': ['Authenticated'] * len(enc_auth) + ['NUMT'] * len(enc_numt)
})
sns.violinplot(data=plot_data_enc, x='Group', y='Value', ax=axes[0, 1],
               palette={'Authenticated': '#2196F3', 'NUMT': '#F44336'},
               inner='box', cut=0)
axes[0, 1].set_title('(B) Effective Number of Codons', fontsize=12, fontweight='bold')
axes[0, 1].set_ylabel('ENC')
axes[0, 1].set_xlabel('')

# Panel C: PCA
for g in ['Technical_Artifacts', 'Nuclear_Pseudogenes', 'Authenticated']:
    mask = rscu_df['group'] == g
    if mask.sum() == 0:
        continue
    axes[1, 0].scatter(pcs[mask, 0], pcs[mask, 1],
                       alpha=0.2, s=6, c=group_colors[g],
                       label=group_labels[g], edgecolors='none')
axes[1, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
axes[1, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
axes[1, 0].set_title('(C) PCA of RSCU Patterns', fontsize=12, fontweight='bold')
axes[1, 0].legend(fontsize=8, markerscale=3)

# Panel D: Top RSCU differences
top_n = 20
top_diff = pd.concat([delta_sorted.head(top_n // 2), delta_sorted.tail(top_n // 2)])
bar_c = ['#2196F3' if v < 0 else '#F44336' for v in top_diff.values]
axes[1, 1].barh(range(len(top_diff)), top_diff.values, color=bar_c)
axes[1, 1].set_yticks(range(len(top_diff)))
axes[1, 1].set_yticklabels([f"{c} ({CODON_TABLE.get(c, '?')})" for c in top_diff.index], fontsize=8)
axes[1, 1].set_xlabel('ΔRSCU (NUMT − Authenticated)')
axes[1, 1].set_title('(D) Top RSCU Differences', fontsize=12, fontweight='bold')
axes[1, 1].axvline(x=0, color='black', linewidth=0.5)

plt.suptitle('Supplementary Figure: NUMT Validation through Codon Usage Bias Analysis',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'Supplementary_NUMT_Validation.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'Supplementary_NUMT_Validation.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: Supplementary_NUMT_Validation.png/pdf")

# ============================================================
# SUMMARY REPORT
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE — SUMMARY")
print("=" * 70)

print(f"""
Sample sizes (unique ASVs with valid sequences):
  Authenticated: {len(auth):,}
  Nuclear_Pseudogenes (NUMTs): {len(numts):,}
  Technical_Artifacts: {len(tech):,}

Key Results:
  1. Nucleotide Composition:""")

for r in comp_results:
    direction = "NUMT > Auth" if r['NUMT_mean'] > r['Auth_mean'] else "Auth > NUMT"
    sig = "***" if r['p_value'] < 0.001 else "ns"
    print(f"     {r['Metric']}: Auth={r['Auth_mean']:.2f}%, NUMT={r['NUMT_mean']:.2f}% ({direction}) {sig}")

print(f"""
  2. Effective Number of Codons:
     Auth ENC: {enc_auth.mean():.2f} ± {enc_auth.std():.2f}
     NUMT ENC: {enc_numt.mean():.2f} ± {enc_numt.std():.2f}
     Direction: {"NUMT > Auth (less bias)" if enc_numt.mean() > enc_auth.mean() else "Auth > NUMT"}
     p = {p_enc:.2e}

  3. RSCU PCA:
     PC1 explains {pca.explained_variance_ratio_[0]:.1%} of variance
     PC2 explains {pca.explained_variance_ratio_[1]:.1%} of variance
     Significant codons: {n_sig}/{len(codon_cols)}

  4. Case Studies: {len(case_studies)} selected

Output Files:
  Figures:
    - {FIG_DIR}NUMT_Validation_NucleotideComposition.png/pdf
    - {FIG_DIR}NUMT_Validation_CodonBias.png/pdf
    - {FIG_DIR}NUMT_Validation_RSCU_PCA.png/pdf
    - {FIG_DIR}NUMT_Validation_RSCU_Difference.png/pdf
    - {FIG_DIR}Supplementary_NUMT_Validation.png/pdf

  Data Tables:
    - {OUTPUT_DIR}nucleotide_composition_comparison.csv
    - {OUTPUT_DIR}amino_acid_composition_comparison.csv
    - {OUTPUT_DIR}rscu_comparison_per_codon.csv
    - {OUTPUT_DIR}per_asv_rscu_matrix.csv
    - {OUTPUT_DIR}numt_case_studies.csv
""")
print("Done!")
