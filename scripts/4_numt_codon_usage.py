"""
Proven NUMT Validation: Codon Usage Bias Analysis
===================================================
Revised analysis following supervisor feedback.

KEY DIFFERENCE FROM SCRIPT 4:
  The previous analysis (Script 4) compared Authenticated ASVs against ALL
  putative Nuclear Pseudogenes (n = 26,416 unique sequences). This revision
  uses ONLY 'proven' NUMTs -- those that co-occur repeatedly alongside the
  SAME authenticated anchor ASV in >=2 independent specimens. This provides
  a much stronger biological argument because these NUMTs are validated as
  stable genomic features rather than sporadic artefacts.

Proven NUMT criteria:
  1. Classified as Nuclear_Pseudogenes in the authentication framework
  2. Read count >= MRCT (4 reads) in at least one specimen
  3. Co-occurs alongside the SAME authenticated Main_ASV in >= 2 independent
     specimens (consistent co-amplification as a stable nuclear insertion)

Analyses performed:
  1. Nucleotide composition (GC/AT content, position-specific GC)
  2. Effective Number of Codons (ENC)
  3. Per-ASV RSCU values + PCA
  4. Amino acid composition (hydrophobic, leucine, internal stops)
  5. Case study panel (3 proven NUMTs from different beetle families)

Outputs:
  - Figure 7 (revised): proven_NUMT_Validation_Figure7.png/pdf
  - Tables: proven_numt_nucleotide_composition.csv
            proven_numt_amino_acid_composition.csv
            proven_numt_rscu_per_codon.csv
            proven_numt_per_asv_rscu_matrix.csv
            proven_numt_case_studies.csv
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # data_analysis/
from config import (CSV_PATH, FASTA_PATH, NUMT_OUT, LMM_OUT,
                    ANCHOR_PAIRS_CSV, NUMT_CONCORDANT_IDS_CSV, NUMT_930_IDS_CSV,
                    LMM_ALL_METRICS_CSV, LMM_GC_CSV, LMM_READCOUNT_CSV, FAMILY_GC_CSV)

import os
import warnings
import pandas as pd
import numpy as np
from collections import Counter
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

warnings.filterwarnings('ignore')

# ============================================================
# PATHS
# ============================================================
CSV_PATH   = str(CSV_PATH)
FASTA_PATH = str(FASTA_PATH)
OUT_DATA   = str(NUMT_OUT) + '/'
OUT_FIG    = str(NUMT_OUT) + '/'   # figures saved alongside data in output/4_numt_analysis/
os.makedirs(OUT_DATA, exist_ok=True)

# ============================================================
# INVERTEBRATE MITOCHONDRIAL GENETIC CODE (NCBI Table 5)
# ============================================================
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'M','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'W','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'S','AGG':'S',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

AA_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa != '*':
        AA_CODONS.setdefault(aa, []).append(codon)

SENSE_CODONS = sorted([c for c, aa in CODON_TABLE.items() if aa != '*'])

# Broad definition for transmembrane COI protein (nonpolar + aromatic side chains)
# A,V,I,L,M,F,W,P (nonpolar aliphatic+aromatic) + G,Y,C gives ~65-68% for COI
HYDROPHOBIC_AAS = set('AVILMFWPGYC')


# ============================================================
# SEQUENCE ANALYSIS FUNCTIONS
# ============================================================

def parse_fasta(path):
    """Parse FASTA file into {id: sequence} dict."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq).upper()
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq).upper()
    return seqs


def gc_content(seq):
    """Overall GC%."""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return 100.0 * gc / len(seq) if len(seq) > 0 else np.nan


def position_gc(seq, pos, frame=1):
    """GC% at codon position pos (1, 2, or 3), respecting reading frame.
    frame=1: skip first base before counting codons (consistent with COI amplicon).
    Without frame correction, pos labels are systematically shifted.
    """
    seq = seq.upper()[frame:]  # apply same frame offset as ENC/RSCU
    bases = [seq[i + pos - 1] for i in range(0, (len(seq) // 3) * 3, 3)
             if i + pos - 1 < len(seq)]
    if not bases:
        return np.nan
    return 100.0 * sum(1 for b in bases if b in 'GC') / len(bases)


def effective_number_codons(seq, frame=1):
    """
    Effective Number of Codons (ENC) using Wright (1990).
    Lower ENC = more biased codon usage (range 20-61).
    Uses frame=1 (offset 1 bp) consistent with the COI amplicon structure.
    F values are clipped to [1/n_syn, 1.0] to prevent numerical instability.
    """
    seq = seq.upper()[frame:]
    codons = [seq[i:i+3] for i in range(0, (len(seq) // 3) * 3, 3)]
    codon_counts = Counter(c for c in codons
                           if c in CODON_TABLE and CODON_TABLE[c] != '*')
    F_by_size = {}
    for aa, syn_codons in AA_CODONS.items():
        n_syn = len(syn_codons)
        if n_syn == 1:
            continue  # single-codon families: F=1 by definition, handled in enc=2
        counts = [codon_counts.get(c, 0) for c in syn_codons]
        n = sum(counts)
        if n < 2:
            continue  # need at least 2 codon observations to estimate F
        p = [c / n for c in counts]
        sum_p2 = sum(pi ** 2 for pi in p)
        F = (n * sum_p2 - 1) / (n - 1)
        F = max(1.0 / n_syn, min(1.0, F))   # clip to valid range
        F_by_size.setdefault(n_syn, []).append(F)

    avg_F = {k: np.mean(v) for k, v in F_by_size.items()}
    enc = 2.0   # contribution from single-codon amino acid families (Met, Trp)
    for n_syn, weight in [(2, 9), (3, 1), (4, 5), (6, 3)]:
        if n_syn in avg_F and avg_F[n_syn] > 0:
            enc += weight / avg_F[n_syn]
    return float(np.clip(enc, 20.0, 61.0))


def calculate_rscu(seq, frame=1):
    """Per-sequence RSCU for all sense codons using COI reading frame (frame=1)."""
    seq = seq.upper()[frame:]
    codons = [seq[i:i+3] for i in range(0, (len(seq) // 3) * 3, 3)]
    # Count only sense codons (exclude stops)
    codon_counts = Counter(c for c in codons
                           if c in CODON_TABLE and CODON_TABLE[c] != '*')
    rscu = {}
    for aa, syn_codons in AA_CODONS.items():
        n_syn = len(syn_codons)
        total = sum(codon_counts.get(c, 0) for c in syn_codons)
        for c in syn_codons:
            rscu[c] = (codon_counts.get(c, 0) * n_syn / total) if total > 0 else np.nan
    return rscu


def translate_seq(seq, frame=1):
    """Translate using invertebrate mt code (Table 5), with reading frame offset."""
    seq = seq.upper()[frame:]
    protein = []
    for i in range(0, (len(seq) // 3) * 3, 3):
        codon = seq[i:i+3]
        protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)


def compute_sequence_metrics(seq, frame=1):
    """
    Return dict of all metrics for one sequence.
    frame=1 corresponds to the COI amplicon reading frame (empirically validated
    across 200 random sequences: 198/200 had fewest stop codons at frame=1).
    """
    seq_raw = seq.upper()
    # GC metrics computed on full sequence (frame-independent)
    gc = gc_content(seq_raw)
    protein = translate_seq(seq_raw, frame=frame)
    sense_aas = [aa for aa in protein if aa not in ('*', 'X')]
    n_aa = len(sense_aas)
    # internal stops = stops within the open reading frame (not counting terminal)
    internal_stops = protein[:-1].count('*')

    hydrophobic_pct = (100.0 * sum(1 for aa in sense_aas if aa in HYDROPHOBIC_AAS) / n_aa
                       if n_aa > 0 else np.nan)
    leucine_pct = (100.0 * sense_aas.count('L') / n_aa
                   if n_aa > 0 else np.nan)

    return {
        'GC_content':      gc,
        'AT_content':      100.0 - gc,
        'gc_pos1':         position_gc(seq_raw, 1, frame=frame),
        'gc_pos2':         position_gc(seq_raw, 2, frame=frame),
        'gc_pos3':         position_gc(seq_raw, 3, frame=frame),
        'ENC':             effective_number_codons(seq_raw, frame=frame),
        'internal_stops':  internal_stops,
        'hydrophobic_pct': hydrophobic_pct,
        'leucine_pct':     leucine_pct,
    }


def effect_size_rb(u_stat, n1, n2):
    """Rank-biserial correlation as effect size for Mann-Whitney U."""
    return 1 - (2 * u_stat) / (n1 * n2)


# ============================================================
# 1. LOAD DATA
# ============================================================
print("=" * 70)
print("PROVEN NUMT VALIDATION: CODON USAGE BIAS ANALYSIS")
print("Using only NUMTs validated by co-occurrence in >= 2 specimens")
print("=" * 70)

print("\n[1/7] Loading data...")
df = pd.read_csv(CSV_PATH, low_memory=False)
print(f"  Total records: {len(df):,}")
print(f"  ASV classification distribution:")
for cat, cnt in df['asv_classification'].value_counts().items():
    print(f"    {cat}: {cnt:,}")

print("\n  Parsing FASTA sequences...")
fasta = parse_fasta(FASTA_PATH)
print(f"  FASTA sequences: {len(fasta):,}")

# ============================================================
# 2. IDENTIFY PROVEN NUMTs
# ============================================================
print("\n[2/7] Identifying Proven NUMTs...")

numt_all = df[df['asv_classification'] == 'Nuclear_Pseudogenes'].copy()
numt_filt = numt_all[numt_all['reads'] >= 4].copy()  # above MRCT
print(f"  All Nuclear_Pseudogenes records: {len(numt_all):,}")
print(f"  Above MRCT=4: {len(numt_filt):,}  ({numt_filt['asv_id'].nunique():,} unique ASVs)")

# Co-occurrence: count unique specimens per (asv_id, Main_ASV) pair
cooccur = (numt_filt
           .groupby(['asv_id', 'Main_ASV'])['project_sample_id']
           .nunique()
           .reset_index())
cooccur.columns = ['asv_id', 'Main_ASV', 'n_specimens']

proven_pairs = cooccur[cooccur['n_specimens'] >= 2].copy()
proven_asv_ids_all = set(proven_pairs['asv_id'])

# Apply three-stage filter (manuscript-concordant pipeline):
# 1. Family-concordance: NUMT morph family == anchor morph family == NUMT phylo family == anchor phylo family → n=1,016
# 2. Pair-specific phylodist <0.50 BL: for each (NUMT, anchor) pair, min phylodist across co-occurrence specimens → n=966
# 3. Subfamily concordance: if both NUMT and anchor have subfamily info, they must match → n=956
# Pre-filtered family-concordance list; pair-specific phylodist and subfamily applied inline
_fam_concordant = pd.read_csv(str(NUMT_CONCORDANT_IDS_CSV))
_final_930 = pd.read_csv(str(NUMT_930_IDS_CSV))
proven_asv_ids_fam = proven_asv_ids_all & set(_fam_concordant['asv_id'])
proven_asv_ids = proven_asv_ids_all & set(_final_930['asv_id'])  # n=930 final dataset (10% p-dist cutoff)

print(f"\n  Proven NUMTs (same anchor, >=2 specimens):")
print(f"    All co-occurrence NUMTs: {len(proven_asv_ids_all):,}")
print(f"    After family-concordance filter (n=1,016): {len(proven_asv_ids_fam):,}")
print(f"    After pair-specific phylodist + subfamily concordance + p-dist ≤10% (n=930): {len(proven_asv_ids):,}")
print(f"    Removed (cross-family contaminants): {len(proven_asv_ids_all) - len(proven_asv_ids_fam):,}")
print(f"    Removed (pair phylodist ≥0.50 BL + cross-subfamily + p-dist >10%): {len(proven_asv_ids_fam) - len(proven_asv_ids):,}")
print(f"    Specimen co-occurrence distribution:")
for n, cnt in proven_pairs['n_specimens'].value_counts().sort_index().items():
    print(f"      {n} specimens: {cnt} NUMTs")

# Get representative records for proven NUMTs (one row per asv_id)
proven_meta = (numt_filt[numt_filt['asv_id'].isin(proven_asv_ids)]
               .sort_values('reads', ascending=False)
               .drop_duplicates(subset='asv_id', keep='first')
               .copy())
proven_meta = proven_meta.merge(proven_pairs[['asv_id','n_specimens']], on='asv_id', how='left')

# ============================================================
# 3. COMPUTE SEQUENCE METRICS
# ============================================================
print("\n[3/7] Computing sequence metrics from FASTA...")

# --- Corrected comparison (per supervisor feedback, Comment #7) ---
# The comparison must be between proven NUMTs and the authenticated ASVs
# from the SAME specimens — not all authenticated ASVs.
# Using the specific anchor ASVs that co-occur with proven NUMTs gives
# the most direct paired comparison and avoids between-specimen confounding
# from species that do not carry NUMTs at all.

# Load anchor ASV IDs directly from the pre-built pair table (canonical source).
# This ensures exact consistency with Table_S_NUMT_Anchor_Pairs_n930.csv.
# NOTE: deriving anchors inline via proven_pairs + phylodist filter inadvertently
# includes 2 extra anchors (uniq1574, uniq222) — the "bad" high-pdist partners of
# uniq398 and uniq80 whose pairs were explicitly removed in the 10% p-dist cutoff.
# Loading from the pair table is the only way to apply all four filter stages correctly.
_pair_table_n930 = pd.read_csv(str(ANCHOR_PAIRS_CSV))
anchor_asv_ids = set(_pair_table_n930['Anchor_asv_id'].unique())
print(f"  Unique anchor ASVs from Table_S_NUMT_Anchor_Pairs_n930.csv: {len(anchor_asv_ids):,}")

# Authenticated ASVs: only the matched anchor ASVs
auth_df = (df[
    (df['asv_classification'] == 'Authenticated') &
    (df['asv_id'].isin(anchor_asv_ids))
].drop_duplicates(subset='asv_id', keep='first').copy())
print(f"  Matched authenticated ASVs retained: {len(auth_df):,}  "
      f"(from specimens with proven NUMTs only)")

def compute_metrics_for_group(id_list, label):
    records = []
    missing = 0
    for asv_id in id_list:
        seq = fasta.get(asv_id)
        if seq is None or len(seq) < 300:
            missing += 1
            continue
        m = compute_sequence_metrics(seq)
        m['asv_id'] = asv_id
        records.append(m)
    print(f"  {label}: {len(records):,} sequences computed ({missing} missing/short)")
    return pd.DataFrame(records)

auth_metrics = compute_metrics_for_group(auth_df['asv_id'].tolist(), 'Authenticated')
numt_metrics = compute_metrics_for_group(list(proven_asv_ids), 'Proven NUMTs')

# Merge with metadata for proven NUMTs
numt_metrics = numt_metrics.merge(
    proven_meta[['asv_id','Main_ASV','project_sample_id','family',
                 'country','Phylogenetic_distance','percentage_reads',
                 'reads','n_specimens']],
    on='asv_id', how='left')

# Deduplicate: keep one row per unique NUMT ASV (primary anchor = highest reads)
# This ensures all analyses use exactly n=948 unique proven NUMT sequences (3-stage filter)
numt_metrics = numt_metrics.drop_duplicates('asv_id').copy()

print(f"\n  Initial sample sizes (before p-distance filter):")
print(f"    Authenticated ASVs: {len(auth_metrics):,}")
print(f"    Proven NUMTs:       {len(numt_metrics):,}  (unique ASV sequences)")

# ============================================================
# 3b. LOAD ALIGNED P-DISTANCES FROM PAIR TABLE
# ============================================================
# Use aligned (Needleman-Wunsch) p-distances from the pre-computed pair table.
# These correctly handle indels (e.g. 2SO16793 has 3-bp indel; positional p-dist=54%
# but aligned p-dist=0.72%). The pair table uses the minimum p-dist per NUMT
# (for NUMTs with multiple valid anchors, the closest anchor is used for reporting).
print("\n[3b] Loading aligned p-distances from pair table (n=930 dataset, ≤10% cutoff)...")

_pair_table = pd.read_csv(str(ANCHOR_PAIRS_CSV))
# Use minimum aligned p-dist per NUMT (closest anchor, for NUMTs with 2+ valid anchors)
_pdist_min = (_pair_table.groupby('NUMT_asv_id')['Seq_pdist_pct']
              .min().reset_index()
              .rename(columns={'NUMT_asv_id': 'asv_id', 'Seq_pdist_pct': 'numt_to_anchor_pdist_pct'}))
numt_metrics = numt_metrics.merge(_pdist_min, on='asv_id', how='left')
# Store as fraction (0–1) to match legacy column name used downstream
numt_metrics['numt_to_anchor_pdist'] = numt_metrics['numt_to_anchor_pdist_pct'] / 100.0

pdist_all = numt_metrics['numt_to_anchor_pdist'].dropna()
n_total_before = len(numt_metrics)
print(f"  Aligned p-distance loaded for {len(pdist_all):,} / {n_total_before:,} proven NUMTs")
print(f"  Full distribution: min={pdist_all.min()*100:.2f}%  "
      f"median={pdist_all.median()*100:.2f}%  "
      f"mean={pdist_all.mean()*100:.2f}%  "
      f"max={pdist_all.max()*100:.2f}%")

bins   = [0, 0.01, 0.05, 0.10, 1.0]
labels = ['<1%', '1–5%', '5–10%', '>10%']
counts = pd.cut(pdist_all, bins=bins, labels=labels, right=True).value_counts().sort_index()
print("\n  Distance bins (all NUMTs, aligned):")
for lbl, cnt in counts.items():
    print(f"    {lbl}: {cnt:,}  ({100*cnt/len(pdist_all):.1f}%)")

print(f"\n  *** Final sample sizes for ALL analyses ***")
print(f"    Authenticated ASVs (matched anchors): {len(auth_metrics):,}")
print(f"    Proven NUMTs:                         {len(numt_metrics):,}")

# ============================================================
# 4. NUCLEOTIDE COMPOSITION ANALYSIS
# ============================================================
print("\n[4/7] Nucleotide Composition Analysis...")

comp_metrics = {
    'GC_content': 'GC Content (%)',
    'AT_content': 'AT Content (%)',
    'gc_pos1':    'GC at 1st Codon Position (%)',
    'gc_pos2':    'GC at 2nd Codon Position (%)',
    'gc_pos3':    'GC at 3rd Codon Position (%)',
}

comp_results = []
for col, label in comp_metrics.items():
    a_v = auth_metrics[col].dropna()
    n_v = numt_metrics[col].dropna()
    u, p = stats.mannwhitneyu(a_v, n_v, alternative='two-sided')
    r  = effect_size_rb(u, len(a_v), len(n_v))
    comp_results.append({
        'Metric': label, 'Column': col,
        'Auth_mean': a_v.mean(), 'Auth_sd': a_v.std(), 'Auth_median': a_v.median(),
        'ProvenNUMT_mean': n_v.mean(), 'ProvenNUMT_sd': n_v.std(), 'ProvenNUMT_median': n_v.median(),
        'U_statistic': u, 'p_value': p, 'rank_biserial_r': r,
        'n_Auth': len(a_v), 'n_ProvenNUMT': len(n_v),
    })
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    print(f"  {label}:")
    print(f"    Auth:        {a_v.mean():.2f} ± {a_v.std():.2f}")
    print(f"    Proven NUMT: {n_v.mean():.2f} ± {n_v.std():.2f}   p={p:.2e}  r={r:.3f} {sig}")

comp_df = pd.DataFrame(comp_results)
comp_df.to_csv(f'{OUT_DATA}proven_numt_nucleotide_composition.csv', index=False)

# ============================================================
# 5. ENC ANALYSIS
# ============================================================
print("\n[5/7] ENC (Effective Number of Codons) Analysis...")

enc_auth = auth_metrics['ENC'].dropna()
enc_numt = numt_metrics['ENC'].dropna()
u_enc, p_enc = stats.mannwhitneyu(enc_auth, enc_numt, alternative='two-sided')
r_enc = effect_size_rb(u_enc, len(enc_auth), len(enc_numt))
sig_enc = '***' if p_enc < 0.001 else '**' if p_enc < 0.01 else '*' if p_enc < 0.05 else 'ns'
print(f"  Auth ENC:        {enc_auth.mean():.2f} ± {enc_auth.std():.2f}")
print(f"  Proven NUMT ENC: {enc_numt.mean():.2f} ± {enc_numt.std():.2f}")
print(f"  Mann-Whitney U={u_enc:.0f}, p={p_enc:.2e}, r={r_enc:.3f} {sig_enc}")

# Amino acid composition
aa_results = []
for col, label in [('hydrophobic_pct', 'Hydrophobic Amino Acids (%)'),
                   ('leucine_pct',     'Leucine (%)'),
                   ('internal_stops',  'Internal Stop Codons')]:
    a_v = auth_metrics[col].dropna()
    n_v = numt_metrics[col].dropna()
    u, p = stats.mannwhitneyu(a_v, n_v, alternative='two-sided')
    r  = effect_size_rb(u, len(a_v), len(n_v))
    aa_results.append({
        'Metric': label, 'Column': col,
        'Auth_mean': a_v.mean(), 'Auth_sd': a_v.std(),
        'ProvenNUMT_mean': n_v.mean(), 'ProvenNUMT_sd': n_v.std(),
        'U_statistic': u, 'p_value': p, 'rank_biserial_r': r,
    })
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    print(f"  {label}:")
    print(f"    Auth:        {a_v.mean():.3f} ± {a_v.std():.3f}")
    print(f"    Proven NUMT: {n_v.mean():.3f} ± {n_v.std():.3f}   p={p:.2e} {sig}")

aa_df = pd.DataFrame(aa_results)
aa_df.to_csv(f'{OUT_DATA}proven_numt_amino_acid_composition.csv', index=False)

# ============================================================
# 6. RSCU CALCULATION + PCA
# ============================================================
print("\n[6/7] RSCU Calculation and PCA...")

def compute_rscu_matrix(id_list, label, max_n=None):
    records = []
    ids_to_use = id_list if max_n is None else id_list[:max_n]
    for asv_id in ids_to_use:
        seq = fasta.get(asv_id)
        if seq is None or len(seq) < 300:
            continue
        rscu = calculate_rscu(seq)
        rscu['asv_id'] = asv_id
        rscu['group'] = label
        records.append(rscu)
    return pd.DataFrame(records)

# Use all proven NUMTs (n=923, ≤10% p-dist) and matched anchor authenticated ASVs
rscu_auth = compute_rscu_matrix(auth_df['asv_id'].tolist(), 'Authenticated')
rscu_numt = compute_rscu_matrix(list(proven_asv_ids), 'Proven NUMT')
rscu_df = pd.concat([rscu_auth, rscu_numt], ignore_index=True)
print(f"  RSCU computed: Auth={len(rscu_auth)}, Proven NUMT={len(rscu_numt)}")

rscu_df.to_csv(f'{OUT_DATA}proven_numt_per_asv_rscu_matrix.csv', index=False)

# PCA
codon_cols = [c for c in SENSE_CODONS if c in rscu_df.columns]
X = rscu_df[codon_cols].fillna(rscu_df[codon_cols].mean()).values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
pca = PCA(n_components=5)
pcs = pca.fit_transform(X_scaled)
print(f"  PCA: PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}")

# RSCU per-codon comparison (Mann-Whitney + BH correction)
auth_mask_r = rscu_df['group'] == 'Authenticated'
numt_mask_r = rscu_df['group'] == 'Proven NUMT'
auth_mean_rscu = rscu_df.loc[auth_mask_r, codon_cols].mean()
numt_mean_rscu = rscu_df.loc[numt_mask_r, codon_cols].mean()
delta_rscu     = numt_mean_rscu - auth_mean_rscu

codon_pvals = {}
for c in codon_cols:
    a_v = rscu_df.loc[auth_mask_r, c].dropna()
    n_v = rscu_df.loc[numt_mask_r, c].dropna()
    if len(a_v) > 5 and len(n_v) > 5:
        _, p = stats.mannwhitneyu(a_v, n_v, alternative='two-sided')
        codon_pvals[c] = p
    else:
        codon_pvals[c] = 1.0

pvals_arr = np.array([codon_pvals[c] for c in codon_cols])
reject, pvals_bh, _, _ = multipletests(pvals_arr, alpha=0.05, method='fdr_bh')
n_sig = int(reject.sum())
print(f"  Significant codons (BH-adjusted p<0.05): {n_sig}/{len(codon_cols)}")

rscu_comp = pd.DataFrame({
    'codon': codon_cols,
    'amino_acid': [CODON_TABLE.get(c, '?') for c in codon_cols],
    'Auth_mean_RSCU': auth_mean_rscu.values,
    'ProvenNUMT_mean_RSCU': numt_mean_rscu.values,
    'delta_RSCU_NUMT_minus_Auth': delta_rscu.values,
    'p_value_raw': pvals_arr,
    'p_adjusted_BH': pvals_bh,
    'significant_BH005': reject,
})
rscu_comp.sort_values('delta_RSCU_NUMT_minus_Auth', inplace=True)
rscu_comp.to_csv(f'{OUT_DATA}proven_numt_rscu_per_codon.csv', index=False)

delta_sorted = pd.Series(delta_rscu.values, index=codon_cols).sort_values()

# ============================================================
# 7. NUMT-TO-ANCHOR P-DISTANCE SUMMARY (aligned, n=930)
# ============================================================
# Aligned p-distances loaded from pair table in section 3b.
# All 930 NUMTs are ≤10% by design (pre-filtered).
print("\n[7/7] p-distance summary (n=930 proven NUMTs, aligned, ≤10% cutoff applied)...")

pdist_vals = numt_metrics['numt_to_anchor_pdist'].dropna()
n_computed = len(pdist_vals)
print(f"  p-distance computed for: {n_computed:,} proven NUMTs (full dataset)")
print(f"  min    = {pdist_vals.min():.4f} ({pdist_vals.min()*100:.2f}%)")
print(f"  Q25    = {pdist_vals.quantile(0.25):.4f} ({pdist_vals.quantile(0.25)*100:.2f}%)")
print(f"  median = {pdist_vals.median():.4f} ({pdist_vals.median()*100:.2f}%)")
print(f"  Q75    = {pdist_vals.quantile(0.75):.4f} ({pdist_vals.quantile(0.75)*100:.2f}%)")
print(f"  max    = {pdist_vals.max():.4f} ({pdist_vals.max()*100:.2f}%)")
print(f"  mean   = {pdist_vals.mean():.4f}  SD = {pdist_vals.std():.4f}")

# Save distribution to CSV
numt_pdist_out = numt_metrics[['asv_id','Main_ASV','family','n_specimens',
                                'percentage_reads','numt_to_anchor_pdist']].copy()
numt_pdist_out.to_csv(f'{OUT_DATA}proven_numt_anchor_pdist.csv', index=False)
print(f"\n  Saved: proven_numt_anchor_pdist.csv")

# ============================================================
# 7b. EXPORT EXTENDED LMM DATA (all metrics, per-occurrence)
# ============================================================
print("\n[7b] Exporting extended LMM data for within-specimen analysis...")

# Per-ASV metric lookup (computed from FASTA) — deduplicate by asv_id first
LMM_COLS = ['GC_content', 'gc_pos1', 'gc_pos2', 'gc_pos3',
            'ENC', 'hydrophobic_pct', 'leucine_pct']

numt_metric_lookup = (numt_metrics
                      .drop_duplicates('asv_id')
                      .set_index('asv_id')
                      [LMM_COLS]
                      .to_dict('index'))
auth_metric_lookup = (auth_metrics
                      .drop_duplicates('asv_id')
                      .set_index('asv_id')
                      [LMM_COLS]
                      .to_dict('index'))

# Proven NUMT occurrences: one row per (asv_id, specimen_id) — n=930 NUMTs (3-stage filter + p-dist ≤10%)
numt_occ = (numt_filt[numt_filt['asv_id'].isin(proven_asv_ids)]
            [['asv_id', 'project_sample_id']]
            .drop_duplicates()
            .copy())
numt_occ.columns = ['asv_id', 'specimen_id']

# Anchor ASV occurrences: one row per (asv_id, specimen_id)
# Use only the matched anchor ASVs that have valid FASTA sequences (same set as auth_metrics)
valid_anchor_ids = set(auth_metrics['asv_id'])
auth_occ = (df[(df['asv_classification'] == 'Authenticated') &
               (df['asv_id'].isin(valid_anchor_ids))]
            [['asv_id', 'project_sample_id']]
            .drop_duplicates()
            .copy())
auth_occ.columns = ['asv_id', 'specimen_id']

def attach_metrics(occ_df, lookup, group_label):
    rows = []
    for _, row in occ_df.iterrows():
        m = lookup.get(row['asv_id'])
        if m is None:
            continue  # no sequence in FASTA
        rows.append({
            'asv_id':          row['asv_id'],
            'specimen_id':     row['specimen_id'],
            'group':           group_label,
            'GC_content':      m['GC_content'],
            'gc_pos1':         m['gc_pos1'],
            'gc_pos2':         m['gc_pos2'],
            'gc_pos3':         m['gc_pos3'],
            'ENC':             m['ENC'],
            'hydrophobic_pct': m['hydrophobic_pct'],
            'leucine_pct':     m['leucine_pct'],
        })
    return pd.DataFrame(rows)

lmm_numt_df = attach_metrics(numt_occ, numt_metric_lookup, 'Proven_NUMT')
lmm_auth_df = attach_metrics(auth_occ, auth_metric_lookup, 'Authenticated')
lmm_all = pd.concat([lmm_numt_df, lmm_auth_df], ignore_index=True)

lmm_all.to_csv(str(LMM_ALL_METRICS_CSV), index=False)
print(f"  Proven NUMT rows:   {len(lmm_numt_df):,}")
print(f"  Authenticated rows: {len(lmm_auth_df):,}")
print(f"  Total rows:         {len(lmm_all):,}")
print(f"  Saved: {LMM_ALL_METRICS_CSV}")

# ============================================================
# FIGURE 7 (REVISED): 4-PANEL COMPOSITE
# ============================================================
print("\nGenerating Figure 7 (revised)...")

PALETTE = {'Authenticated': '#2196F3', 'Proven NUMT': '#E53935'}
ALPHA = 0.55

fig = plt.figure(figsize=(16, 14))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.35)

# Within-specimen LMM p-values (from additional_analyses.R, Analysis 2b)
# Used for significance labels on violin plots — within-specimen is the primary analysis
LMM_P_GC  = 0.139   # GC Content:  LMM LRT p (n=930, anchor n=330 corrected, ns)
LMM_P_ENC = 0.406   # ENC:         LMM LRT p (n=930, anchor n=330 corrected, ns)

def lmm_sig_label(p):
    return '****' if p < 0.0001 else '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

# ---- Panel A: GC Content violin ----
ax_a = fig.add_subplot(gs[0, 0])
plot_gc = pd.DataFrame({
    'GC Content (%)': pd.concat([auth_metrics['GC_content'].dropna(),
                                  numt_metrics['GC_content'].dropna()], ignore_index=True),
    'Group': (['Authenticated'] * auth_metrics['GC_content'].notna().sum() +
               ['Proven NUMT']  * numt_metrics['GC_content'].notna().sum())
})
sns.violinplot(data=plot_gc, x='Group', y='GC Content (%)', ax=ax_a,
               palette=PALETTE, inner='box', cut=0, linewidth=0.8)
ax_a.set_title('(A) GC Content', fontsize=12, fontweight='bold')
ax_a.set_xlabel('')
ax_a.set_ylabel('GC Content (%)')
ax_a.text(0.5, 0.97, lmm_sig_label(LMM_P_GC), ha='center', va='top',
          transform=ax_a.transAxes, fontsize=13, fontweight='bold', color='black')
ax_a.text(0.5, 0.03, 'within-specimen LMM', ha='center', va='bottom',
          transform=ax_a.transAxes, fontsize=7, fontstyle='italic', color='grey')
ax_a.set_ylim(bottom=15)

# ---- Panel B: ENC violin ----
ax_b = fig.add_subplot(gs[0, 1])
plot_enc = pd.DataFrame({
    'ENC': pd.concat([enc_auth, enc_numt], ignore_index=True),
    'Group': (['Authenticated'] * len(enc_auth) + ['Proven NUMT'] * len(enc_numt))
})
sns.violinplot(data=plot_enc, x='Group', y='ENC', ax=ax_b,
               palette=PALETTE, inner='box', cut=0, linewidth=0.8)
ax_b.set_title('(B) Effective Number of Codons (ENC)', fontsize=12, fontweight='bold')
ax_b.set_xlabel('')
ax_b.set_ylabel('ENC')
ax_b.text(0.5, 0.97, lmm_sig_label(LMM_P_ENC), ha='center', va='top',
          transform=ax_b.transAxes, fontsize=13, fontweight='bold', color='black')
ax_b.text(0.5, 0.03, 'within-specimen LMM', ha='center', va='bottom',
          transform=ax_b.transAxes, fontsize=7, fontstyle='italic', color='grey')
ax_b.set_ylim(bottom=15)

# ---- Panel C: PCA scatter ----
ax_c = fig.add_subplot(gs[1, 0])
for g, col in [('Authenticated', '#2196F3'), ('Proven NUMT', '#E53935')]:
    mask = (rscu_df['group'] == g).values
    ax_c.scatter(pcs[mask, 0], pcs[mask, 1],
                 c=col, alpha=ALPHA, s=12, label=g, edgecolors='none')
ax_c.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=10)
ax_c.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=10)
ax_c.set_title('(C) PCA of Relative Synonymous Codon Usage', fontsize=12, fontweight='bold')
ax_c.legend(fontsize=9, markerscale=2)

# ---- Panel D: NUMT-to-anchor p-distance histogram ----
ax_d = fig.add_subplot(gs[1, 1])
pdist_all_vals = numt_metrics['numt_to_anchor_pdist'].dropna() * 100

ax_d.hist(pdist_all_vals, bins=30, color='#E53935', alpha=0.75,
          edgecolor='white', linewidth=0.4)
ax_d.axvline(pdist_all_vals.median(), color='black', linewidth=1.2, linestyle='--',
             label=f'Median = {pdist_all_vals.median():.1f}%')
ax_d.axvline(pdist_all_vals.mean(), color='grey', linewidth=1.0, linestyle=':',
             label=f'Mean = {pdist_all_vals.mean():.1f}%')
ax_d.set_xlim(0, 30)
ax_d.set_xlabel('Pairwise p-distance to anchor ASV (%)', fontsize=10)
ax_d.set_ylabel('Number of proven NUMTs', fontsize=10)
ax_d.set_title('(D) NUMT Divergence from Anchor ASV', fontsize=11, fontweight='bold')
ax_d.legend(fontsize=8)

# ---- Figure caption info ----
fig.suptitle(
    f'Figure 7. Codon usage bias analysis comparing matched Authenticated ASVs (n = {len(auth_metrics):,}) '
    f'and Proven NUMTs (n = {len(numt_metrics):,}, ≤10% p-distance threshold)\n'
    f'Significance labels (Panels A–B) from within-specimen LMM (Table 3); '
    f'cross-specimen Mann-Whitney in Supplementary Table S8',
    fontsize=10, fontweight='bold', y=0.998)

plt.savefig(f'{OUT_FIG}proven_NUMT_Validation_Figure7.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f'{OUT_FIG}proven_NUMT_Validation_Figure7.pdf',
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: proven_NUMT_Validation_Figure7.png / .pdf")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE — KEY RESULTS")
print("=" * 70)
print(f"""
Sample Sizes:
  Authenticated ASVs (matched anchor ASVs of proven NUMTs): {len(auth_metrics):,}
  Proven NUMTs (co-occur >=2 specimens): {len(numt_metrics):,}
  NOTE: Authenticated ASVs are restricted to the direct anchor sequences
        of proven NUMTs (not all authenticated ASVs) to avoid
        between-specimen confounding from species without NUMTs.

Nucleotide Composition:""")
for r in comp_results:
    d = 'NUMT>Auth' if r['ProvenNUMT_mean'] > r['Auth_mean'] else 'Auth>NUMT'
    s = '***' if r['p_value'] < 0.001 else '**' if r['p_value'] < 0.01 else '*' if r['p_value'] < 0.05 else 'ns'
    print(f"  {r['Metric'][:30]:30s}  Auth={r['Auth_mean']:.2f}  NUMT={r['ProvenNUMT_mean']:.2f}  {d}  p={r['p_value']:.2e}  r={r['rank_biserial_r']:.3f} {s}")

print(f"""
ENC:
  Authenticated: {enc_auth.mean():.2f} ± {enc_auth.std():.2f}
  Proven NUMT:   {enc_numt.mean():.2f} ± {enc_numt.std():.2f}
  p={p_enc:.2e}, r={r_enc:.3f} {sig_enc}

Amino Acid Composition:""")
for r in aa_results:
    s = '***' if r['p_value'] < 0.001 else '**' if r['p_value'] < 0.01 else '*' if r['p_value'] < 0.05 else 'ns'
    print(f"  {r['Metric'][:35]:35s}  Auth={r['Auth_mean']:.3f}  NUMT={r['ProvenNUMT_mean']:.3f}  p={r['p_value']:.2e}  r={r['rank_biserial_r']:.3f} {s}")

print(f"""
RSCU:
  PCA PC1: {pca.explained_variance_ratio_[0]:.1%}  PC2: {pca.explained_variance_ratio_[1]:.1%}
  Significant codons (BH-FDR<0.05): {n_sig}/{len(codon_cols)}

Output files:
  {OUT_DATA}proven_numt_nucleotide_composition.csv
  {OUT_DATA}proven_numt_amino_acid_composition.csv
  {OUT_DATA}proven_numt_rscu_per_codon.csv
  {OUT_DATA}proven_numt_per_asv_rscu_matrix.csv
  {OUT_DATA}proven_numt_case_studies.csv
  {OUT_FIG}proven_NUMT_Validation_Figure7.png/pdf
""")
print("Done!")
