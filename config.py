# config.py — path configuration for all scripts
# Auto-detects data_analysis/ from this file's location. No manual edits needed.
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent   # = data_analysis/

# ── Inputs (large files; download from Dryad before running) ──────────
INPUT_DIR        = ANALYSIS_DIR / "input" / "raw_data"
CSV_PATH         = INPUT_DIR / "ASV_Authentication_Results_030925.csv"
FASTA_PATH       = INPUT_DIR / "ASV_sequences_64544.fasta"
FASTA_PATH_66595 = INPUT_DIR / "ASV_sequences_66595.fasta"

# ── Reference data (small; included in repo) ──────────────────────────
REF_DIR          = ANALYSIS_DIR / "input" / "reference"
ANCHOR_PAIRS_CSV        = REF_DIR / "Table_S_NUMT_Anchor_Pairs_n930.csv"
NUMT_CONCORDANT_IDS_CSV = REF_DIR / "proven_numt_family_concordant_ids.csv"
NUMT_930_IDS_CSV        = REF_DIR / "proven_numt_930_ids.csv"

# ── Outputs (numbered to match scripts) ───────────────────────────────
OUTPUT_DIR = ANALYSIS_DIR / "output"
SEQ_OUT    = OUTPUT_DIR / "1_sequences_analysis"
CLASS_OUT  = OUTPUT_DIR / "2_classification_analysis"
STATS_OUT  = OUTPUT_DIR / "3_statistical_analysis"
NUMT_OUT   = OUTPUT_DIR / "4_numt_analysis"
LMM_OUT    = OUTPUT_DIR / "5_lmm_analysis"
DADA2_OUT  = OUTPUT_DIR / "6_dada2_comparison"

# ── Cross-script intermediates (script 4 produces → script 5 reads) ───
LMM_ALL_METRICS_CSV = NUMT_OUT / "lmm_numt_all_metrics_data.csv"
LMM_GC_CSV          = NUMT_OUT / "lmm_numt_gc_data.csv"
LMM_READCOUNT_CSV   = LMM_OUT  / "lmm_readcount_data.csv"
FAMILY_GC_CSV       = LMM_OUT  / "family_gc_content.csv"
