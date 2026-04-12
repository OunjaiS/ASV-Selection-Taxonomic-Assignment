# config.R — path configuration for R scripts
# Auto-detects data_analysis/ from this file's location. No manual edits needed.

# If calling script already set ANALYSIS_DIR, use it; otherwise auto-detect.
if (!exists("ANALYSIS_DIR")) {
  config_dir <- tryCatch({
    # Works when sourced interactively or from source()
    calls <- sys.calls()
    src_calls <- Filter(function(x) deparse(x[[1]]) == "source", rev(calls))
    if (length(src_calls) > 0) {
      p <- tryCatch(as.character(src_calls[[1]][[2]]), error = function(e) NULL)
      if (!is.null(p) && nchar(p) > 0) dirname(normalizePath(p)) else getwd()
    } else {
      dirname(normalizePath(sys.frame(sys.nframe())$ofile))
    }
  }, error = function(e) getwd())
  ANALYSIS_DIR <- normalizePath(config_dir)
}

# ── Inputs ────────────────────────────────────────────────────────────
INPUT_DIR <- file.path(ANALYSIS_DIR, "input", "raw_data")
REF_DIR   <- file.path(ANALYSIS_DIR, "input", "reference")

# ── Outputs (numbered to match scripts) ───────────────────────────────
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "output")
NUMT_OUT   <- file.path(OUTPUT_DIR, "4_numt_analysis")
LMM_OUT    <- file.path(OUTPUT_DIR, "5_lmm_analysis")

# ── Cross-script intermediates ────────────────────────────────────────
LMM_ALL_METRICS_CSV <- file.path(NUMT_OUT, "lmm_numt_all_metrics_data.csv")
LMM_GC_CSV          <- file.path(NUMT_OUT, "lmm_numt_gc_data.csv")
LMM_READCOUNT_CSV   <- file.path(LMM_OUT,  "lmm_readcount_data.csv")
FAMILY_GC_CSV       <- file.path(LMM_OUT,  "family_gc_content.csv")

cat(sprintf("[config.R] ANALYSIS_DIR = %s\n", ANALYSIS_DIR))
