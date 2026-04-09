# 00_config.R — Shared Configuration

# --- Paths -------------------------------------
DATA_DIR       <- "data"
RESULTS_DIR    <- "results"
RAW_DATA_FILE  <- file.path(DATA_DIR, "raw.rds")

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# --- Disease Stage Order ----------------------
DISEASE_STAGES <- c("NORM", "HO", "STEA", "NASH", "FIB", "CIRR")
