# Convenience script to source all core functions at once
#
# Usage:
#   source("scripts/core_load_all.R")
#
# This loads all 8 core function files:
#   - core_00_import.R
#   - core_01_clean.R
#   - core_02_standards.R
#   - core_03_mass.R
#   - core_04_blank.R
#   - core_05_linearity.R
#   - core_06_scale.R
#   - core_07_output.R
#
# ============================================================================

# Determine the directory where this script is located
core_dir <- dirname(sys.frame(1)$ofile %||% "scripts/core")

# If running from scripts directory, assume core is a subdirectory
if (!dir.exists(core_dir)) {
  core_dir <- file.path("scripts", "core")
}

# If core directory doesn't exist, try parent directory
if (!dir.exists(core_dir)) {
  core_dir <- file.path("..", "core")
}

# List all core files
core_files <- list.files(
  core_dir,
  pattern = "^core_[0-9]{2}.*\\.R$",
  full.names = TRUE
)

# Source them in order
if (length(core_files) == 0) {
  warning(paste("No core files found in:", core_dir))
  warning("Make sure core/*.R files are in scripts/core/ or adjust path")
} else {
  for (file in core_files) {
    source(file, local = FALSE)
  }
  
  cat("✓ Loaded", length(core_files), "core function files\n")
  cat("  Directory:", core_dir, "\n")
  cat("  Files:\n")
  for (file in core_files) {
    cat("    -", basename(file), "\n")
  }
  cat("\n")
}

