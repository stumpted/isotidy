# Core IRMS Functions

Universal, instrument-agnostic processing functions organized by category.

## Files

| File | Purpose |
|------|---------|
| `core_00_import.R` | Convert raw data to canonical format using adapter functions. |
| `core_01_clean.R` | Remove test runs and filter by peak/compound. Identify blanks. |
| `core_02_standards.R` | Identify standards and merge with reference database. |
| `core_03_mass.R` | Build and apply mass calibration (optional, instrument-dependent). |
| `core_04_blank.R` | Calculate blank statistics and apply blank correction (multiple methods). |
| `core_05_linearity.R` | Build and apply linearity correction (optional, instrument-dependent). |
| `core_06_scale.R` | Build and apply scale calibration (VPDB, etc.) - converts to reference scale. |
| `core_07_output.R` | Prepare output dataframe and generate quality summaries. |

## Loading

Load all files:
```r
source("scripts/core_load_all.R")
```

Load individual file:
```r
source("scripts/core/core_04_blank.R")
```

## Dependencies

- tidyverse (used by all files)
- No interdependencies between core files (safe to load selectively)

## Usage

All functions follow the same pattern:
- Input: canonical dataframe + configuration
- Output: modified dataframe + optional model objects
- Designed to chain together in pipeline sequence

See EA pipeline (`ea_pipeline.R`) for example usage.
