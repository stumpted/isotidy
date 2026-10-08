# R/ layout

Updated 2026-10-06

One file per topic. Files are named `<area>-<topic>.R`, and the prefix groups them in directory listings (R packages require a flat `R/` folder, so there are no subfolders).

## Files

### Configuration and import

| File | Functions |
|---|---|
| `config-schema.R` | `load_peripheral_schemas`, `get_peripheral_schema`, `validate_peripheral_schemas`, `validate_peripheral_schema` |
| `config-experiment.R` | `load_IRMS_config`, `validate_config`, `create_config_template`, `get_config_element`, `get_isotope_system` |
| `normalize.R` | `normalize_data` (generic entry point: raw data to canonical format via an adapter function) |

### Processing stages (peripheral-agnostic)

| File | Functions |
|---|---|
| `clean.R` | `clean_data` |
| `standards.R` | `load_lab_standards`, `identify_standards`, `extract_standards` |
| `correct-mass.R` | `prepare_mass_standards`, `build_mass_model`, `apply_mass_correction` |
| `correct-blank.R` | `identify_blanks`, `calculate_blank_stats`, `apply_blank_correction` |
| `correct-linearity.R` | `build_linearity_model`, `apply_linearity_correction` |
| `correct-scale.R` | `build_scale_model`, `apply_scale_correction` |
| `output.R` | `prepare_output`, `summarize_calibrations`, `print_results_summary` |

### EA-IRMS (peripheral-specific)

| File | Functions |
|---|---|
| `ea-adapter.R` | `adapt_ea_data` and its helpers: `validate_ea_raw_data`, `apply_ea_validation_rules`, `add_ea_specific_columns`, `validate_ea_canonical_data` |
| `ea-pipeline.R` | `process_ea` |
| `ea-save.R` | `save_ea_results` |

### GC-IRMS (peripheral-specific)

| File | Functions |
|---|---|
| `gc-adapter.R` | `adapt_gc_data`, GC method-based peak assignment, external-standard correction |
| `gc-derivatization.R` | Estimate date- and compound-specific NACME values from AA mix runs and correct sample peaks |
| `gc-pipeline.R` | `process_gc`: adapt, assign compounds, apply F8 offsets, and apply NACME correction when carbon-count metadata is available |
| `gc-save.R` | `save_gc_results`: corrected peak CSV, offset/derivatization diagnostics, and HTML report |

The GC workflow can now apply F8 external-standard correction followed by the empirical NACME derivatization correction. It does not apply a NorLeu internal-standard correction.

### Package-level

| File | Purpose |
|---|---|
| `isotidy-package.R` | Package docs and the only imports: `%>%` (dplyr), `.data` and `%||%` (rlang) |

## Pipeline order

`process_ea()` runs the stages in this order:

1. `adapt_ea_data` converts raw EA data to the canonical format
2. `clean_data`
3. `identify_standards` and `identify_blanks`
4. `extract_standards`
5. Mass: `prepare_mass_standards`, `build_mass_model`, `apply_mass_correction`
6. Blank: `calculate_blank_stats`, `apply_blank_correction`
7. Linearity: `build_linearity_model`, `apply_linearity_correction`
8. Scale: `build_scale_model`, `apply_scale_correction`
9. `prepare_output`

Standards are re-extracted after each correction so later models see the corrected deltas.

## Adding a peripheral

Add a schema entry for the peripheral, then create `<peripheral>-adapter.R`, `<peripheral>-pipeline.R` and `<peripheral>-save.R`, mirroring the `ea-*` files. The adapter must produce the canonical columns that `normalize_data()` checks for. The stage files above are meant to be reused as-is.

## Conventions

- **Qualify calls** as `pkg::fn()`. Do not use `library()`. Operators and the `.data` pronoun are the only imports, and they live in `isotidy-package.R`.
- **Column names in dplyr verbs** use `.data$col`, or `all_of("col")` in tidyselect. Do not use bare names.
- **New dependencies** go in `Imports` in `DESCRIPTION`. `rmarkdown` is in `Suggests`.
- **Comments** are single `#` lines that explain intent or why. Skip dividers and comments that restate the next line.
- **No top-level code.** Files contain only function definitions, so load order doesn't matter. If that changes, add a `Collate` field.
- **Docs and exports** come from roxygen. After changing `@export` or `@importFrom` tags, run `devtools::document()` to regenerate `NAMESPACE`. Only functions tagged `@export` are public. Currently that is the config, schema, adapter and `save_ea_results` functions. The stage functions and `process_ea` are internal.
