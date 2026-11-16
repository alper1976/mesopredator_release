Minimal input manifest
- Metadata_TnT - metadata_main.csv and species_table_n.csv: source Mesopredator_check.zip / Data/Inputs.
- område og art_data.csv: source Mesopredator_check.zip / Data/Inputs, renamed from ZIP filename encoding to the exact filename required by the master script. File bytes unchanged.
- summary_17652/19661/19662/19663.csv: source Mesopredator_check.zip / ices_cod_age_diagnostic; placed in legacy fallback path.
- DATRAS_age_weight_summary.csv: user-provided DATRAS summary.
- Original_vs_Master_Methods_References.docx: latest coauthor methods/reference change log.
- Master script: v1.4 source with ONLY age-1-to-age-2 survival assignment and comment updated.
- Portable paths: root resolved from script location; all local input/output paths inside package.\n- Root detection now also supports RStudio interactive line execution (active saved editor or package working directory).
