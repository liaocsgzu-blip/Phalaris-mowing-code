# Nature Communications code package
## Figures 1-5 and Supplementary Figures S1-S3

This archive contains the English-only R scripts corresponding to the final main and supplementary figures.

## Code-to-figure mapping

1. `01_Fig1A_long_term_trajectories.R` - Figure 1A
2. `02_Fig1BC_yield_predictor_importance.R` - Figure 1B-C
3. `03_Fig2AB_random_forest_Boruta_key_taxa.R` - Figure 2A-B
4. `04_Fig2C_key_taxa_profile_contribution.R` - Figure 2C
5. `05_Fig3A-D_RT_taxa_regressions.R` - Figure 3A-D
6. `06_Fig3E_RS_RT_BS_association_network.R` - Figure 3E
7. `07_Fig3FG_root_metabolite_pathway_analysis.R` - Figure 3F-G
8. `08_Fig4_and_FigS2_RT_functional_coexpression_PLS.R` - Figure 4 and Figure S2A-E
9. `09_Fig5_SEM_and_shared_variance.R` - Figure 5A-D
10. `10_FigS1_NMDS_and_ternary_enrichment.R` - Figure S1A-B
11. `12_FigS3_extended_yield_path_sensitivity_SEM.R` - Figure S3A-B

## Standardized input names

The code archive does not include source data. The upload-ready scripts use the following standardized English input names where applicable:

- `main_plant_soil_data.xlsx`
- `species_abundance.csv`
- `metabolomics_positive.xlsx`
- `metabolomics_negative.xlsx`
- `sample_mapping.csv`
- `metagenome_selected_functions.csv`
- `metatranscriptome_selected_functions.csv`
- `rt_species_abundance.csv`
- `plant_traits_yield.csv`
- `root_metabolome_summary.csv`
- `carbon_use_rs_bs.csv`
- `genus_abundance.csv`
- `nmds_coordinates.csv`

For `main_plant_soil_data.xlsx`, the standardized sheet names used in the scripts are `Plant`, `Soil`, and `CarbonUse`.

The statistical analysis logic, model structures, thresholds, random seeds, standardization procedures, and feature definitions were retained from the author-supplied scripts. Packaging edits were restricted to English-language standardization, portable file handling, final figure numbering, and removal of legacy output-name conflicts.

Source data are not included in this archive.
