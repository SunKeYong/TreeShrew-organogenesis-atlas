# Code index and contributors

Code from Keyong Sun and Ou Lin (欧霖) is combined by figure. These tables list all 79 supplied code/configuration files. Original filenames are retained. The source-to-figure mapping is available in [FILE_MAP.tsv](FILE_MAP.tsv).

## Figure 1: Atlas construction and annotation

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig1_AllCells_UMAP.ipynb](../Figure_1/Fig1_AllCells_UMAP.ipynb) | Jupyter notebook | Keyong Sun |
| [Fig1_CellType_define.R](../Figure_1/Fig1_CellType_define.R) | R | Keyong Sun |
| [Fig1_Process_fastq_files.sh](../Figure_1/Fig1_Process_fastq_files.sh) | Shell | Keyong Sun |
| [Fig1_SubLibrary_QC.ipynb](../Figure_1/Fig1_SubLibrary_QC.ipynb) | Jupyter notebook | Keyong Sun |

## Figure 2: cCRE catalogue and regulatory topics

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig2_Motif_Enrichment.R](../Figure_2/Fig2_Motif_Enrichment.R) | R | Keyong Sun |
| [Fig2_PeakCallling_Embedding.ipynb](../Figure_2/Fig2_PeakCallling_Embedding.ipynb) | Jupyter notebook | Keyong Sun |
| [Fig2_VISTA.R](../Figure_2/Fig2_VISTA.R) | R | Keyong Sun |
| [Fig2_cCREs_UMAP.R](../Figure_2/Fig2_cCREs_UMAP.R) | R | Keyong Sun |
| [Fig2_deeptools_and_liftover_for_GREAT_analysis.sh](../Figure_2/Fig2_deeptools_and_liftover_for_GREAT_analysis.sh) | Shell | Keyong Sun |
| [Fig2_liftover.sh](../Figure_2/Fig2_liftover.sh) | Shell | Keyong Sun |
| [pycisTopic/create_cistopic.py](../Figure_2/pycisTopic/create_cistopic.py) | Python | Ou Lin |
| [pycisTopic/plot_heatmap.R](../Figure_2/pycisTopic/plot_heatmap.R) | R | Ou Lin |
| [pycisTopic/transferModel_annotation_runUmap.ipynb](../Figure_2/pycisTopic/transferModel_annotation_runUmap.ipynb) | Jupyter notebook | Ou Lin |

## Figure 3: Sequence-to-function models

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig3_PeakCallling_For_ChromBPNet.ipynb](../Figure_3/Fig3_PeakCallling_For_ChromBPNet.ipynb) | Jupyter notebook | Keyong Sun |
| [Motif_cooccurrence/motif_cooccurrence_allpeaks.ipynb](../Figure_3/Motif_cooccurrence/motif_cooccurrence_allpeaks.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_cooccurrence/plot_cooccurrence_allpeaks.R](../Figure_3/Motif_cooccurrence/plot_cooccurrence_allpeaks.R) | R | Ou Lin |
| [QC_metrics/get_observed_counts.py](../Figure_3/QC_metrics/get_observed_counts.py) | Python | Ou Lin |
| [QC_metrics/get_observed_counts_rep.py](../Figure_3/QC_metrics/get_observed_counts_rep.py) | Python | Ou Lin |
| [QC_metrics/interpretation.R](../Figure_3/QC_metrics/interpretation.R) | R | Ou Lin |
| [QC_metrics/make_counts_scatter_plots.ipynb](../Figure_3/QC_metrics/make_counts_scatter_plots.ipynb) | Jupyter notebook | Ou Lin |
| [QC_metrics/random_split_bed.sh](../Figure_3/QC_metrics/random_split_bed.sh) | Shell | Ou Lin |
| [model_training_and_analysis/MotifPlotting.ipynb](../Figure_3/model_training_and_analysis/MotifPlotting.ipynb) | Jupyter notebook | Ou Lin |
| [model_training_and_analysis/Train_bias_model/make_splits.ipynb](../Figure_3/model_training_and_analysis/Train_bias_model/make_splits.ipynb) | Jupyter notebook | Ou Lin |
| [model_training_and_analysis/Train_bias_model/preprocess_trainBias.sh](../Figure_3/model_training_and_analysis/Train_bias_model/preprocess_trainBias.sh) | Shell | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_0.json](../Figure_3/model_training_and_analysis/Train_bias_model/splits/fold_0.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_1.json](../Figure_3/model_training_and_analysis/Train_bias_model/splits/fold_1.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_2.json](../Figure_3/model_training_and_analysis/Train_bias_model/splits/fold_2.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_3.json](../Figure_3/model_training_and_analysis/Train_bias_model/splits/fold_3.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_4.json](../Figure_3/model_training_and_analysis/Train_bias_model/splits/fold_4.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/inputs.yaml](../Figure_3/model_training_and_analysis/Train_use_snakefile_example/inputs.yaml) | YAML | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/run.sh](../Figure_3/model_training_and_analysis/Train_use_snakefile_example/run.sh) | Shell | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/setup_environment.sh](../Figure_3/model_training_and_analysis/Train_use_snakefile_example/setup_environment.sh) | Shell | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/snakefile](../Figure_3/model_training_and_analysis/Train_use_snakefile_example/snakefile) | Snakemake | Ou Lin |
| [model_training_and_analysis/finemo.sh](../Figure_3/model_training_and_analysis/finemo.sh) | Shell | Ou Lin |
| [model_training_and_analysis/modicso.sh](../Figure_3/model_training_and_analysis/modicso.sh) | Shell | Ou Lin |
| [plot_predicted/predict_prebw_marker_region.sh](../Figure_3/plot_predicted/predict_prebw_marker_region.sh) | Shell | Ou Lin |
| [plot_predicted/tangermeme.markerpeak.ipynb](../Figure_3/plot_predicted/tangermeme.markerpeak.ipynb) | Jupyter notebook | Ou Lin |
| [solo_motif_footprinting/test_solo_motif_footprinting_all_celltypes_complete.py](../Figure_3/solo_motif_footprinting/test_solo_motif_footprinting_all_celltypes_complete.py) | Python | Ou Lin |

## Figure 4: Cross-species cCRE comparisons

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig4_Figure_plot1.R](../Figure_4/Fig4_Figure_plot1.R) | R | Keyong Sun |
| [Fig4_Figure_plot2.R](../Figure_4/Fig4_Figure_plot2.R) | R | Keyong Sun |
| [Fig4_Human_Peak_Process.R](../Figure_4/Fig4_Human_Peak_Process.R) | R | Keyong Sun |
| [Fig4_Mouse_MOPA_Peak_Process.R.R](../Figure_4/Fig4_Mouse_MOPA_Peak_Process.R.R) | R | Keyong Sun |
| [Fig4_TE.sh](../Figure_4/Fig4_TE.sh) | Shell | Keyong Sun |
| [Fig4_TE_TreeShrew_UCSC.sh](../Figure_4/Fig4_TE_TreeShrew_UCSC.sh) | Shell | Keyong Sun |
| [Fig4_TE_mm10_Dfam_3.8.sh](../Figure_4/Fig4_TE_mm10_Dfam_3.8.sh) | Shell | Keyong Sun |
| [Fig4_Treeshrew_Peak_Process.R](../Figure_4/Fig4_Treeshrew_Peak_Process.R) | R | Keyong Sun |
| [Fig4_liftover.sh](../Figure_4/Fig4_liftover.sh) | Shell | Keyong Sun |

## Figure 5: Immune-cell comparisons

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig5_Immune_Concervation_Peaks.R](../Figure_5/Fig5_Immune_Concervation_Peaks.R) | R | Keyong Sun |
| [Fig5_Immune_Conservation_liftover.sh](../Figure_5/Fig5_Immune_Conservation_liftover.sh) | Shell | Keyong Sun |
| [Fig5_Immune_Human.R](../Figure_5/Fig5_Immune_Human.R) | R | Keyong Sun |
| [Fig5_Immune_Plot.R](../Figure_5/Fig5_Immune_Plot.R) | R | Keyong Sun |
| [Fig5_Immune_Plot2.R](../Figure_5/Fig5_Immune_Plot2.R) | R | Keyong Sun |
| [Fig5_Immune_Treeshrew.R](../Figure_5/Fig5_Immune_Treeshrew.R) | R | Keyong Sun |
| [cross_species_model_prediction/B_cell_Ery_ortho_predict/tangermeme.markerpeak.ipynb](../Figure_5/cross_species_model_prediction/B_cell_Ery_ortho_predict/tangermeme.markerpeak.ipynb) | Jupyter notebook | Ou Lin |
| [cross_species_model_prediction/Hs_TS_immune_prebw.sh](../Figure_5/cross_species_model_prediction/Hs_TS_immune_prebw.sh) | Shell | Ou Lin |
| [cross_species_model_prediction/filter_peaks_not_expand2114bp.ipynb](../Figure_5/cross_species_model_prediction/filter_peaks_not_expand2114bp.ipynb) | Jupyter notebook | Ou Lin |
| [cross_species_model_prediction/get_observed_counts.py](../Figure_5/cross_species_model_prediction/get_observed_counts.py) | Python | Ou Lin |
| [cross_species_model_prediction/get_predicted_counts.py](../Figure_5/cross_species_model_prediction/get_predicted_counts.py) | Python | Ou Lin |
| [motif_cooccurrence/motif_cooccurrence_allpeaks_OBO.ipynb](../Figure_5/motif_cooccurrence/motif_cooccurrence_allpeaks_OBO.ipynb) | Jupyter notebook | Ou Lin |
| [motif_cooccurrence/plot_cooccurrence_all_peaks_OBO.R](../Figure_5/motif_cooccurrence/plot_cooccurrence_all_peaks_OBO.R) | R | Ou Lin |

## Figure 6: Enhancer candidate analyses

| File | Type | Code contributor |
| --- | --- | --- |
| [Motif_screening/01Motif_cooccurrence/01a-liftover.HsToTS.sh](../Figure_6/Motif_screening/01Motif_cooccurrence/01a-liftover.HsToTS.sh) | Shell | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/01b-liftover.sh](../Figure_6/Motif_screening/01Motif_cooccurrence/01b-liftover.sh) | Shell | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/01c-filter_TS_bed.sh](../Figure_6/Motif_screening/01Motif_cooccurrence/01c-filter_TS_bed.sh) | Shell | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/01d-make_predict_region_bed.sh](../Figure_6/Motif_screening/01Motif_cooccurrence/01d-make_predict_region_bed.sh) | Shell | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/02-finemo.2.sh](../Figure_6/Motif_screening/01Motif_cooccurrence/02-finemo.2.sh) | Shell | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/03a-Hs_motif_cooccurrence.ipynb](../Figure_6/Motif_screening/01Motif_cooccurrence/03a-Hs_motif_cooccurrence.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/03b-TS_motif_cooccurrence.ipynb](../Figure_6/Motif_screening/01Motif_cooccurrence/03b-TS_motif_cooccurrence.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/03c-plot_vista_motif_density_per_1000bp.ipynb](../Figure_6/Motif_screening/01Motif_cooccurrence/03c-plot_vista_motif_density_per_1000bp.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_screening/01Motif_cooccurrence/03d-plot_motif_cooccurrence.R](../Figure_6/Motif_screening/01Motif_cooccurrence/03d-plot_motif_cooccurrence.R) | R | Ou Lin |
| [Motif_screening/02Filter_marker_peak_with_motif/TS_C3_1_motif_cooccurrence.ipynb](../Figure_6/Motif_screening/02Filter_marker_peak_with_motif/TS_C3_1_motif_cooccurrence.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_screening/02Filter_marker_peak_with_motif/TS_C3_2_motif_cooccurrence.ipynb](../Figure_6/Motif_screening/02Filter_marker_peak_with_motif/TS_C3_2_motif_cooccurrence.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_screening/03Select_peaks_with_high_accessibility/plot_bw_region.R](../Figure_6/Motif_screening/03Select_peaks_with_high_accessibility/plot_bw_region.R) | R | Ou Lin |
| [plot_predicted_enhancer/TS_enhancers_for_validation/tangermeme.markerpeak.ipynb](../Figure_6/plot_predicted_enhancer/TS_enhancers_for_validation/tangermeme.markerpeak.ipynb) | Jupyter notebook | Ou Lin |
| [plot_predicted_enhancer/hs3087_ortho_across_mammals/tangermeme_custom_sequence.ipynb](../Figure_6/plot_predicted_enhancer/hs3087_ortho_across_mammals/tangermeme_custom_sequence.ipynb) | Jupyter notebook | Ou Lin |

## Shared environments

| File | Code contributor |
| --- | --- |
| [environments/chrombpnet/chrombpnetv8.yaml](../environments/chrombpnet/chrombpnetv8.yaml) | Ou Lin |
| [environments/finemo/Finemo.yaml](../environments/finemo/Finemo.yaml) | Ou Lin |
| [environments/modisco/modisco-lite.yaml](../environments/modisco/modisco-lite.yaml) | Ou Lin |
| [environments/pycistopic/pycisTopic.yaml](../environments/pycistopic/pycisTopic.yaml) | Ou Lin |
| [environments/tangermeme/tangermeme.yaml](../environments/tangermeme/tangermeme.yaml) | Ou Lin |
