# Analysis guide

Code contributed by Keyong Sun and Ou Lin (欧霖) is combined under Figure_1 through Figure_6. Use this as a navigation map. Steps describe how the supplied scripts relate to one another; intermediate data generation is not completely automated by this repository.

## Figure 1: atlas construction and cell annotation

1. [`Fig1_Process_fastq_files.sh`](../Figure_1/Fig1_Process_fastq_files.sh): trim reads, correct barcodes, attach barcodes, align reads and construct fragment files. The script references additional demultiplexing helpers and barcode tables that must be supplied separately.
2. [`Fig1_SubLibrary_QC.ipynb`](../Figure_1/Fig1_SubLibrary_QC.ipynb): calculate sublibrary QC metrics, including TSS enrichment and doublet scores, and export cell metadata.
3. [`Fig1_AllCells_UMAP.ipynb`](../Figure_1/Fig1_AllCells_UMAP.ipynb): build tile matrices, spectral embeddings, UMAPs and Leiden clusters.
4. [`Fig1_CellType_define.R`](../Figure_1/Fig1_CellType_define.R): cell-type annotation examples using ArchR and reference information.

## Figure 2: cCRE catalogue and regulatory topics

- [`Fig2_PeakCallling_Embedding.ipynb`](../Figure_2/Fig2_PeakCallling_Embedding.ipynb): cell-type peak calling, merged peak sets, CPM bigWig export and pseudo-cell peak embeddings.
- [`Figure_2/pycisTopic/`](../Figure_2/pycisTopic/): construct a cisTopic object, run Mallet topic models, inspect/annotate models and plot topic heatmaps. Start with `create_cistopic.py`, then the relevant sections of `transferModel_annotation_runUmap.ipynb` and `plot_heatmap.R` after preparing their inputs.
- [`Fig2_cCREs_UMAP.R`](../Figure_2/Fig2_cCREs_UMAP.R) and [`Fig2_Motif_Enrichment.R`](../Figure_2/Fig2_Motif_Enrichment.R): cCRE visualisation and motif enrichment.
- [`Fig2_liftover.sh`](../Figure_2/Fig2_liftover.sh), [`Fig2_deeptools_and_liftover_for_GREAT_analysis.sh`](../Figure_2/Fig2_deeptools_and_liftover_for_GREAT_analysis.sh) and [`Fig2_VISTA.R`](../Figure_2/Fig2_VISTA.R): coordinate conversion, signal summaries and external annotation analyses.

## Figure 3: sequence-to-function models

1. [`Fig3_PeakCallling_For_ChromBPNet.ipynb`](../Figure_3/Fig3_PeakCallling_For_ChromBPNet.ipynb): prepare per-cell-type peak sets for model training.
2. [`Train_bias_model/`](../Figure_3/model_training_and_analysis/Train_bias_model/): chromosome splits and bias-model training examples. Five supplied split JSON files are retained in `splits/`.
3. [`Train_use_snakefile_example/`](../Figure_3/model_training_and_analysis/Train_use_snakefile_example/): example model-training configuration and Snakefile. Edit `inputs.yaml`, environment activation and resource settings before a dry run. Relative helper-script paths require an appropriate working directory.
4. [`QC_metrics/`](../Figure_3/QC_metrics/): observed counts, pseudo-replicate counts and model QC plots.
5. [`model_training_and_analysis/`](../Figure_3/model_training_and_analysis/), [`plot_predicted/`](../Figure_3/plot_predicted/), [`Motif_cooccurrence/`](../Figure_3/Motif_cooccurrence/) and [`solo_motif_footprinting/`](../Figure_3/solo_motif_footprinting/): model interpretation, TF-MoDISco, Fi-NeMo, prediction visualisation, motif co-occurrence and footprinting examples.

## Figure 4: cross-species cCRE comparisons

- Process peak annotations using `Fig4_Treeshrew_Peak_Process.R`, `Fig4_Human_Peak_Process.R` and `Fig4_Mouse_MOPA_Peak_Process.R.R` in [`Figure_4/`](../Figure_4/).
- [`Fig4_liftover.sh`](../Figure_4/Fig4_liftover.sh) provides reciprocal coordinate-conversion examples. Adapt source/target assemblies and chromosome conventions for each comparison.
- `Fig4_TE.sh`, `Fig4_TE_TreeShrew_UCSC.sh` and `Fig4_TE_mm10_Dfam_3.8.sh` prepare repeat annotations.
- `Fig4_Figure_plot1.R` and `Fig4_Figure_plot2.R` generate comparative plots from analysis objects created in preceding steps. These plot scripts assume objects are already present in the R session.

## Figure 5: immune comparisons

- [`Figure_5/`](../Figure_5/): `Fig5_Immune_Conservation_liftover.sh` and `Fig5_Immune_Concervation_Peaks.R` prepare comparative peak sets; `Fig5_Immune_Human.R` and `Fig5_Immune_Treeshrew.R` perform species-specific analyses; `Fig5_Immune_Plot.R` and `Fig5_Immune_Plot2.R` generate plots.
- [`Figure_5/`](../Figure_5/): cross-species model predictions, observed/predicted count extraction and motif co-occurrence analyses, including B-cell and erythroid examples.

## Figure 6: enhancer candidate analyses

- [`Motif_screening/`](../Figure_6/Motif_screening/): follow the numbered subdirectories for coordinate conversion and motif co-occurrence, marker-peak filtering, and accessibility-based inspection.
- [`plot_predicted_enhancer/`](../Figure_6/plot_predicted_enhancer/): sequence prediction examples for tree shrew validation candidates and an orthologous region across mammals.

See the [complete code index](CODE_INDEX.md) for every supplied source and configuration file.
