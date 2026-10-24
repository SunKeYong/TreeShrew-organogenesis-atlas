# Code provenance and acknowledgements

## Source collections

Keyong Sun supplied the `sky_Code/` collection as his curated analysis code. The `Fig2_pycistopic/`, `Fig3_models/`, `Fig5_immune_cells/` and `Fig6_enhancers/` collections were supplied as collaborator-curated code. `environments/` contains the accompanying environment configurations.

This separation records the source collections provided for upload. It does not assert sole authorship or establish a complete author list. Collaborator names, roles and publication citation can be added once confirmed. Inline comments, acknowledgements, original filenames and references have been retained.

## Third-party references retained in the supplied code

| Reference | Files that acknowledge or refer to it |
| --- | --- |
| [Lan-lab/SPATAC-seq](https://github.com/Lan-lab/SPATAC-seq/) | `sky_Code/Fig1_Process_fastq_files.sh` |
| [GreenleafLab/HDMA ChromBPNet use cases](https://github.com/GreenleafLab/HDMA/blob/main/code/05-misc/04-ChromBPNet_use_cases.ipynb) | Prediction notebooks under `Fig3_models/plot_predicted/`, `Fig5_immune_cells/cross_species_model_prediction/` and `Fig6_enhancers/plot_predicted_enhancer/` |
| [kundajelab/chrombpnet-figures](https://github.com/kundajelab/chrombpnet-figures/tree/main/src/figures_factory/Fig1) | `Fig3_models/QC_metrics/make_counts_scatter_plots.ipynb` (source acknowledgement) |
| [LlorensLab/zamboni_et_al](https://github.com/LlorensLab/zamboni_et_al/blob/main/code/explore/ml_interpretations.R) | `Fig3_models/QC_metrics/interpretation.R` (source acknowledgement) |
| [kundajelab/chrombpnet](https://github.com/kundajelab/chrombpnet) | Model preparation, interpretation and prediction shell scripts |

External helper files and model weights referenced by these analyses have not been copied into this repository. Consult upstream attribution and licence files before redistributing or adapting their code. No new repository-wide licence is applied by this upload.
