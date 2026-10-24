# Code provenance and acknowledgements

## Contributors

- **Keyong Sun**: supplied the original `sky_Code` collection, now distributed across Figures 1–5 according to the figure prefixes in the filenames.
- **Ou Lin (欧霖)**: supplied the pycisTopic, model, immune-cell and enhancer code collections, now incorporated into Figures 2, 3, 5 and 6, together with the shared environment specifications.

Both contributors' code appears together within the relevant figure. The [code index](CODE_INDEX.md) and the per-figure README files identify the source contributor for each file. [FILE_MAP.tsv](FILE_MAP.tsv) records all original and current paths. Contributor labels record the supplied collections; original upstream acknowledgements are retained and should not be interpreted as sole authorship of every line or as a publication author order.

The reorganisation moves files without changing analysis code or configuration contents. The original local source collection is untouched, and the previous Git history retains the earlier layout.

## Third-party references retained in the supplied code

| Reference | Files that acknowledge or refer to it |
| --- | --- |
| [Lan-lab/SPATAC-seq](https://github.com/Lan-lab/SPATAC-seq/) | `Figure_1/Fig1_Process_fastq_files.sh` |
| [GreenleafLab/HDMA ChromBPNet use cases](https://github.com/GreenleafLab/HDMA/blob/main/code/05-misc/04-ChromBPNet_use_cases.ipynb) | Prediction notebooks under `Figure_3/plot_predicted/`, `Figure_5/cross_species_model_prediction/` and `Figure_6/plot_predicted_enhancer/` |
| [kundajelab/chrombpnet-figures](https://github.com/kundajelab/chrombpnet-figures/tree/main/src/figures_factory/Fig1) | `Figure_3/QC_metrics/make_counts_scatter_plots.ipynb` (source acknowledgement) |
| [LlorensLab/zamboni_et_al](https://github.com/LlorensLab/zamboni_et_al/blob/main/code/explore/ml_interpretations.R) | `Figure_3/QC_metrics/interpretation.R` (source acknowledgement) |
| [kundajelab/chrombpnet](https://github.com/kundajelab/chrombpnet) | Model preparation, interpretation and prediction shell scripts |

External helper files and model weights referenced by these analyses have not been copied into this repository. Consult upstream attribution and licence files before redistributing or adapting their code. No new repository-wide licence is applied by this upload.
