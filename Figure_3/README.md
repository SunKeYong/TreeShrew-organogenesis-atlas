# Figure 3: Sequence-to-function models

Peak preparation, bias and ChromBPNet training, QC, interpretation, prediction plots, motif co-occurrence and footprinting.

Code contributed by Keyong Sun and Ou Lin (欧霖).

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig3_PeakCallling_For_ChromBPNet.ipynb](Fig3_PeakCallling_For_ChromBPNet.ipynb) | Jupyter notebook | Keyong Sun |
| [Motif_cooccurrence/motif_cooccurrence_allpeaks.ipynb](Motif_cooccurrence/motif_cooccurrence_allpeaks.ipynb) | Jupyter notebook | Ou Lin |
| [Motif_cooccurrence/plot_cooccurrence_allpeaks.R](Motif_cooccurrence/plot_cooccurrence_allpeaks.R) | R | Ou Lin |
| [QC_metrics/get_observed_counts.py](QC_metrics/get_observed_counts.py) | Python | Ou Lin |
| [QC_metrics/get_observed_counts_rep.py](QC_metrics/get_observed_counts_rep.py) | Python | Ou Lin |
| [QC_metrics/interpretation.R](QC_metrics/interpretation.R) | R | Ou Lin |
| [QC_metrics/make_counts_scatter_plots.ipynb](QC_metrics/make_counts_scatter_plots.ipynb) | Jupyter notebook | Ou Lin |
| [QC_metrics/random_split_bed.sh](QC_metrics/random_split_bed.sh) | Shell | Ou Lin |
| [model_training_and_analysis/MotifPlotting.ipynb](model_training_and_analysis/MotifPlotting.ipynb) | Jupyter notebook | Ou Lin |
| [model_training_and_analysis/Train_bias_model/make_splits.ipynb](model_training_and_analysis/Train_bias_model/make_splits.ipynb) | Jupyter notebook | Ou Lin |
| [model_training_and_analysis/Train_bias_model/preprocess_trainBias.sh](model_training_and_analysis/Train_bias_model/preprocess_trainBias.sh) | Shell | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_0.json](model_training_and_analysis/Train_bias_model/splits/fold_0.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_1.json](model_training_and_analysis/Train_bias_model/splits/fold_1.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_2.json](model_training_and_analysis/Train_bias_model/splits/fold_2.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_3.json](model_training_and_analysis/Train_bias_model/splits/fold_3.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_bias_model/splits/fold_4.json](model_training_and_analysis/Train_bias_model/splits/fold_4.json) | Chromosome-split JSON | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/inputs.yaml](model_training_and_analysis/Train_use_snakefile_example/inputs.yaml) | YAML | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/run.sh](model_training_and_analysis/Train_use_snakefile_example/run.sh) | Shell | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/setup_environment.sh](model_training_and_analysis/Train_use_snakefile_example/setup_environment.sh) | Shell | Ou Lin |
| [model_training_and_analysis/Train_use_snakefile_example/snakefile](model_training_and_analysis/Train_use_snakefile_example/snakefile) | Snakemake | Ou Lin |
| [model_training_and_analysis/finemo.sh](model_training_and_analysis/finemo.sh) | Shell | Ou Lin |
| [model_training_and_analysis/modicso.sh](model_training_and_analysis/modicso.sh) | Shell | Ou Lin |
| [plot_predicted/predict_prebw_marker_region.sh](plot_predicted/predict_prebw_marker_region.sh) | Shell | Ou Lin |
| [plot_predicted/tangermeme.markerpeak.ipynb](plot_predicted/tangermeme.markerpeak.ipynb) | Jupyter notebook | Ou Lin |
| [solo_motif_footprinting/test_solo_motif_footprinting_all_celltypes_complete.py](solo_motif_footprinting/test_solo_motif_footprinting_all_celltypes_complete.py) | Python | Ou Lin |

See the [analysis guide](../docs/ANALYSIS_GUIDE.md), [environment setup](../environments/README.md), [data requirements](../docs/DATA_AND_REFERENCES.md) and [reproducibility notes](../docs/REPRODUCIBILITY.md) before execution. The table lists files for navigation, rather than a single automatic execution sequence; refer to the guide for the relevant workflow stages.
