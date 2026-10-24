# Figure 3: Sequence-to-function models

Peak preparation, bias and ChromBPNet training, QC, interpretation, prediction plots, motif co-occurrence and footprinting.

Code contributed by Keyong Sun and Lin Ou.

| File | Type | Code contributor |
| --- | --- | --- |
| [Fig3_PeakCallling_For_ChromBPNet.ipynb](Fig3_PeakCallling_For_ChromBPNet.ipynb) | Jupyter notebook | Keyong Sun |
| [Motif_cooccurrence/motif_cooccurrence_allpeaks.ipynb](Motif_cooccurrence/motif_cooccurrence_allpeaks.ipynb) | Jupyter notebook | Lin Ou |
| [Motif_cooccurrence/plot_cooccurrence_allpeaks.R](Motif_cooccurrence/plot_cooccurrence_allpeaks.R) | R | Lin Ou |
| [QC_metrics/get_observed_counts.py](QC_metrics/get_observed_counts.py) | Python | Lin Ou |
| [QC_metrics/get_observed_counts_rep.py](QC_metrics/get_observed_counts_rep.py) | Python | Lin Ou |
| [QC_metrics/interpretation.R](QC_metrics/interpretation.R) | R | Lin Ou |
| [QC_metrics/make_counts_scatter_plots.ipynb](QC_metrics/make_counts_scatter_plots.ipynb) | Jupyter notebook | Lin Ou |
| [QC_metrics/random_split_bed.sh](QC_metrics/random_split_bed.sh) | Shell | Lin Ou |
| [model_training_and_analysis/MotifPlotting.ipynb](model_training_and_analysis/MotifPlotting.ipynb) | Jupyter notebook | Lin Ou |
| [model_training_and_analysis/Train_bias_model/make_splits.ipynb](model_training_and_analysis/Train_bias_model/make_splits.ipynb) | Jupyter notebook | Lin Ou |
| [model_training_and_analysis/Train_bias_model/preprocess_trainBias.sh](model_training_and_analysis/Train_bias_model/preprocess_trainBias.sh) | Shell | Lin Ou |
| [model_training_and_analysis/Train_bias_model/splits/fold_0.json](model_training_and_analysis/Train_bias_model/splits/fold_0.json) | Chromosome-split JSON | Lin Ou |
| [model_training_and_analysis/Train_bias_model/splits/fold_1.json](model_training_and_analysis/Train_bias_model/splits/fold_1.json) | Chromosome-split JSON | Lin Ou |
| [model_training_and_analysis/Train_bias_model/splits/fold_2.json](model_training_and_analysis/Train_bias_model/splits/fold_2.json) | Chromosome-split JSON | Lin Ou |
| [model_training_and_analysis/Train_bias_model/splits/fold_3.json](model_training_and_analysis/Train_bias_model/splits/fold_3.json) | Chromosome-split JSON | Lin Ou |
| [model_training_and_analysis/Train_bias_model/splits/fold_4.json](model_training_and_analysis/Train_bias_model/splits/fold_4.json) | Chromosome-split JSON | Lin Ou |
| [model_training_and_analysis/Train_use_snakefile_example/inputs.yaml](model_training_and_analysis/Train_use_snakefile_example/inputs.yaml) | YAML | Lin Ou |
| [model_training_and_analysis/Train_use_snakefile_example/run.sh](model_training_and_analysis/Train_use_snakefile_example/run.sh) | Shell | Lin Ou |
| [model_training_and_analysis/Train_use_snakefile_example/setup_environment.sh](model_training_and_analysis/Train_use_snakefile_example/setup_environment.sh) | Shell | Lin Ou |
| [model_training_and_analysis/Train_use_snakefile_example/snakefile](model_training_and_analysis/Train_use_snakefile_example/snakefile) | Snakemake | Lin Ou |
| [model_training_and_analysis/finemo.sh](model_training_and_analysis/finemo.sh) | Shell | Lin Ou |
| [model_training_and_analysis/modicso.sh](model_training_and_analysis/modicso.sh) | Shell | Lin Ou |
| [plot_predicted/predict_prebw_marker_region.sh](plot_predicted/predict_prebw_marker_region.sh) | Shell | Lin Ou |
| [plot_predicted/tangermeme.markerpeak.ipynb](plot_predicted/tangermeme.markerpeak.ipynb) | Jupyter notebook | Lin Ou |
| [solo_motif_footprinting/test_solo_motif_footprinting_all_celltypes_complete.py](solo_motif_footprinting/test_solo_motif_footprinting_all_celltypes_complete.py) | Python | Lin Ou |

See the [analysis guide](../docs/ANALYSIS_GUIDE.md), [environment setup](../environments/README.md), [data requirements](../docs/DATA_AND_REFERENCES.md) and [reproducibility notes](../docs/REPRODUCIBILITY.md) before execution. The table lists files for navigation, rather than a single automatic execution sequence; refer to the guide for the relevant workflow stages.
