# Environment setup

Environment specifications contributed by **Lin Ou**. These shared environments support the figure-specific analyses in `Figure_1/` through `Figure_6/`.

Use separate environments for the relevant analysis stages. These are supplied research environment specifications; they have not been freshly installed or solved during this upload. They retain original package versions and build identifiers, including Linux/GPU-specific dependencies. They are not guaranteed to solve on macOS or on every CUDA installation.

| Environment file | Supplied name | Main use |
| --- | --- | --- |
| [`pycistopic/pycisTopic.yaml`](pycistopic/pycisTopic.yaml) | `pycisTopic` | cisTopic/Mallet analysis |
| [`chrombpnet/chrombpnetv8.yaml`](chrombpnet/chrombpnetv8.yaml) | `chrombpnetv8` | ChromBPNet training and interpretation |
| [`modisco/modisco-lite.yaml`](modisco/modisco-lite.yaml) | `modisco-lite` | TF-MoDISco motif discovery |
| [`finemo/Finemo.yaml`](finemo/Finemo.yaml) | `Finemo` | Fi-NeMo motif occurrence analysis |
| [`tangermeme/tangermeme.yaml`](tangermeme/tangermeme.yaml) | `tangermeme` | Model predictions and sequence/motif analyses |

For example, from the repository root:

```bash
conda env create -f environments/chrombpnet/chrombpnetv8.yaml
conda activate chrombpnetv8
```

Inspect the environment first and match CUDA/cuDNN/TensorFlow/PyTorch requirements to your actual GPU host. The supplied ChromBPNet YAML does not by itself install all external tool source checkouts used in the scripts. Obtain the corresponding ChromBPNet implementation, Snakemake, Mallet and any external helpers required for the chosen workflow; retain their actual versions in run records.

The original machine-specific `prefix:` fields have been removed, and the ChromBPNet YAML's surrounding shell heredoc has been removed. All listed dependency/version/build entries are otherwise retained. Some exports may reference distribution-specific packages or editable source locations; review those before installation.

The original code also requires SnapATAC2/Scanpy and multiple R/Bioconductor packages, for which complete environment exports were not supplied. [`observed_imports.tsv`](observed_imports.tsv) inventories literal imports detected in the provided code. It is not a requirements file and does not include every dynamic import, command-line executable, custom package, data package or transitive dependency.

See [reproducibility notes](../docs/REPRODUCIBILITY.md) and [data requirements](../docs/DATA_AND_REFERENCES.md).
