# Tree shrew organogenesis atlas

Analysis code for the tree shrew organogenesis chromatin accessibility atlas, including single-cell processing, candidate cis-regulatory element (cCRE) analysis, sequence-to-function modelling and cross-species comparisons.

[Browse the atlas](https://sunkeyong.github.io/sun-lab-website/tree-shrew-atlas.html) · [Analysis guide](docs/ANALYSIS_GUIDE.md) · [Code index](docs/CODE_INDEX.md) · [Environment setup](environments/README.md)

## Repository layout

| Directory | Contents | Source collection |
| --- | --- | --- |
| [`sky_Code/`](sky_Code/) | Figures 1–5: read processing, QC, annotation, peak analysis and comparative analyses | Code curated by Keyong Sun |
| [`Fig2_pycistopic/`](Fig2_pycistopic/) | cisTopic construction, topic model analysis and visualisation | Collaborator-provided code |
| [`Fig3_models/`](Fig3_models/) | ChromBPNet training, interpretation, model QC and motif analyses | Collaborator-provided code |
| [`Fig5_immune_cells/`](Fig5_immune_cells/) | Human–tree shrew immune model prediction and motif comparisons | Collaborator-provided code |
| [`Fig6_enhancers/`](Fig6_enhancers/) | Enhancer candidate screening and sequence prediction examples | Collaborator-provided code |
| [`environments/`](environments/) | Five supplied Conda environment specifications | Supplied environment configurations |

Original directory names and script filenames have been retained so that figure references and contributor collections remain recognisable. The two source collections complement one another; the folder names do not establish individual authorship of every line. See [provenance and acknowledgements](docs/PROVENANCE.md).

## Getting started

```bash
git clone https://github.com/SunKeYong/TreeShrew-organogenesis-atlas.git
cd TreeShrew-organogenesis-atlas
```

1. Start with the [analysis guide](docs/ANALYSIS_GUIDE.md) to select a figure or workflow.
2. Prepare the inputs listed in [data and references](docs/DATA_AND_REFERENCES.md).
3. Select the appropriate environment using the [environment guide](environments/README.md).
4. Edit input/output paths, sample lists and resource settings in a working copy before running a script or notebook. Restart notebook kernels and execute the relevant sections in order.

This repository is a curated collection of research scripts and notebooks, rather than a single automated workflow. Many scripts retain the original analysis paths and interactive session assumptions. Raw sequencing data, reference genomes, matrices, trained models and large result files are not bundled here. The supplied Snakemake example covers a subset of model training, not the entire study.

## Validation and reproducibility

The initial upload checks file integrity, standalone Python and R syntax, shell syntax, notebook structure and Python-cell syntax, and YAML formatting. Notebook outputs and execution counts have been cleared for a compact source-code archive. Full analyses have not been rerun against the study data or GPUs as part of this upload. Environment exports retain original version/build pins and have not been re-solved.

See [reproducibility notes](docs/REPRODUCIBILITY.md) for required path edits, external helper code, example-only cells and the exact scope of validation.

## Attribution and citation

Please retain the original inline references and consult [PROVENANCE.md](docs/PROVENANCE.md) for third-party code acknowledgements. The manuscript citation and archival DOI can be added when available. A repository-wide software licence has not yet been selected; this upload does not relicense contributor or third-party code.
