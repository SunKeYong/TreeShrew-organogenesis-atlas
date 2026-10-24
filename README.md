# Tree shrew organogenesis atlas

Analysis code for the tree shrew organogenesis chromatin accessibility atlas, including single-cell processing, candidate cis-regulatory element (cCRE) analysis, sequence-to-function modelling and cross-species comparisons.

**Code contributors: Keyong Sun and Lin Ou (欧霖).** Code from both contributors is organised together in manuscript figure order.

[Browse the atlas](https://sunkeyong.github.io/sun-lab-website/tree-shrew-atlas.html) · [Analysis guide](docs/ANALYSIS_GUIDE.md) · [Code index and contributors](docs/CODE_INDEX.md) · [Environment setup](environments/README.md)

## Code organised by figure

| Figure | Directory | Analyses | Code contributors |
| --- | --- | --- | --- |
| 1 | [`Figure_1/`](Figure_1/) | Atlas construction, QC, embedding and cell-type annotation | Keyong Sun |
| 2 | [`Figure_2/`](Figure_2/) | cCRE catalogue, peak embeddings, pycisTopic, motifs and annotations | Keyong Sun; Lin Ou |
| 3 | [`Figure_3/`](Figure_3/) | Sequence-to-function model preparation, training, QC and interpretation | Keyong Sun; Lin Ou |
| 4 | [`Figure_4/`](Figure_4/) | Cross-species cCRE conservation, divergence and repeat analyses | Keyong Sun |
| 5 | [`Figure_5/`](Figure_5/) | Immune-cell comparative analysis and cross-species predictions | Keyong Sun; Lin Ou |
| 6 | [`Figure_6/`](Figure_6/) | Enhancer candidate screening and sequence predictions | Lin Ou |

Shared environment specifications are in [`environments/`](environments/); analysis, data and provenance documentation is in [`docs/`](docs/). Within each figure, subdirectories describe analysis stages rather than contributor names. Each figure README and the [complete index](docs/CODE_INDEX.md) identify the source contributor for each file. [FILE_MAP.tsv](docs/FILE_MAP.tsv) records the original-to-current path mapping for all 79 supplied source/configuration files.

## Getting started

```bash
git clone https://github.com/SunKeYong/TreeShrew-organogenesis-atlas.git
cd TreeShrew-organogenesis-atlas
```

1. Select a figure in order using the [analysis guide](docs/ANALYSIS_GUIDE.md).
2. Prepare the [data and reference inputs](docs/DATA_AND_REFERENCES.md).
3. Choose the appropriate environment using the [environment guide](environments/README.md).
4. Edit input/output paths, sample lists and resource settings in a working copy before execution. Restart notebook kernels and run the relevant sections in order.

This is a curated collection of research scripts and notebooks. Many scripts retain original analysis paths and interactive-session assumptions. Raw sequencing data, reference genomes, matrices, trained models and large results are not bundled here. The Snakemake example covers a subset of model training, rather than the entire study.

## Validation and reproducibility

The initial source upload was statically checked for Python/R/shell syntax, notebook structure and Python-cell syntax, YAML formatting and file integrity. Notebook outputs and execution counts were cleared. The figure reorganisation preserves the bytes of all 79 source/configuration files; only their locations, navigation and contributor documentation have changed.

Full analyses have not been rerun against study data or GPUs, and environment exports have not been re-solved. See [reproducibility notes](docs/REPRODUCIBILITY.md) for required local edits and external helper code.

## Attribution and citation

See [code provenance and acknowledgements](docs/PROVENANCE.md) for contributions from Keyong Sun and Lin Ou and the third-party references retained in the supplied code. The manuscript citation and archival DOI can be added when available. A repository-wide software licence has not yet been selected; this upload does not relicense contributor or third-party code.
