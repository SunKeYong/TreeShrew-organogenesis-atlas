# Reproducibility notes

## Scope of the initial upload

The code has been organised and statically checked. No study analysis, model training, data download, dependency installation or command embedded in a notebook was executed during preparation. Static parsing does not establish numerical reproducibility or guarantee that a workflow will run with newly installed software.

## Required edits before execution

- Replace original `/home/...`, `/miniconda3/...`, `~/...`, `path/...`, sample-list placeholders and output directories with paths on your machine or cluster. `.libPaths()` calls must point to an appropriate local R installation.
- In Python, `os.chdir("~/...")` and `sys.path.append("~/...")` do not expand a home directory automatically. Replace these strings with absolute paths or use `os.path.expanduser(...)` in your working copy. Quoted `~` paths in shell variables and paths passed through configuration files also need explicit handling.
- In `Fig1_SubLibrary_QC.ipynb`, the illustrative `sample[10]` cell requires a sample list with at least 11 entries; it does not work with the single placeholder sample. Update it or skip that inspection cell after preparing your real sample list.
- Some R scripts are consecutive analysis/plotting sections that rely on objects loaded or created earlier in an interactive session. Restart R and prepare those objects deliberately; a syntax check alone cannot resolve these dependencies.
- `ss2_py.py`, `BC_process.py`, barcode correction tables and the custom tree shrew BSgenome package are external inputs, not files supplied by this code collection.
- Several prediction notebooks and the footprinting script import `tangermeme_utils` from an HDMA checkout. Obtain the appropriate helper code separately and update `sys.path`; it is not included as a standard pip package in this repository.
- Web uploads store shell scripts without executable mode. Invoke standalone scripts with `bash`; before using the Snakemake wrapper, run `chmod +x setup_environment.sh run.sh` in `Figure_3/model_training_and_analysis/Train_use_snakefile_example/`.
- The ChromBPNet Snakefile references `./setup_environment.sh` relative to its configured working directory. Make the wrapper available there or update the rule paths. The example uses GPU resources and a 64-GB rule memory request; adapt resource settings to your cluster. Its training rule references precomputed nonpeak BED files; prepare those explicitly before training. It is not a complete end-to-end study workflow.
- Random sampling, model seeds, versions, hardware and input files influence results. Some supplied snippets do not pin every random seed; preserve the actual run metadata when reproducing a figure.

## Source cleanup performed

1. Retained all 79 supplied non-temporary code/configuration files and original filenames. Code from Keyong Sun and Ou Lin (欧霖) is combined under `Figure_1/` through `Figure_6/`; shared configurations remain in `environments/`. The path mapping is recorded in `FILE_MAP.tsv`.
2. Removed notebook execution outputs/counts and transient execution/widget metadata; retained analysis source and order.
3. Converted two cells containing only the prose `Tree shrew sequences` into Markdown headings in the immune ortholog-prediction notebook.
4. Commented pasted console-table lines in `Fig4_Figure_plot1.R` and `Fig4_Figure_plot2.R` so R can parse them. Values remain visible as comments.
5. Removed the shell heredoc wrapper from `environments/chrombpnet/chrombpnetv8.yaml` and machine-specific `prefix:` entries from four Conda exports. Dependency/version/build entries are retained.
6. Excluded `.DS_Store` and empty `.Rhistory` files. Added documentation and ignore rules for future generated outputs.

No scientific thresholds, genomic intervals, model parameters or analytical code paths were intentionally changed. Original local source files were left untouched. The later figure reorganisation preserves every prepared source/configuration file byte for byte.

## Validation performed

- SHA-256 inventory of source and prepared files, with only the documented cleanup differences allowed.
- Syntax parsing for standalone Python and R scripts, and `bash -n` for shell scripts.
- JSON/notebook structure checks, cleared-output checks and syntax parsing of Python notebook cells. IPython timing/line magics and shell escape lines are excluded from plain-Python parsing; this is not a notebook execution test.
- YAML parsing for the five Conda specifications and the model-training configuration.
- Relative documentation links, duplicate cell IDs (where present), file inclusion and a targeted scan for credential-like strings.

The Snakefile has not been validated by a Snakemake dry run, and environments have not been installed or solved. Linux/GPU exports may need platform-specific changes. Additional environment specifications for the initial SnapATAC2/R analyses were not supplied; observed imports are listed in `environments/observed_imports.tsv` as an inventory, not a lockfile.
