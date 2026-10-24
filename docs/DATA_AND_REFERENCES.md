# Data and reference inputs

This code upload contains scripts, notebooks, small chromosome-split JSON files and environment specifications. It does not contain the large analysis inputs or model weights. The interactive [tree shrew atlas](https://sunkeyong.github.io/sun-lab-website/tree-shrew-atlas.html) is a separate resource; it is not a complete data package for every analysis in this repository.

| Analysis | Inputs to prepare |
| --- | --- |
| Read processing and atlas QC | FASTQ reads, barcode whitelists/correction tables, demultiplexing helpers, genome FASTA/GTF and alignment indices |
| Cell annotation | Fragment BED files with cell barcodes, cell metadata, reference annotations and the required ArchR/R objects |
| Peak calling and topic analysis | Per-cell-type fragment files, cell labels, peak BED/CSV files, pseudo-cell peak matrices/H5AD files, Mallet executable |
| Model training and interpretation | Fragment/peak BED files, genome FASTA and index, chromosome sizes and splits, trained bias model, model weights, motif resources and interpretation results |
| Cross-species analysis | Assembly-matched peak sets, forward/reverse liftOver chains, chromosome-name maps, gene/repeat/conservation annotations and species-specific cell labels |
| Figure plotting | Intermediate CSV/RData/H5AD/BigWig/HDF5 objects and any R objects referenced from previous interactive steps |

## Assembly and coordinate conventions

The tree shrew read-processing script explicitly references NCBI assembly `GCA_033439345.1`; other scripts use locally prepared `genome.fa`/`genes.gtf` paths and chromosome-name conversion tables. Human and mouse comparisons use hg38 and mm10 references in the supplied code. Match every FASTA, GTF, peak set, chain file and chromosome-split JSON to the exact assembly and naming convention used in that workflow. Do not mix accession-style chromosome names and renamed `chr` identifiers without the corresponding map.

Do not assume that a genome file prepared for the atlas website is interchangeable with every manuscript-analysis reference. Verify sequences, assembly version and chromosome names before reuse. Reference genomes and the custom R package `BSgenome.tupBel1.UCSC.TreeShrew` are not included in this upload.

## External sources referenced by the code

- [SPATAC-seq processing reference](https://github.com/Lan-lab/SPATAC-seq/)
- [NCBI tree shrew assembly GCA_033439345.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_033439345.1/)
- [TreeShrewDB](https://www.treeshrewdb.org/)
- [UCSC assembly resources](https://hgdownload.soe.ucsc.edu/hubs/GCA/033/439/345/GCA_033439345.1/)
- [UCSC hg38 liftOver resources](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/)
- [UCSC mm10 liftOver resources](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/)
- [Dfam 3.8 annotations](https://www.dfam.org/releases/Dfam_3.8/annotations/mm10/)

These links record references present in the supplied code. Availability, downloaded versions and the exact data accession list for the publication should be recorded when the study data release is finalised.
