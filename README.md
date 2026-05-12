# SCI bulk RNA-seq analysis

Code used for the analysis part of my Master's thesis project at Osaka University / Karolinska Institutet.

This repo is mainly just a cleaned version of the analysis scripts used for the RNA-seq and qPCR figures/results. I removed the unfinished notebooks, raw data, intermediate files, and large reference objects before making it public.

Most scripts are meant to be run manually from the repository root in roughly numerical order.

The repository includes:
- differential expression analysis
- PCA/QC
- pathway enrichment
- Reactome/community analyses
- program scoring
- qPCR validation analysis

Raw sequencing data and private metadata are not included.

Example usage:

```bash
Rscript scripts/01_limma_voom_de.R
python scripts/02_pca_sample_identity.py


