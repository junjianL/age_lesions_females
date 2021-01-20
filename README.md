### Disentangling tumor- from age-associated DNA methylation changes occurring in the colonic mucosa of women

## Code and data analysis

This repository contains the code for the following paper:

* S Orjuela, HR Parker, S Sajibu, F Cereatti, M Sauter, F Buffoli, MD Robinson and G Marra: [Disentangling tumor- from age-associated DNA methylation changes occurring in the colonic mucosa of women]().

All files outside the `scripts/` folder are necessary to run a Snakemake workflow with:

`snakemake --use-conda -cores 1` similar to the [ARMOR pipeline](https://github.com/csoneson/ARMOR)

Before running the pipeline, raw FASTQ files must be downloaded from [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB27521) (lesion females) and [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB40783) (healthy females). The location of the FASTQ files must be specified in `config.yaml`. Please download the files from the "Submitted FTP" section, because the file names match the names in `metadata_*.txt`

*Folder description*: 
* `scripts/`
R scripts with code for DMR detection analysis, and code for plotting most of figures in paper.
* `data/` .RData objects with DMR lists mentioned in the publication.
