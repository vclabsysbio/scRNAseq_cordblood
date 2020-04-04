scRNAseq_cordblood
=======
Script and samples for scRNA-seq cordblood analysis

## Under construction
This is a collection of script used for cordblood dataset analysis. There will be mainly Rscripts and some of shell script.

### Dependencies
This project requires the following R packages:
- [Suerat (3.0)](https://cole-trapnell-lab.github.io/monocle3/)
- [doubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder/blob/master/README.md) 
- [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/)

### Project outlines
The bioinformatic analyses of this work can be divided into 3 main steps:
1. data preprocessing: import 10x files (mapped by CellRanger 3.0), preprocess, and remove doublet
2. data integration
3. downstream analysis: find diferentially expressed genes, visualise gene expression, draw trajectory
