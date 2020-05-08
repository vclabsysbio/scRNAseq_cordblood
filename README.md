scRNAseq_cordblood
=======
Script and samples for scRNA-seq cordblood analysis

## Under construction
This is a collection of script used for cordblood dataset analysis. There will be mainly Rscripts and some of shell script.

### Dependencies
This project requires the following R packages:
- [Suerat (3.1.4)](https://satijalab.org/seurat/)
- [doubletFinder (2.0.2)](https://github.com/chris-mcginnis-ucsf/DoubletFinder/blob/master/README.md) 
- [Monocle3 (0.2.1)](https://cole-trapnell-lab.github.io/monocle3/)
- [tidyverse (1.3.0)](https://www.tidyverse.org/)
- [velocyto.R (0.6)](https://github.com/velocyto-team/velocyto.R)
- [biomaRt (2.42.0)](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [ggplot2 (3.3.0)](https://ggplot2.tidyverse.org/)
- [gridExtra (2.3)](https://cran.r-project.org/web/packages/gridExtra/index.html)
- grid
- [sctransform (0.2.1)](https://github.com/ChristophH/sctransform)
- [future (1.16.0)](https://github.com/HenrikBengtsson/future)
- [fastcluster (1.1.25)](http://danifold.net/fastcluster.html)
- [ComplexHeatmap (2.2.0)](https://github.com/jokergoo/ComplexHeatmap)
- [circlize (0.4.8)](https://jokergoo.github.io/circlize_book/book/)
- [magritter (1.5)](https://cran.r-project.org/web/packages/magrittr/index.html)
- [uwot (0.1.5)](https://github.com/jlmelville/uwot)

### Project outlines
The bioinformatic analyses of this work can be divided into 3 main steps:
1. data preprocessing: import 10x files (mapped by CellRanger 3.0), preprocess, and remove doublet
2. data integration
3. downstream analysis: find diferentially expressed genes, visualise gene expression, draw trajectory
