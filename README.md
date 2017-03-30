# stageR

This is the repository for the stageR package. stageR allows user-friendly automated stage-wise analysis of high-throughput genomic data.

To install the package from the GitHub repository in R please use

```
library(devtools)
if(!all(c("BiocStyle","BiocInstaller") %in% installed.packages()[,1])){
source("https://bioconductor.org/biocLite.R")
biocLite("BiocStyle")
}
install_github("statOmics/stageR")
```

The repository containing all code required to reproduce the analyses in the paper can be found at http://www.github.com/statOmics/stageWiseTestingPaper.

