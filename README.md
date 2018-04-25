# stageR

This is the repository for the stageR package. stageR allows user-friendly automated stage-wise analysis of high-throughput genomic data.

To install the package from the GitHub repository in R please use

```
library(devtools)
install_github("statOmics/stageR")
```

Note, that we discovered a bug in the ‘holm’ and ‘user’ corrections for the ‘stageWiseAdjustment’ function that was introduced when making changes for Bioconductor submission. The bug was present in version 0.99.08 until version 1.1.1. It is also present in the Bioconductor 3.6 release branch, hence it is preferable to install the package through Github.



The repository containing all code required to reproduce the analyses in the paper can be found at http://www.github.com/statOmics/stageWiseTestingPaper.

