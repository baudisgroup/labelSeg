# LabelSeg

This is an R package to call different levels of SCNA from segment profiles.

You can install this package from GitHub using:

```r
install.packages("devtools")
devtools::install_github("baudisgroup/LabelSeg")
```

## Example use

```{r}
# load library and example data
library(LabelSeg)
segfile <- system.file("extdata", "example.seg", package = "LabelSeg")
data <- read.table(segfile,header=T,sep='\t')

# run LabelSeg
labelled.data <- labelseg(data=data, genome="hg38")
```
