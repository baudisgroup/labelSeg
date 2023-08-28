# labelSeg

This is an R package to call different levels of SCNA from segment profiles.

You can install this package from GitHub using:

```r
install.packages("devtools")
devtools::install_github("baudisgroup/labelSeg")
```

## Example use

```{r}
# load library and example data
library(labelSeg)
segfile <- system.file("extdata", "example.seg", package = "labelSeg")
segdata <- read.table(segfile,header=T,sep='\t')

# run labelSeg
labeled.data <- labelseg(data=segdata)
```
