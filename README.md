# LabelSeg

This is an R package to call different levels of SCNA from segment profiles.

You can install this package from GitHub using:

```r
install.packages("devtools")
devtools::install_github("baudisgroup/LabelSeg")
```

## Example use

```{r}
labelseg(data=segment_data, genome="hg38",baseshift = 'n')
```
