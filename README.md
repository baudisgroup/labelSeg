# labelSeg

This is an R package designed to identify and label different levels of Copy Number Alterations (CNA) in segmented profiles.

## Installation

To install the package, you can use the `devtools` package as follows:

```r
install.packages("devtools")
devtools::install_github("baudisgroup/labelSeg")
```

## Input

The package requires segmentation data in the form of a *.seg* file. This file should contain six columns to specify the following information:

1. Sample identifier
2. Chromosome
3. Segment start position
4. Segment end position
5. Number of markers in the segment (optional)
6. LogR values

## Output 

After processing your input data, labelSeg will add an additional column called "label" to your original segment data. This column indicates the SCNA levels for each segment. The labels are as follows:

* "+2" and "-2" represent high-level focal duplication and deletion.
* "+1" and "-1" represent low-level (broad) duplication and deletion.
* "0" represents segments with no relative copy changes.

## Example use

### Load library and example data

```{r}
library(labelSeg)
data(example.seg)
```
### Run 

```{r}
labeled.seg <- labelseg(data=example.seg, genome="hg38")
```
