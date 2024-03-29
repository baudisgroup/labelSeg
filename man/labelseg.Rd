% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labelseg.R
\name{labelseg}
\alias{labelseg}
\title{Annotate segment profiles}
\usage{
labelseg(
  data,
  genome = "hg38",
  baseshift = "n",
  method = "dbscan",
  len.sep = TRUE,
  minPts = 1,
  lb.dup = 0.15,
  hb.dup = 0.7,
  lb.del = 0.15,
  hb.del = 1.5,
  shiftnum = 1,
  labeled = FALSE
)
}
\arguments{
\item{data}{CNA segment data with logR values}

\item{genome}{Genome version of input data. Available options are "hg38","hg19" and "hg18". Default is "hg38".}

\item{baseshift}{Direction of baseline shift. Available options are  "n" (no shift), "h" (shift higher), and "l" (shift lower). Default is "n".}

\item{method}{Clustering method used in this algorithm. Available options are "dbscan", "hdbscan", and "optics". Default is "dbscan".}

\item{len.sep}{Logic value to determine whether separate clustering by segment length is applied. Default is TRUE.}

\item{minPts}{Integer to specify the minimum size of clusters. Default is 1.}

\item{lb.dup}{Lower bound for calling low-level duplication. Default is 0.15.}

\item{hb.dup}{Higher bound for calling low-level duplication. Default is 0.7.}

\item{lb.del}{Lower bound for calling low-level deletion. Default is 0.15.}

\item{hb.del}{Higher bound for calling low-level deletion. Default is 1.5.}

\item{shiftnum}{Number of target baseline clusters to shift. Only used if the parameter \code{baseshift} is "h" or "l".}

\item{labeled}{Logic value to determine whether the input data includes labels. If set to TRUE, the baseline will be determined by the original callings. Default is FALSE.}
}
\value{
Segment data with labels. In terms of SCNA labels, "-2" represents high-level focal deletion. "-1" represents
low-level deletion. "+2" represents high-level focal duplication. "+1" represents low-level duplication.
"0" represents no copy changes.
}
\description{
This function assigns segments SCNA labels to call different levels of SCNA.
It estimates calling thresholds based on logR and segment length of individual samples.
}
\examples{
data(example.seg)
labeled.seg <- labelseg(data=example.seg)
}
