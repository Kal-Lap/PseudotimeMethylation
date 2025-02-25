% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectSite.R
\name{select_site_pcc}
\alias{select_site_pcc}
\title{Select site by Pearson's Correlation Coefficient}
\usage{
select_site_pcc(meth_data, meta_data, pcc_thresh = 0.5)
}
\arguments{
\item{meth_data}{a data frame of methylation data}

\item{meta_data}{a data frame of meta data}

\item{pcc_thresh}{a numeric value of PCC threshold}
}
\value{
a vector of index that satisfy the PCC threshold
}
\description{
Select site with larger Pearson's Correlation Coefficient (PCC)
with chronological age than the given threshold
}
\examples{
\dontrun{select_site_pcc(meth_data, meta_data)
select_site_pcc(meth_data, meta_data, pcc_thresh = 0.002)}
}
