% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommPath-main.R
\name{diffAllPath}
\alias{diffAllPath}
\title{To find differentially activated pathways in each identity class}
\usage{
diffAllPath(object, method = "t.test", only.posi = FALSE, only.sig = FALSE)
}
\arguments{
\item{object}{CommPath object}

\item{method}{Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'}

\item{only.posi}{only save the information of pathways showing up regulation}

\item{only.sig}{only save the information of pathways showing significant differences}
}
\value{
Data frame including the statistic result comparing the pathway enrichment sorces between cells in each cluster and all other clusters, the significant recetor and ligand in the pathways, and the corresponding up stream identity class and ligand
}
\description{
To find differentially activated pathways in each identity class
}
