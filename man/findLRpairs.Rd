% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommPath-main.R
\name{findLRpairs}
\alias{findLRpairs}
\title{To find marker ligands and marker receptors}
\usage{
findLRpairs(object, logFC.thre = 0, p.thre = 0.05)
}
\arguments{
\item{object}{CommPath object}

\item{logFC.thre}{logFC threshold, marker genes with a logFC > logFC.thre will be considered}

\item{p.thre}{p threshold, marker genes with a adjust p value < p.thre will be considered}
}
\value{
List containing the ligand-receptor interaction information
}
\description{
To find marker ligands and marker receptors
}
