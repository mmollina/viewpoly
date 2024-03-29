% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_qtl.R
\name{select_haplo}
\alias{select_haplo}
\title{Plot selected haplotypes}
\usage{
select_haplo(
  input.haplo,
  probs,
  selected_mks,
  effects.data,
  exclude.haplo = NULL
)
}
\arguments{
\item{input.haplo}{character vector with selected haplotypes. It contains the information: "Trait:<trait ID>_LG:<linkage group ID_Pos:<QTL position>"}

\item{probs}{data.frame with first column (named `ind`) as individuals ID and next columns named with markers ID and containing the genotype probability at each marker}

\item{selected_mks}{data.frame with: LG - linkage group ID; mk - marker ID; pos - position in linkage map (cM)}

\item{effects.data}{output object from \code{data_effects} function}

\item{exclude.haplo}{character vector with haplotypes to be excluded. It contains the information: "Trait:<trait ID>_LG:<linkage group ID_Pos:<QTL position>"}
}
\value{
ggplot graphic
}
\description{
Plot selected haplotypes
}
\keyword{internal}
