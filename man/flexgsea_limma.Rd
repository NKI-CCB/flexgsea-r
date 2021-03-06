% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexgsea_limma.R
\name{flexgsea_limma}
\alias{flexgsea_limma}
\title{Score genes with limma with a design matrix.}
\usage{
flexgsea_limma
}
\description{
Scores genes by their coefficients in from limma. 
The argument \code{y} should be a design matrix, as produced by
\code{model.matrix}. 
For RNA-seq data, use the
\code{EList} object produced by the \code{\link[limma]{voom}} function.
Do not call directly, but give as the \option{gene.score.fn}  argument to
\code{\link{flexgsea}}.
}
\seealso{
Other gene scoring functions: \code{\link{flexgsea_lm}},
  \code{\link{flexgsea_s2n}}
}
\concept{gene scoring functions}
