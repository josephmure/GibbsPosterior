\name{GibbsPosterior-package}
\alias{GibbsPosterior-package}
\alias{GibbsPosterior}
\docType{package}
\title{
  Samples from the Gibbs reference posterior distribution.
}
\description{
  In anisotropic Kriging models, estimating the correlation lengths
  is a difficult problem. The Gibbs reference posterior distribution
  is an Objective Bayesian solution. This package allows the user to
  sample from the Gibbs reference posterior distribution on correlation
  lengths.
}
\details{
  The GibbsPosteriorC functions samples from the Gibbs reference 
  posterior distribution on the correlation lengths for Kriging models
  with Matern anisotropic (tensorized or geometric) kernel.
  It does this using a Metropolis-within-Gibbs algorithm.
}
\author{
Joseph MURE <joseph.mure@edf.fr>
}
\references{
     J.Muré (2019), _Optimal compromise between incompatible 
     conditional probability distributions, with application 
     to Objective Bayesian Kriging_, _ESAIM: Probability and 
     Statistics_, *23*, 271-309.
}
\keyword{ Gibbs reference posterior }
}
