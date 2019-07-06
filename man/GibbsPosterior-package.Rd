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
##\references{
##  This optional section can contain literature or other references for
##  background information.
##}
##\keyword{ package }
##\seealso{
##  Optional links to other man pages
##}
##\examples{
##  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
##  }
}
