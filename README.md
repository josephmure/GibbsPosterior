# Gibbs reference posterior distribution

This R implementation consists in 2 packages that need to be installed in this order:

1. `GibbsPosteriorRequired`
2. `GibbsPosterior`

As `GibbsPosterior` uses C++ code, the package [Rcpp](https://cran.r-project.org/web/packages/Rcpp/) is a prerequisite.
Moreover, `GibbsPosterior` makes use of the [GNU Scientific Library](https://www.gnu.org/software/gsl/) so the package [RcppGSL](https://cran.r-project.org/web/packages/RcppGSL/) is also required.

## Installation

If it is not already done, install [Rcpp](https://cran.r-project.org/web/packages/Rcpp/) and [RcppGSL](https://cran.r-project.org/web/packages/RcppGSL/).

Move to the `GibbsPosterior/` folder and execute the following commands:

```
R CMD build GibbsPosteriorRequired
R CMD INSTALL GibbsPosteriorRequired
cd ..
R CMD build GibbsPosterior
R CMD INSTALL GibbsPosterior
```

## Usage

The `GibbsPosteriorC` functions samples from the Gibbs reference posterior distribution on the correlation lengths for Kriging models with Matérn anisotropic (tensorized or geometric) kernel.
It does this using a Metropolis-within-Gibbs algorithm, as detailed in the article [**Optimal compromise between incompatible conditional probability distributions, with application to Objective Bayesian Kriging**](https://www.esaim-ps.org/articles/ps/abs/2019/02/ps170094/ps170094.html)

It expects the following arguments, in this order:

1. The dimension of the input space.
2. The smoothness of the Matérn kernel.
3. The number of points that must be sampled from the Gibbs reference posterior distribution.
4. Thinning factor (the length of the Gibbs sampling Markov chain is number of points to be generated * thinning factor).
5. Standard deviation used in the Metropolis part of the algorithm.
6. Dimension of the trend vector space.
7. Seed of the random generator.
8. Vector of two strings. The first string must be "REML": it corresponds to an option that is now deprecated. The second string must be either "tensorise" if the anisotropic Matérn kernel should be tensorized or "geometrique" if it should be anisotropic geometric.
9. Matrix of inputs. Every row corresponds to a point, every column to a dimension of the input space.
10. Vector of outputs. Its size shoulf be the same as the number of rows of the matrix of inputs.
11. Trend matrix. Each columns contains the values of one basis function of the trend vector space at the rows in the matrix of inputs.
12. Starting point of the Gibbs algorithm (i.e. the initial values of the correlation lengths).

`GibbsPosteriorC' produces an R object whose "posterior" attribute contains the sample from the Gibbs reference posterior distribution.

An exemple is given in the script `test_basique.r` in the `tests/` folder.
