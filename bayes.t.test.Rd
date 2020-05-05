\name{bayes.t.test}
\alias{bayes.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bayesttest
%%  ~~function to do ... ~~
}
\description{
The function provides a Bayesian t-test by running a Markov-Chain-Monte Carlo algorithm to produce posterior draws of the means and standard deviations (or variances) of both groups. The posterior distributions of the difference of means, difference of variances or standard deviations and the effect size are also obtained via the underlying Gibbs sampler for the two-component Gaussian mixture model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bayes.t.test(n, plot, firstComp, secondComp, hyperpars, ci, burnin, sd, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Number of iterations of the single-block Gibbs sampler used for inference in the two-component Gaussian mixture
}
  \item{plot}{
  A string specifying which plots to include for detailed analysis. "none" is used if only the data frame should be returned. "all" is used for detailed plots for an analysis based on the region of practical equivalence (ROPE).
}
  \item{firstComp}{
  A numeric vector representing the values of the first group
}
  \item{secondComp}{
  A numeric vector representing the values of the second group
}
  \item{hyperpars}{
  A string specifying the set of hyperparameters used for the prior distributions. "raftery" uses the standard hyperparameters of Raftery (1996), and "custom" uses modified parameters which are designed to stay weakly informative.
}
  \item{ci}{
  A real number between zero and one representing the credible level used for posterior credible intervals
}
  \item{burnin}{
  An integer which represents the burn-in used for the single-block Gibbs sampler
}
  \item{sd}{
  A string selecting if posterior inference is desired for the variances or standard deviations. "var" approximates the posterior distribution of the variances, "sd" of the standard deviations.
}
  \item{q}{
  A real number for influencing the hyperparameters selected. q=0.1 as weakly informative standard is recommended.
}
}
\value{
  If the option "none" is used for argument plot, a dataframe including the posterior mode, posterior expectation as well as the posterior credible intervals corresponding to the selected credible level will be returned for the difference of means, difference of standard deviations respectively variances and the effect size. If the option "all" is used for argument plot, plots of the posterior distributions of the difference of means, difference of standard deviations respectively variances and the effect size will be returned, including a posterior analysis of the effect size via the region of practical equivalence (ROPE).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
For details, see: https://arxiv.org/abs/1906.07524v1
}
\author{
Riko Kelter
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
x=rnorm(100,mean=0,sd=1)
y=rnorm(100,mean=1,sd=1)
bayes.t.test(n=2500,plot="all",firstComp=x,secondComp=y,
hyperpars="custom",ci=0.95,burnin=1250,sd="sd",q=0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.