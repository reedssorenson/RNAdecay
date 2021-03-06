% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_optimization.R
\name{mod_optimization}
\alias{mod_optimization}
\title{model optimization for fitting exponential decay models to normalized data}
\usage{
mod_optimization(
  gene,
  data,
  alpha_bounds,
  beta_bounds,
  models,
  group,
  mod,
  file_only = TRUE,
  path = "modeling_results"
)
}
\arguments{
\item{gene}{geneID from \code{data} to be modeled}

\item{data}{decay data data.frame with columns named 'geneID', 'treatment',
't.decay', 'rep', 'value.'}

\item{alpha_bounds}{vector of length 2 with lower and upper bounds for
alpha}

\item{beta_bounds}{vector of length 2 with lower and upper bounds for beta}

\item{models}{vector spceifying which models to run optimization on (e.g.,
c('mod1', 'mod239'))}

\item{group}{grouping matrix for alphas or betas}

\item{mod}{data.frame specifying alpha and beta group pairs for each model}

\item{file_only}{logical; should output only be written to file (TRUE) or
also return a data.frame of the results (FALSE)}

\item{path}{specify folder for output to be written}
}
\value{
returns (if \code{file_only = FALSE}) and writes to \code{path} a
  data frame of model optimization results for \code{models} one row for
  each for \code{gene} using values for it found in \code{data}, the columns
  of the data frame are:
geneID, mod (model), model estimates [alpha_treatment1, ...,
  alpha_treatmentn, beta_treatment1, ..., beta_treatmentn, sigma2],
 logLik (maximum log likelihood), nPar (number of parameters in the model),
 nStarts (number of parameter starting value sets (of 50) that converged on
 a maximum likelihood peak), J (number of parameter starting value sets that
 converged on the highest - within 1e-4 - maximum likelihood of all
 parameter starting value sets), range.LL (range of maximum likelihoods
 values reached by algorithm convergence from all parameter starting value
 sets), nUnique.LL (number of unique maximum likelihoods values reached by
 algorithm convergence from all parameter starting value sets), C.alpha (sum
 of all coefficients of variation for each column of alpha estimates),
 C.beta (sum of all coefficients of variation for each column of beta
 estimates), C.tot (C.alpha+C.beta), AICc (calculated from the single
 highest maximum likelihood of all parameter starting value sets), AICc_est
 (calculated from the log likelihood value computed by using the mean of
 each parameter from all optimizations that converged on the highest
 maximum likelihood of all starting parameter value sets.)
}
\description{
The mod_optimization function finds the estimates of model parameters by
  maximum likelihood, for a single gene on a specified list of models, and
  saves a tab delimited text file of the results named,'
  [geneID]_results.txt'.
The function does the following for each gene:
(1) it calculates log likelihood for each point in a 2 dimensional grid of
  evenly spaced alpha and beta values within the alpha and beta bounds
  specified using the null model (in which all treatment alphas are
  equivalent and all betas are equivalent).
(2) it calculates log likelihood for each point in a 1 dimensional range of
  evenly spaced alpha values within the alpha bounds using the single
  exponential null model (in which all treatment alphas are equivalent).
(3) For each of the grid points with the highest log likelihood from steps
  (1) and (2) 25 starting parameter value sets that are normally
  distributed around these points are generated.
(4) Parameter values are optimized for maximum likelihood using each of
  these 50 starting parameter sets using pre-compiled C++ functions loaded
 from dynamically linked libraries stored in the package on all models
 specified in the models argument.
(5) evaluates parameter estimates of all 50 optimizations based on the
  reported maximum liklihood upon convergence. Only parameter estimates
  that converged on the same and highest maximum likelihood are returned.
(6) returns the optimized parameter estimates, with model selection
  statistics.
}
\examples{

mod_optimization(gene = 'Gene_BooFu',
                data = data.frame(geneID=rep('Gene_BooFu',30),
                            treatment=c(rep('WT',15),rep('mut',15)),
                            t.decay=rep(c(0,7.5,15,30,60),6),
                            rep=rep(paste0('rep',c(rep(1,5),rep(2,5),rep(3,5))),2),
                            value= c(0.9173587, 0.4798672, 0.3327807, 0.1990708, 0.1656554,
                                     0.9407511, 0.7062988, 0.3450886, 0.3176824, 0.2749946,
                                     1.1026497, 0.6156978, 0.4563346, 0.2865779, 0.1680075,
                                     0.8679866, 0.6798788, 0.2683555, 0.5120951, 0.2593122,
                                     1.1348219, 0.8535835, 0.6423996, 0.5308946, 0.4592902,
                                     1.1104068, 0.5966838, 0.3949790, 0.3742632, 0.2613560)),
                alpha_bounds = c(1e-4,0.75),
                beta_bounds = c(1e-3,0.075),
                models = 'mod1',
                group = t(matrix(c(1,2,1,1,NA,NA),nrow=2,
                          dimnames=list(c('treat1','treat2'),c('mod1','mod2','mod3')))),
                mod = as.data.frame(t(matrix(c(1,1,1,2,1,3,2,1,2,2,2,3),nrow=2,
                        dimnames=list(c('a','b'),paste0('mod',1:6))))),
                file_only = FALSE,
                path = paste0(tempdir(),"/modeling results"))

}
