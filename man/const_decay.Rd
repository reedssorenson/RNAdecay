% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decay_functions.R
\name{const_decay}
\alias{const_decay}
\alias{decaying_decay}
\title{exponential decay functions}
\usage{
const_decay(t, a)

decaying_decay(t, par)
}
\arguments{
\item{t}{time (in minutes)}

\item{a}{alpha (in per time, thus in per minute when time is in minutes)}

\item{par}{vector of length 2 containing alpha (par[1]) and beta (par[2]) values; alpha=initial decay rate, beta=decay of decay rate (both in per time, thus in per minute when time is in minutes)}
}
\value{
returns abundance after time \code{t}  at alpha initial decay rate and beta decay of decay rate relative to an initial abundance of 1
}
\description{
Constant decay rate function (const_decay(), case when betas=0) e^-a*t; decaying decay rate function (decaying_decay())  e^(-(a/b)*(1-e^(-b*t))). Functions are normalized so at t=0 the function is 1.
}
\examples{
const_decay(10,log(2)/10) ## returns 0.5
decaying_decay(10,c(log(2)/10,0.01)) ##returns 0.5170495
}
