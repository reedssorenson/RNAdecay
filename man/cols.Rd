% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cols.R
\name{cols}
\alias{cols}
\title{Indexes column names of a data.frame matching multiple patterns (i.e., multigrep)}
\usage{
cols(patterns, df, w = NA, x = NA, y = NA, z = NA)
}
\arguments{
\item{patterns}{character vector or vector of regular expressions passed to grep pattern argument}

\item{df}{a dataframe with column names to index}

\item{w, x, y, z}{(for backwards compatibility) separate arguments for patterns, if used patterns argument will be ignored}
}
\value{
returns a vector of integer indices of the column names of \code{df} that match to all of \code{patterns}
}
\description{
Identifies dataframe column names that have all of the pattern arguments .
}
\details{
Be aware that column data labels that are part of another data label are not advisable (e.g. mut1, mut2, mut1.mut2; cols(df,'mut1') will return indices for both 'mut1' and 'mut1.mut2' labeled columns
}
\examples{
cols(df=data.frame(xyz=1:5,zay=6:10,ybz=11:15,tuv=16:20),patterns = c('y','z')) ## returns 1 2 3
cols(df=data.frame(xyz=1:5,zay=6:10,ybz=11:15,tuv=16:20), w = 'y', x = 'z') ## returns 1 2 3

}
