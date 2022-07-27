\name{goodTuring}
\alias{goodTuring}
\alias{goodTuringPlot}
\alias{goodTuringProportions}

\title{Good-Turing Frequency Estimation}

\description{
Non-parametric empirical Bayes estimates of the frequencies of observed (and unobserved) species.
}

\usage{goodTuring(x, conf=1.96)
goodTuringPlot(x)
goodTuringProportions(counts)}

\arguments{ 
	\item{x}{numeric vector of non-negative integers, representing the observed frequency of each species.}
	\item{conf}{confidence factor, as a quantile of the standard normal distribution, used to decide for what values the log-linear relationship between frequencies and frequencies of frequencies is acceptable.}
	\item{counts}{matrix of counts}
}

\value{
\code{goodTuring} returns a list with components
	\item{count}{observed frequencies, i.e., the unique positive values of \code{x}}
	\item{n}{frequencies of frequencies}
	\item{n0}{frequency of zero, i.e., number of zeros found in \code{x}}
	\item{proportion}{estimated proportion of each species given its count}
	\item{P0}{estimated combined proportion of all undetected species}
\code{goodTuringProportions} returns a matrix of proportions of the same size
as \code{counts}.
}

\details{
Observed counts are assumed to be Poisson distributed.
Using an non-parametric empirical Bayes strategy, the algorithm evaluates the posterior expectation of each species mean given its observed count.
The posterior means are then converted to proportions.
In the empirical Bayes step, the counts are smoothed by assuming a log-linear relationship between frequencies and frequencies of frequencies.
The fundamentals of the algorithm are from Good (1953).
Gale and Sampson (1995) proposed a simplied algorithm with a rule for switching between the observed and smoothed frequencies, and it is Gale and Sampson's simplified algorithm that is implemented here.
The number of zero values in \code{x} is not used as part of the algorithm, but is returned by this function.

Sampson gives a C code version on his webpage at
\url{http://www.grsampson.net/RGoodTur.html}
which gives identical results to this function.

\code{goodTuringPlot} plots log-probability (i.e., log frequencies of frequencies) versus log-frequency.

\code{goodTuringProportions} runs \code{goodTuring} on each column of data, then
uses the results to predict the proportion of each gene in each library.
}

\references{
Gale, WA, and Sampson, G (1995).
Good-Turing frequency estimation without tears.
\emph{Journal of Quantitative Linguistics} 2, 217-237.

Good, IJ (1953).
The population frequencies of species and the estimation of population parameters.
\emph{Biometrika} 40, 237-264.
}

\author{Aaron Lun and Gordon Smyth, adapted from Sampson's C code from \url{http://www.grsampson.net/RGoodTur.html}}

\examples{
#  True means of observed species
lambda <- rnbinom(10000,mu=2,size=1/10)
lambda <- lambda[lambda>1]

#  Oberved frequencies
Ntrue <- length(lambda)
x <- rpois(Ntrue, lambda=lambda)
freq <- goodTuring(x)
goodTuringPlot(x)
}

\concept{Normalization}
