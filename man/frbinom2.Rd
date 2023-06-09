% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frbinom2.R
\name{dfrbinom2}
\alias{dfrbinom2}
\alias{pfrbinom2}
\alias{qfrbinom2}
\alias{rfrbinom2}
\title{Fractional binomial distribution II}
\usage{
dfrbinom2(x, size, h, c, la = c/2, start = FALSE)

pfrbinom2(x, size, h, c, la = c/2, start = FALSE)

qfrbinom2(p, size, h, c, la = c/2, start = FALSE)

rfrbinom2(n, size, h, c, la = c/2, start = FALSE)
}
\arguments{
\item{x}{A numeric vector specifying values of the fractional binomial random variable II at which the pmf or cdf is computed.}

\item{size}{A number specifying the total number of trials.}

\item{h}{A number specifying the strength of the dependency among trials.}

\item{c}{A number specifying the dispersion of distributions.}

\item{la}{A number related to the probability of success in a trial; the default value is \code{c}/2.}

\item{start}{logical; if TRUE, the starting point is changed after the first success in the generalized Bernoulli process II.  The default is FALSE.}

\item{p}{A numeric vector specifying probabilities at which quantiles of the fractional binomial distribution II are computed.}

\item{n}{A number of random variables to be simulated.}
}
\value{
A numeric vector of random variables (rfrbinom2) or pmf (dfrbinom2) or cdf (pfrbinom2) or quantile (qfrbinom2) of the fractional binomial distribution II.
}
\description{
Generating random variables and computing density, cumulative distribution, and quantiles of the fractional binomial distribution II with the parameters \code{size}, \code{h}, \code{c}, \code{la}.
}
\details{
In the fractional binomial distribution II,  the number of successes is counted in the generalized Bernoulli process II (GBP II).
In GBP II, each trial has the constant probability of success \code{la}*\code{size}^\{2\code{h}-2\}, and the strength of the dependency among the trials is determined by the parameters, \code{h} and \code{c}.
The parameters \code{h}=\eqn{H}, \code{c}=\eqn{c}, \code{la}=\eqn{\lambda} should satisfy
\eqn{H \in (0.5,1)}, \eqn{0< c< 2^{2H-2}}, and  \eqn{0<\lambda<c}.
With the number of trials (\code{size}) =\eqn{n}, the mean of the fractional binomial random variable II is \eqn{E(X)=n\lambda^{2H-1}},
and the k-th moment is asymptotically proportional to \eqn{n^{(2H-1)k}} for \eqn{k\geq 2}.
}
\section{Functions}{
\itemize{
\item \code{dfrbinom2()}: The pmf of fractional binomial distribution II.

\item \code{pfrbinom2()}: The cdf of fractional binomial distribution II.

\item \code{qfrbinom2()}: The quantiles  of fractional binomial distribution II.

\item \code{rfrbinom2()}: The generation of fractional binomial random variables II.

}}
\examples{
# 10 random variables of a fractional binomial distribution II.
rfrbinom2(n=10, size=50, h=.8, c=.2)
# The probability that the fractional binomial random variable II equals 22.
dfrbinom2(x=22, size=50, h=.8, c=.2)
# The probability that the fractional binomial random variable II is less than or equal to 22.
pfrbinom2(x=22, size=50, h=.8, c=.2)
 # The 80th percentile of the fractional binomial distribution.
qfrbinom2(p=.8, size=50, h=.8, c=.2)
}
\references{
Lee, J. (2023). Generalized Bernoulli process and fractional Poisson process. arXiv:2209.01516.
}
