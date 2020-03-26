\name{TB_AR_test}
\alias{TB_AR_test}

\title{Computes Tiao-Box Test for Autocorrelation.}

\description{
This function computes Tiao-Box test for autocorrelation, i.e, coefficient of \code{p}-th lag in VAR(p) model. Its null hypothesis is that \code{p}-th lag is not essential. The alternative hypothesis is that it is essential. 
}

\usage{
TB_AR_test(d,p)
}

\arguments{
\item{d}{\code{\link[base]{matrix}} of time-series, assumed to be the stationary VARMA type, columns correspond to time index, and rows to different time-series}
\item{p}{\code{\link[base]{numeric}} indicating a lag length beyond which we are willing to assume that the autocorrelation is essentially zero}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\code{q}, a lag length}
\item{alternative}{alternative hypothesis of the test}
\item{p.value}{p-value}
\item{method}{name of the test}
\item{data.name}{name of the tested time-series}
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
l <- loss(realized=ts,evaluated=forecasts,loss.type="SE")
d <- d_t(l)
TB_AR_test(d=d,p=10)
}

\references{
Tiao, G.C., Box, G.E.P. 1981. Modeling multiple times series with applications. \emph{Journal of the American Statistical Association} \bold{76}, 802--816.
}
