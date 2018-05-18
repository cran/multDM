\name{TB_MA}
\alias{TB_MA}

\title{Checks for a Lag in VMA Process with Tiao-Box Procedure.}

\description{
This function helps to find a lag in stationary VMA process with Tiao-Box procedure, i.e., the lag length beyond which we are willing to assume that the autocorrelation is essentially zero. 
}

\usage{
TB_MA(d,q.max)
}

\arguments{
\item{d}{\code{\link[base]{matrix}} of time-series, assumed to be the stationary VARMA type, columns correspond to time index, and rows to different time-series}
\item{q.max}{\code{\link[base]{numeric}} indicating the maximum number of lag to be considered}
}

\value{\code{\link[base]{numeric}} indicating the found lag length
}

\details{
The function searches for correlations smaller than \eqn{-2n^{-0.5}} or higher than \eqn{2n^{-0.5}}, where \eqn{n} is the lenght of the time-series. 
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
l <- loss(realized=ts,evaluated=forecasts,loss.type="SE")
d <- d_t(l)
TB_MA(d=d,q.max=10)
}

\references{
Tiao, G.C., Box, G.E.P. 1981. Modeling multiple times series with applications. \emph{Journal of the American Statistical Association} \bold{76}, 802--816.
}
