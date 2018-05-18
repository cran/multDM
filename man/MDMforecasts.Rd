\name{MDMforecasts}
\alias{MDMforecasts}

\title{Sample Data.}

\description{Sample artificial data.}

\usage{data(MDMforecasts)}

\format{
\itemize{\code{MDMforecasts} is \code{\link[base]{list}} object such that
\item \code{MDMforecasts$ts} is \code{\link[base]{vector}} of time-series which is of interest to model
\item \code{MDMforecasts$forecasts} is \code{\link[base]{matrix}} of 20 different forecasts of \code{ts} from 20 different forecasting models, each row represents different forecast and time is indexed by columns
}
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
MDM.test(realized=ts,evaluated=forecasts,q=10,statistic="S",loss.type="AE")
}
