\name{print.MDM}
\alias{print.MDM}

\title{Prints \code{MDM} Object.}

\description{The function prints selected outcomes obtained from \code{\link{MDM.selection}}.}

\usage{\method{print}{MDM}(x, ...)}

\arguments{
\item{x}{an object of \code{MDM} class}
\item{...}{not used}
}

\details{
The function prints models with outstanding predictive ability, their mean loss function, statistics corresponding to their loss differentials (they are the number of the models less one), and orders these statistics. It also prints the p-value of the test and the number of eliminated models. If no models with outstanding predictive ability were found, the funtions prints such an information.
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
m <- MDM.selection(realized=ts,evaluated=forecasts,q=10,alpha=0.1,statistic="S",loss.type="AE")
print(m)
}
