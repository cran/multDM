\name{d_t}
\alias{d_t}

\title{Computes Loss Differential.}

\description{
This function computes loss differential, i.e., differences between losses from \eqn{k+1}-th and \eqn{k}-th models.
}

\usage{
d_t(e)
}

\arguments{
\item{e}{\code{\link[base]{matrix}} of loss functions, columns correspond to time index, and rows to different models}
}

\value{\code{\link[base]{matrix}} of loss differentials
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
l <- loss(realized=ts,evaluated=forecasts,loss.type="SE")
d <- d_t(l)
}

\references{
Mariano R.S., Preve, D., 2012. Statistical tests for multiple forecast comparison. \emph{Journal of Econometrics} \bold{169}, 123--130.
}
