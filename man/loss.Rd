\name{loss}
\alias{loss}

\title{Computes Loss Function.}

\description{
This function computes various loss functions for given realized values of time-series and a collection of forecasts. 
}

\usage{
loss(realized,evaluated,loss.type)
}

\arguments{
\item{realized}{\code{\link[base]{vector}} of the real values of the modelled time-series}
\item{evaluated}{\code{\link[base]{matrix}} of the forecasts, columns correspond to time index, rows correspond to different models}
\item{loss.type}{method to compute the loss function, \code{loss.type="SE"} will use squared errors, \code{loss.type="AE"} will use absolute errors, \code{loss.type="SPE"} will use squred proportional error (useful if errors are heteroskedastic), \code{loss.type="ASE"} will use absolute scaled error, if \code{loss.type} will be specified as some \code{\link[base]{numeric}}, then the function of type \code{exp(loss.type*errors)-1-loss.type*errors} will be used (useful when it is more costly to underpredict \code{realized} than to overpredict)}
}

\value{\code{\link[base]{matrix}} with columns corresponding to time index and rows to different models 
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
l <- loss(realized=ts,evaluated=forecasts,loss.type="SE")
}

\references{
Hyndman, R.J., Koehler, A.B. 2006. Another look at measures of forecast accuracy. \emph{International Journal of Forecasting} \bold{22}, 679--688.

Taylor, S. J., 2005. \emph{Asset Price Dynamics, Volatility, and Prediction},  Princeton University Press.

Triacca, U., 2018. \emph{Comparing Predictive Accuracy of Two Forecasts}, \url{http://www.phdeconomics.sssup.it/documents/Lesson19.pdf}.
}
