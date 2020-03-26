\name{MDM.selection}
\alias{MDM.selection}

\title{Selects Models with Outstanding Predictive Ability basing on Multivariate Diebold-Mariano Test.}

\description{
This function selects models with outstanding predictive ability basing on multivariate Diebold-Mariano test \code{\link{MDM.test}}. 
}

\usage{
MDM.selection(realized,evaluated,q,alpha,statistic="Sc",loss.type="SE")
}

\arguments{
\item{realized}{\code{\link[base]{vector}} of the real values of the modelled time-series}
\item{evaluated}{\code{\link[base]{matrix}} of the forecasts, columns correspond to time index, rows correspond to different models}
\item{q}{\code{\link[base]{numeric}} indicating a lag length beyond which we are willing to assume that the autocorrelation of loss differentials is essentially zero}
\item{alpha}{\code{\link[base]{numeric}} indicating a significance level for multivariate Diebold-Mariano tests}
\item{statistic}{\code{statistic="S"} for the basic version of the test, and \code{statistic="Sc"} for the finite-sample correction, if not specified \code{statistic="Sc"} is used}
\item{loss.type}{method to compute the loss function, \code{loss.type="SE"} will use squared errors, \code{loss.type="AE"} will use absolute errors, \code{loss.type="SPE"} will use squred proportional error (useful if errors are heteroskedastic), \code{loss.type="ASE"} will use absolute scaled error, if \code{loss.type} will be specified as some \code{\link[base]{numeric}}, then the function of type \code{exp(loss.type*errors)-1-loss.type*errors} will be used (useful when it is more costly to underpredict \code{realized} than to overpredict), if not specified \code{loss.type="SE"} is used}
}

\value{class \code{MDM} object, \code{\link[base]{list}} of 
\item{outcomes}{\code{\link[base]{matrix}} with mean losses for the selected models, statistics corresponding to losses differentials and ranking of these statistics}
\item{p.value}{\code{\link[base]{numeric}} of p-value from the procedure, i.e., p-value of multivariate Diebold-Mariano test from the last step}
\item{alpha}{\code{alpha}, i.e., the chosen significance level}
\item{eliminated}{\code{\link[base]{numeric}} indicating the number of eliminated models}
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
MDM.selection(realized=ts,evaluated=forecasts,q=10,alpha=0.1,statistic="S",loss.type="AE")
}

\references{
Mariano R.S., Preve, D., 2012. Statistical tests for multiple forecast comparison. \emph{Journal of Econometrics} \bold{169}, 123--130.
}
