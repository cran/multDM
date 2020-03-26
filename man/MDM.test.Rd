\name{MDM.test}
\alias{MDM.test}

\title{Computes Multivariate Diebold-Mariano Test for the Equal Predictive Accuracy of Two or More Non-nested Forecasting Models.}

\description{
This function computes multivariate Diebold-Mariano test for the equal predictive accuracy of two or more non-nested forecasting models. The null hypothesis of this test is that the evaluated forecasts have the same accuracy. The alternative hypothesis is that Equal predictive accuracy (EPA) does not hold. 
}

\usage{
MDM.test(realized,evaluated,q,statistic="Sc",loss.type="SE")
}

\arguments{
\item{realized}{\code{\link[base]{vector}} of the real values of the modelled time-series}
\item{evaluated}{\code{\link[base]{matrix}} of the forecasts, columns correspond to time index, rows correspond to different models}
\item{q}{\code{\link[base]{numeric}} indicating a lag length beyond which we are willing to assume that the autocorrelation of loss differentials is essentially zero}
\item{statistic}{\code{statistic="S"} for the basic version of the test, and \code{statistic="Sc"} for the finite-sample correction, if not specified \code{statistic="Sc"} is used}
\item{loss.type}{method to compute the loss function, \code{loss.type="SE"} will use squared errors, \code{loss.type="AE"} will use absolute errors, \code{loss.type="SPE"} will use squred proportional error (useful if errors are heteroskedastic), \code{loss.type="ASE"} will use absolute scaled error, if \code{loss.type} will be specified as some \code{\link[base]{numeric}}, then the function of type \code{exp(loss.type*errors)-1-loss.type*errors} will be used (useful when it is more costly to underpredict \code{realized} than to overpredict), if not specified \code{loss.type="SE"} is used}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\code{q}, a lag length}
\item{alternative}{alternative hypothesis of the test}
\item{p.value}{p-value}
\item{method}{name of the test}
\item{data.name}{names of the tested objects}
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
MDM.test(realized=ts,evaluated=forecasts,q=10,statistic="S",loss.type="AE")
}

\references{
Mariano R.S., Preve, D., 2012. Statistical tests for multiple forecast comparison. \emph{Journal of Econometrics} \bold{169}, 123--130.
}
