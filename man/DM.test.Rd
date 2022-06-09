\name{DM.test}
\alias{DM.test}

\title{Computes Diebold-Mariano Test for the Equal Predictive Accuracy.}

\description{
This function computes Diebold-Mariano test for the equal predictive accuracy. The null hypothesis of this test is that two forecasts have the same accuracy. The alternative hypothesis can be specified as ''Both forecasts have different accuracy'', ''The first forecast is less accurate than the second forecast'', or ''The first forecast is more accurate than the second forecast''. 
}

\usage{
DM.test(f1,f2,y,loss.type="SE",h,c=FALSE,H1="same")
}

\arguments{
\item{f1}{\code{\link[base]{vector}} of the first forecast}
\item{f2}{\code{\link[base]{vector}} of the second forecast}
\item{y}{\code{\link[base]{vector}} of the real values of the modelled time-series}
\item{loss.type}{method to compute the loss function, \code{loss.type="SE"} will use squared errors, \code{loss.type="AE"} will use absolute errors, \code{loss.type="SPE"} will use squred proportional error (useful if errors are heteroskedastic), \code{loss.type="ASE"} will use absolute scaled error, if \code{loss.type} will be specified as some \code{\link[base]{numeric}}, then the function of type \code{exp(loss.type*errors)-1-loss.type*errors} will be used (useful when it is more costly to underpredict \code{y} than to overpredict), if not specified \code{loss.type="SE"} is used}
\item{h}{\code{\link[base]{numeric}} dentoing that the forecast \code{h}-steps ahead are evaluated, if not specified \code{h=1} is used}
\item{c}{\code{\link[base]{logical}} indicating if Harvey-Leybourne-Newbold correction for small samples should be used, if not specified \code{c=FALSE} is used}
\item{H1}{alternative hypothesis, \code{H1="same"} for ''both forecasts have different accuracy'', \code{H1="more"} for ''the first forecast is more accurate than the second forecast'', \code{H1="less"} for ''the first forecast is less accurate than the second forecast'', if not specified \code{H1="same"} is taken}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\code{h}, forecast horizon used}
\item{alternative}{alternative hypothesis of the test}
\item{p.value}{p-value}
\item{method}{name of the test}
\item{data.name}{names of the tested time-series}
}

\examples{
data(MDMforecasts)
ts <- MDMforecasts$ts
forecasts <- MDMforecasts$forecasts
DM.test(f1=forecasts[,1],f2=forecasts[,2],y=ts,loss="SE",h=1,c=FALSE,H1="same")
}

\references{
Diebold, F.X., Mariano, R. 1995. Comparing predictive accuracy. \emph{Journal of Business and Economic Statistics} \bold{13}, 253--265.

Harvey, D., Leybourne, S., Newbold, P., 1997. Testing the equality of prediction mean squared errors. \emph{International Journal of Forecasting} \bold{13}, 281--291.

Hyndman, R.J., Koehler, A.B. 2006. Another look at measures of forecast accuracy. \emph{International Journal of Forecasting} \bold{22}, 679--688.

Taylor, S. J., 2005. \emph{Asset Price Dynamics, Volatility, and Prediction},  Princeton University Press.

Triacca, U., 2018. \emph{Comparing Predictive Accuracy of Two Forecasts}, \url{https://www.lem.sssup.it/phd/documents/Lesson19.pdf}.
}
