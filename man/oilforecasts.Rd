\name{oilforecasts}
\alias{oilforecasts}

\title{Sample Data from Crude Oil Price Forecasting.}

\description{Forecasts obtained from various methods applied to crude oil price.}

\usage{data(oilforecasts)}

\format{
\itemize{\code{oilforecasts} is \code{\link[base]{matrix}} object such that its rows correspond to forecasts from various methods, i.e.,
\item \code{REALIZED} is the forecasted time-series,
\item \code{DMA.DOW} is the forecast from Dynamic Model Averaging with the dynamic Occam's window,
\item \code{BMA.DOW} is the forecast from Bayesian Model Averaging with the dynamic Occam's window,
\item \code{DMA.1V} is the forecast from Dynamic Model Averaging applied only to one-variable models,
\item \code{BMA.1V} is the forecast from Bayesian Model Averaging applied only to one-variable models,
\item \code{DMS.1V} is the forecast from Dynamic Model Selection applied only to one-variable models,
\item \code{BMS.1V} is the forecast from Bayesian Model Selection applied only to one-variable models,
\item \code{TVP} is the forecast from Time-Varying Parameters regression,
\item \code{LASSO} is the forecast from LASSO regression,
\item \code{RIDGE} is the forecast from RIDGE regression,
\item \code{DYN.EL.NET} is the forecast from the elastic net regression, with the elastic net mixing parameter changing with time index,
\item \code{LARS} is the forecast from the least-angle regression,
\item \code{B.LASSO} is the forecast from the Bayesian LASSO regression,
\item \code{B.RIDGE} is the forecast from the Bayesian RIDGE regression,
\item \code{ARIMA} is the forecast from the best ARIMA model according to AIC,
\item \code{NAIVE} is the naive forecast, i.e., the last observation is taken as a one-step ahead prediction,
\item \code{MA} is the moving average.
}
}

\details{
The data were taken from Juvenal and Petrella (2015). They cover the period between 1971 and 2009, and are in quarterly freqency. Time-series with missing observations were excluded from the original data set, resulting finally in 127 explanatory variables, instead of 150 in the original data set. In particular, the excluded time-series are the ones which start date is after 1971. The dependent time-series is the average world price of oil taken in logarithmic differences. The independent time-series represent various stationarity tranformations of macroeconomic and financial variables of the G7 countries and from the oil market, global economic activity and various commodity prices. The details of the original data set are given in the paper by Juvenal and Petrella (2015). The forecasting with various models, based on this data set, was done by the author of this package, just to provide some more concrete example set of forecasts. The independent variables were taken in the first lags. The forgetting parameters in DMA/DMS models were set to 0.99, resulting in the effective rolling window size of 100. Therefore, such a window was taken for the moving average. LASSO and RIDGE (also in the Bayesian versions), the elastic net, the least-angle regression and ARIMA models were estimated in rolling windows of the size of 100 observations. First 100 observations were excluded, and \code{oilforecasts} consists of the remaining last observations. The estimations were done with the following packages \code{fDMA}, \code{forecast}, \code{glmnet}, \code{lars} and \code{monomvn}. 
}

\examples{
data(oilforecasts)
ts <- oilforecasts[1,]
forecasts <- oilforecasts[-1,]
l <- loss(realized=ts,evaluated=forecasts,loss.type="SE")
d <- d_t(l)
q <- TB_MA(d=d,q.max=4)
MDM.selection(realized=ts,evaluated=forecasts,q=q,alpha=0.1,statistic="Sc",loss.type="SE")
}

\references{
Drachal, K. 2017. fDMA: Dynamic Model Averaging and Dynamic Model Selection for continuous outcomes. \url{https://CRAN.R-project.org/package=fDMA}

Friedman, J., Hastie, T., Tibshirani, R. 2010. Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software} \bold{33}, 1--22.

Gramacy, R.B. 2017. monomvn: Estimation for Multivariate Normal and Student-t Data with Monotone Missingness. \url{https://CRAN.R-project.org/package=monomvn}

Hastie, T., Efron, B. 2013. lars: Least Angle Regression, Lasso and Forward Stagewise. \url{https://CRAN.R-project.org/package=lars}

Hyndman, R.J., Khandakar. Y. 2008. Automatic time series forecasting: the forecast package for R. \emph{Journal of Statistical Software} \bold{26}, 1--22.

Juvenal, L., Petrella, I. 2015. Speculation in the oil market. \emph{Journal of Applied Econometrics} \bold{30}, 612--649.
}
