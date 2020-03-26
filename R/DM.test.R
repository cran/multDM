
### Diebold-Mariano test
###
### f1 - forecast 1 (vector)
### f2 - forecast 2 (vector)
### y - observed values (vector)
### loss.type - "SE" for squared errors,
###             "AE" for absolute errors,
###             "SPE" for squared proportional error
###               (useful if errors are heteroskedastic)
###               see S. J. Taylor, 2005. Asset Price Dynamics, Volatility, and Prediction,
###               Princeton University Press,
###             "ASE" for absolute scaled error (R. J. Hyndman, A. B. Koehler, 2006,
###               Another Look at Measures of Forecast Accuracy,
###               International Journal of Forecasting volume 22, 679-688, 
###             positive numeric value for loss function of type
###               exp(loss*errors)-1-loss*errors
###               (useful when it is more costly to underpredict y than to overpredict)
###               see U. Triacca,  Comparing Predictive Accuracy of Two Forecasts,
###               http://www.phdeconomics.sssup.it/documents/Lesson19.pdf 
### h - dentoes that the forecast h-steps ahead is evaluated
### c - TRUE for using Harvey-Leybourne-Newbold correction for small samples,
###     see D. Harvey, S. Leybourne, P. Newbold, 1997,
###     Testing the Equality of Prediction Mean Squared Errors,
###     International Journal of Forecasting 13, 281-291
### H1 - alternative hypothesis
###      f1 and f2 have the "same" accuracy,
###      f1 is "more" accurate than f2,
###      f1 is "less" accurate than f2

DM.test <- function(f1,f2,y,loss.type="SE",h=1,c=FALSE,H1="same")
  {
    n <- paste(deparse(substitute(f1)),deparse(substitute(f2)),deparse(substitute(y)),sep=" and ")
    n1 <- deparse(substitute(f1))
    n2 <- deparse(substitute(f2))

    e1 <- f1 - y
    e2 <- f2 - y
    if (loss.type=="SE") 
      {
        g1 <- e1^2
        g2 <- e2^2
      }
    if (loss.type=="AE") 
      {
        g1 <- abs(e1)
        g2 <- abs(e2)
      }
    if (loss.type=="SPE") 
      {
        g1 <- ((y-f1)/f1)^2
        g2 <- ((y-f2)/f2)^2
      }
    if (loss.type=="ASE") 
      {
        g1 <- abs(e1[-1])/mean(abs((y-c(NA,y[-length(y)]))[-1]))
        g2 <- abs(e2[-1])/mean(abs((y-c(NA,y[-length(y)]))[-1]))
      }
   if (is.numeric(loss.type)) 
      {
        g1 <- exp(loss.type*e1)-1-loss.type*e1
        g2 <- exp(loss.type*e2)-1-loss.type*e2
      }
    d <- g1-g2
    T <- length(d)
    dbar <- mean(d)
    gammahat <- function(k)
      {
        temp1 <- d-dbar
        temp2 <- rep(NA,abs(k))
        temp2 <- c(temp2,temp1)
        temp2 <- temp2[1:T]
        temp2 <- temp2-dbar
        temp <- temp1*temp2
        temp <- temp[(1+abs(k)):T]
        temp <- sum(temp)/T
        return(temp)
      }
    if (h>1) 
      {
        gdk <- lapply(seq(from=1,to=h-1,by=1),gammahat)
        gdk <- unlist(gdk)
      }
    else
      {
        gdk <- 0
      }
    gdk <- gammahat(0)+2*sum(gdk)
    DM <- dbar/sqrt(gdk/T)

    if (H1=="same") { pval <- 2 * min(pnorm(DM,lower.tail=FALSE),1 - pnorm(DM,lower.tail=FALSE)) }
    if (H1=="less") { pval <- pnorm(DM,lower.tail=FALSE) }
    if (H1=="more") { pval <- 1 - pnorm(DM,lower.tail=FALSE) }

    if (c)
      {
        DM <- DM*sqrt((T+1-2*h+h*(h-1))/T)
        if (H1=="same") { pval <- 2 * min(pt(q=DM,df=T-1,lower.tail=FALSE),1 - pt(q=DM,df=T-1,lower.tail=FALSE)) }
        if (H1=="less") { pval <- pt(q=DM,df=T-1,lower.tail=FALSE) }
        if (H1=="more") { pval <- 1 - pt(q=DM,df=T-1,lower.tail=FALSE) }
      }
    
    names(DM) <- "statistic"
    names(h) <- "forecast horizon"
    if (H1=="same") { alt <- "Forecast f1 and f2 have different accuracy." }
    if (H1=="less") { alt <- "Forecast f1 is less accurate than f2." }
    if (H1=="more") { alt <- "Forecast f1 is more accurate than f2." }

    ret <- list(DM,h,paste(alt),pval,"Diebold-Mariano test",n)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }
