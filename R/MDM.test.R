
### computes multivariate DM test for equal predictive accuracy (EPA) 
### of two or more non-nested forecasting models,
### as in R.S. Mariano and D. Preve (2012) 
### Statistical tests for multiple forecast comparison,
### Journal of Econometrics 169, 123-130
###
### d - kxP matrix with observations of d_t
###       k+1 - number of considered models
###       P - lenght of time-series
###       d_t = g(y.hat_1t - y.real_t ) - g(y.hat_2t - y.real_t), i.e.,
###       d_jt = g(e_jt) - g(e_j+1,t), j = 1,...,k,
###       where g is the loss function
###
### q - order of the VMA representation of d_t,
###     i.e., a lag length beyond which we are willing 
###     to assume that the correlation between d_t and d_t-h 
###     is essentially zero
###
### statistic - "S" for the basic version,
###             "Sc"  for the finite-sample correction
###
### the null hypothesis of EPA is E(d_t)=0


.in.MDM.test <- function(d,q,statistic)
  {
    n <- deparse(substitute(d))
    
    G <- function(d,h)
      {
        SCM <- matrix(0,nrow(d),nrow(d))
        for (t in (h+1):ncol(d))
          {
            SCM <- SCM + (d[,t,drop=FALSE] - dbar) %*% t(d[,t-h,drop=FALSE] - dbar)
          }
        SCM <- SCM / ncol(d)
        
        return(SCM)
      }

    O <- function(d,q)
      {
        SLRV <- G(d,0)
        if (q>0)
          {
            for (h in 1:q)
              {
                TEMP <- G(d,h)
                SLRV <- SLRV + TEMP + t(TEMP)
              }
          }
          
        return(SLRV)
      }

    dbar <- rowMeans(d)
    c <- 1 - (1 + 2 * q) / ncol(d) + q * (q+1) / (ncol(d))^2
    Om <- O(d,q)
    s <- as.vector(dbar) / sqrt(as.vector(abs(diag(Om))) / ncol(d))  
    S <- ncol(d) * t(dbar) %*% solve(Om) %*% dbar
    if (statistic=="Sc") 
      {
        S <- c * S
      }
      
    pval <- pchisq(q=S,df=nrow(d),lower.tail=TRUE)
    
    names(S) <- "statistic"
    names(q) <- "lag length"

    ret <- list(S,q,"Equal predictive accuracy does not hold.",pval,"multivariate Diebold-Mariano test",n)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    ret <- list(s,ret)
    return(ret)
  }


MDM.test <- function(realized,evaluated,q,statistic="Sc",loss.type="SE")
  {
    n <- paste(deparse(substitute(realized)),deparse(substitute(evaluated)),sep=" and ")
    l <- loss(realized=realized,evaluated=evaluated,loss.type=loss.type)
    d <- d_t(l)
    out <- .in.MDM.test(d=d,q=q,statistic=statistic)[[2]]
    out$data.name <- n
    return(out)
  }


##################################################################


### computes loss function for a set of forecasts
###
### realized - vector of realized values
###
### evaluated - (k+1)xP matrix of forecasts
###               P - length of time-series
###               k+1 - number of models
### 
### loss.type - "SE" for squared errors,
###             "AE" for absolute errors,
###              "SPE" for squared proportional error
###                (useful if errors are heteroskedastic)
###                see S. J. Taylor, 2005. Asset Price Dynamics, Volatility, and Prediction,
###                Princeton University Press,
###             "ASE" for absolute scaled error (R. J. Hyndman, A. B. Koehler, 2006,
###               Another Look at Measures of Forecast Accuracy,
###               International Journal of Forecasting volume 22, 679-688, 
###             positive numeric value for loss function of type
###              exp(loss*errors)-1-loss*errors
###              (useful when it is more costly to underpredict y than to overpredict)
###              see U. Triacca,  Comparing Predictive Accuracy of Two Forecasts,
###              http://www.phdeconomics.sssup.it/documents/Lesson19.pdf 


loss <- function(realized,evaluated,loss.type)
  {
    e <- evaluated
    for (i in 1:nrow(evaluated))
      {
        e[i,] <- realized - as.vector(evaluated[i,])
      }
    
    if (loss.type=="SE") 
      {
        e <- e^2
      }
    if (loss.type=="AE") 
      {
        e <- abs(e)
      }
    if (loss.type=="SPE") 
      {
        for (i in 1:nrow(e))
          {
            e[i,] <- (as.vector(e[i,]) / as.vector(evaluated[i,]))^2
          }
      }
    if (loss.type=="ASE") 
      {
       for (i in 1:nrow(e))
          {
            e[i,] <- abs(as.vector(e[i,])) / mean(abs(realized-c(NA,realized[-length(realized)]))[-1])   
          }
       e <- e[,-1,drop=FALSE]
      }
    if (is.numeric(loss.type)) 
      {
        e <- exp(loss.type*e)-1-loss.type*e
      }
    
    return(e)
  }


### computes d_t based on loss function

d_t <- function(e)
  {
    for (i in 1:(nrow(e)-1))
      {
        e[i,] <- (as.vector(e[i,]) - as.vector(e[i+1,]))
      }   
    e <- e[-nrow(e),,drop=FALSE] 

    return(e)
  }


MDM.selection <- function(realized,evaluated,q,alpha,statistic="Sc",loss.type="SE")
  {
    p <- 0
    e <- evaluated
    models <- (1:nrow(e))
    if (is.null(rownames(e))) { rownames(e) <- seq(1:nrow(e)) }
    n.mods <- nrow(e)
    
    while(p<alpha && length(models)>1)
      {
        d <- loss(realized=realized,evaluated=e,loss.type=loss.type)
        d <- d_t(d)
        mdm <- .in.MDM.test(d=d,q=q,statistic=statistic)
        p <- mdm[[2]]$p.value
        j <- which.max(abs(mdm[[1]]))
        if (mdm[[1]][j]>0)
          {
            models <- models[-j]
            j.drop <- j
          }
        else
          {
            models <- models[-(j+1)]
            j.drop <- j+1
          }
        e <- e[-j.drop,,drop=FALSE]
        models <- (1:nrow(e))
      }
    
    n.mods <- n.mods - nrow(e)
    ret <- matrix(NA,ncol=3,nrow=nrow(e))
    rownames(ret) <- rownames(e)
    if (j.drop==(nrow(e)+1))
      {
        ret[,2] <- mdm[[1]][-j.drop]
        ret[,1] <- sort(ret[,2],index.return=TRUE)$ix
      }
    else
      {
        ret[-nrow(ret),2] <- mdm[[1]][-j.drop]
        ret[-nrow(ret),1] <- sort(ret[-nrow(ret),2],index.return=TRUE)$ix
      }
    ret[,3] <- rowMeans(loss(realized=realized,evaluated=e,loss.type=loss.type))
    colnames(ret) <- c("Rank",statistic,"Mean loss")

    ret <- list(ret,as.numeric(p),alpha,n.mods)
    names(ret) <- c("outcomes","p.value","alpha","eliminated")
    class(ret) <- "MDM"
    return(ret)
  }
 
  
print.MDM <- function(x, ...)
  {
    if (x[[2]]>x[[3]])
      {
        cat("####################################################")
        cat("\n")
        cat("Models with outstanding predictive ability:")
        cat("\n")
        cat("\n")
        print(round(x[[1]],digits=4),na.print="")
        cat("\n")
        cat("p-value: ")
        cat(round(as.numeric(x[[2]]),digits=4))
        cat("\n")
        cat("\n")
        cat("Number of eliminated models: ",x[[4]])
        cat("\n")
        cat("####################################################")
      }
    else
      {
        cat("No models with outstanding predictive ability were found.")
      }
  }

