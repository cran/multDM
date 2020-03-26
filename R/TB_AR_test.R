
### Tiao-Box test of the lag in AR(p)
###
### G.C. Tiao, G.E.P. Box, (1981),
### Modeling multiple times series with applications,
### Journal of the American Statistical Association 76, 802-816.
###
### d - kxP matrix of observed values
###       P - lenght of time-series
###
### p - lag order tested

TB_AR_test <- function(d,p)
  {
    n <- deparse(substitute(d))

    m1 <- VAR(t(d),p=p,output=FALSE)$Sigma
    m2 <- VAR(t(d),p=p-1,output=FALSE)$Sigma
    
    m1 <- det(m1)
    m2 <- det(m2)
    U <- m1 / m2
    
    N <- ncol(d) - p - 1
    M <- -(N - 0.5 - p * nrow(d)) * log(U)
    pval <- pchisq(q=M,df=(nrow(d))^2,lower.tail=TRUE)

    names(M) <- "statistic"
    names(p) <- "lag length"

    ret <- list(M,p,"Higher lag is essential.",pval,"Tiao-Box likelihood test",n)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    return(ret)
  }
