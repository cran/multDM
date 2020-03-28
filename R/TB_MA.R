
### Tiao-Box procedure to determine the lag in MA(q)
###
### G.C. Tiao, G.E.P. Box, (1981),
### Modeling multiple times series with applications,
### Journal of the American Statistical Association 76, 802-816.
###
### d - kxP matrix of observed values
###       P - lenght of time-series
###
### q.max - maximum number of lag to be considered

TB_MA <- function(d,q.max)
  {
    n <- deparse(substitute(d))

    dbar <- rowMeans(d)
    sds <- t(apply(d,1,sd))
    sds <- crossprod(sds)
    
    G <- function(h)
      {
        SCM <- matrix(0,nrow(d),nrow(d))
        for (t in (h+1):ncol(d))
          {
            SCM <- SCM + (d[,t,drop=FALSE] - dbar) %*% t(d[,t-h,drop=FALSE] - dbar)
          }
        SCM <- SCM / ncol(d)
        SCM <- SCM / sds
        
        return(SCM)
      }
    
    m1 <- lapply(seq(from=1,to=q.max),G)
    b <- 2*ncol(d)^(-0.5)
    
    C.check <- function(i)
      {
        return(any(m1[[i]] > b | m1[[i]] < (-b)))
      }
    
    m2 <- unlist(lapply(seq(from=1,to=length(m1)),C.check))
    
    q <- which(m2==TRUE)
    if (length(q) > 0)
      {
        q <- max(which(m2==TRUE))
      }
    else
      {
        q <- 0
      }

    return(q)
  }
