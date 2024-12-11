## mgcv functions with outer.ok = TRUE added

#' @exportS3Method
smooth.construct.bs2.smooth.spec <- function(object,data,knots) {
## a B-spline constructor method function
  ## get orders: m[1] is spline order, 3 is cubic. m[2] is order of derivative in penalty.
  if (length(object$p.order)==1) m <- c(object$p.order,max(0,object$p.order-1)) 
  else m <- object$p.order  # m[1] - basis order, m[2] - penalty order
  if (is.na(m[1])) if (is.na(m[2])) m <- c(3,2) else m[1] <- m[2] + 1
  if (is.na(m[2])) m[2] <- max(0,m[1]-1)
  object$m <- object$p.order <- m
  if (object$bs.dim<0) object$bs.dim <- max(10,m[1]) ## default
  nk <- object$bs.dim - m[1] + 1  # number of interior knots
  if (nk<=0) stop("basis dimension too small for b-spline order")
  if (length(object$term)!=1) stop("Basis only handles 1D smooths")
  x <- data[[object$term]]    # find the data
  k <- knots[[object$term]]
  if (is.null(k)) { xl <- min(x);xu <- max(x) } else
  if (length(k)==2) { 
    xl <- min(k);xu <- max(k);
    if (xl>min(x)||xu<max(x)) stop("knot range does not include data")
  }
  if (!is.null(k)&&length(k)==4&&length(k)<nk+2*m[1]) {
    ## 4 knots supplied: lower prediction limit, lower data limit,
    ##   upper data limit, upper prediction limit
    k <- sort(k)
    dx <- (k[4]-k[1])/(nk-1)
    ko <- c(k[1]-dx*m[1],k[4]+dx*m[1]) ## limits for outer knots
    k <- c(seq(ko[1],k[1],length=m[1]+1),
       seq(k[2],k[3],length=max(0,nk-2)),
       seq(k[4],ko[2],length=m[1]+1))
    
  } else if (is.null(k)||length(k)==2) {
    xr <- xu - xl # data limits and range
    xl <- xl-xr*0.001;xu <- xu+xr*0.001;dx <- (xu-xl)/(nk-1) 
    k <- seq(xl-dx*m[1],xu+dx*m[1],length=nk+2*m[1])   
  } else {
    if (length(k)!=nk+2*m[1]) 
    stop(paste("there should be ",nk+2*m[1]," supplied knots"))
  }
  if (is.null(object$deriv)) object$deriv <- 0 
  object$X <- splines::spline.des(k,x,m[1]+1,x*0+object$deriv, outer.ok = TRUE)$design # get model matrix
  if (!is.null(k)) {
    if (sum(colSums(object$X)==0)>0) warning("there is *no* information about some basis coefficients")
  }  
  if (length(unique(x)) < object$bs.dim) warning("basis dimension is larger than number of unique covariates")
 
  ## now construct derivative based penalty. Order of derivate
  ## is equal to m, which is only a conventional spline in the 
  ## cubic case...
  
  object$knots <- k; 
  class(object) <- "Bspline2.smooth"  # Give object a class
  k0 <- k[m[1]+1:nk] ## the interior knots
  object$D <- object$S <- list()
  m2 <- m[2:length(m)] ## penalty orders
  if (length(unique(m2))<length(m2)) stop("multiple penalties of the same order is silly")
  for (i in 1:length(m2)) { ## loop through penalties
    object$deriv <- m2[i] ## derivative order of current penalty
    pord <- m[1]-m2[i] ## order of derivative polynomial 0 is step function
    if (pord<0) stop("requested non-existent derivative in B-spline penalty") 
    h <- diff(k0) ## the difference sequence...
    ## now create the sequence at which to obtain derivatives
    if (pord==0) k1 <- (k0[2:nk]+k0[1:(nk-1)])/2 else {
      h1 <- rep(h/pord,each=pord)
      k1 <- cumsum(c(k0[1],h1)) 
    } 
    dat <- data.frame(k1);names(dat) <- object$term 
    D <- Predict.matrix.Bspline2.smooth(object,dat) ## evaluate basis for mth derivative at the k1
    object$deriv <- NULL ## reset or the smooth object will be set to evaluate derivs in prediction! 
    if (pord==0) { ## integrand is just a step function...
      object$D[[i]] <- sqrt(h)*D
    } else { ## integrand is a piecewise polynomial...
      P <- solve(matrix(rep(seq(-1,1,length=pord+1),pord+1)^rep(0:pord,each=pord+1),pord+1,pord+1))
      i1 <- rep(1:(pord+1),pord+1)+rep(1:(pord+1),each=pord+1) ## i + j
      H <- matrix((1+(-1)^(i1-2))/(i1-1),pord+1,pord+1)
      W1 <- t(P)%*%H%*%P
      h <- h/2 ## because we map integration interval to to [-1,1] for maximum stability
      ## Create the non-zero diagonals of the W matrix... 
      ld0 <- rep(sdiag(W1),length(h))*rep(h,each=pord+1)
      i1 <- c(rep(1:pord,length(h)) + rep(0:(length(h)-1) * (pord+1),each=pord),length(ld0))
      ld <- ld0[i1] ## extract elements for leading diagonal
      i0 <- 1:(length(h)-1)*pord+1
      i2 <- 1:(length(h)-1)*(pord+1)
      ld[i0] <- ld[i0] + ld0[i2] ## add on extra parts for overlap
      B <- matrix(0,pord+1,length(ld))
      B[1,] <- ld
      for (k in 1:pord) { ## create the other diagonals...
        diwk <- sdiag(W1,k) ## kth diagonal of W1
        ind <- 1:(length(ld)-k)
        B[k+1,ind] <- (rep(h,each=pord)*rep(c(diwk,rep(0,k-1)),length(h)))[ind]  
      }
      ## ... now B contains the non-zero diagonals of W
      B <- bandchol(B) ## the banded cholesky factor.
      ## Pre-Multiply D by the Cholesky factor...
      D1 <- B[1,]*D
      for (k in 1:pord) {
        ind <- 1:(nrow(D)-k)
        D1[ind,] <- D1[ind,] + B[k+1,ind] * D[ind+k,]
      }
      object$D[[i]] <- D1
    }
    object$S[[i]] <- crossprod(object$D[[i]])
  }
  object$rank <- object$bs.dim-m2  # penalty rank 
  object$null.space.dim <- min(m2)    # dimension of unpenalized space 
 
  object
} ### end of B-spline constructor

#' @exportS3Method
Predict.matrix.Bspline2.smooth <- function(object,data) {
  object$mono <- 0
  object$m <- object$m - 1 ## for consistency with p-spline defn of m
  Predict.matrix.pspline2.smooth(object,data)
}



#' @exportS3Method
Predict.matrix.pspline2.smooth <- function(object,data)
# prediction method function for the p.spline smooth class
{ ##require(splines)
  m <- object$m[1]+1
  ## find spline basis inner knot range...
  ll <- object$knots[m+1];ul <- object$knots[length(object$knots)-m]
  m <- m + 1
  x <- data[[object$term]]
  n <- length(x)
  ind <- x<=ul & x>=ll ## data in range
  if (is.null(object$deriv)) object$deriv <- 0 
  if (sum(ind)==n) { ## all in range
    X <- splines::spline.des(object$knots,x,m,rep(object$deriv,n, outer.ok = TRUE))$design
  } else { ## some extrapolation needed 
    ## matrix mapping coefs to value and slope at end points...
    D <- splines::spline.des(object$knots,c(ll,ll,ul,ul),m,c(0,1,0,1), outer.ok = TRUE)$design
    X <- matrix(0,n,ncol(D)) ## full predict matrix
    nin <- sum(ind)
    if (nin>0) X[ind,] <- 
         splines::spline.des(object$knots,x[ind],m,rep(object$deriv,nin), outer.ok = TRUE)$design ## interior rows
    ## Now add rows for linear extrapolation (of smooth itself)...
    if (object$deriv<2) { ## under linear extrapolation higher derivatives vanish.
      ind <- x < ll 
      if (sum(ind)>0) X[ind,] <- if (object$deriv==0) cbind(1,x[ind]-ll)%*%D[1:2,] else 
                                 matrix(D[2,],sum(ind),ncol(D),byrow=TRUE)
      ind <- x > ul
      if (sum(ind)>0) X[ind,] <- if (object$deriv==0) cbind(1,x[ind]-ul)%*%D[3:4,] else 
                                 matrix(D[4,],sum(ind),ncol(D),byrow=TRUE)
    }
  }
  if (object$mono==0) X else X %*% object$B
} ## Predict.matrix.pspline.smooth

