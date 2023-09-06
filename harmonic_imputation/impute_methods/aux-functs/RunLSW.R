#Load required R.matlab package

library(R.matlab)
library(locits)
library(lpacf)
library(MLmetrics)

lacv.fc.leftwp<-function(x, h=0, filter.number = 1, family = "DaubExPhase", binwidth=0, smooth.dev = var, AutoReflect = TRUE, lag.max = NULL, WPsmooth.type = "RM", ...)
  {
  # modification of the lacv function in the locits package
  # modified to do non-2^J 
  # left facing code is in the below but is commented out - thus this is still a central window calculation
  # modified to do automatic bandwidth selection
  # Also removed return of most of the list which we do not require and changed default to Haar
  
  # Added h parameter for forecast horizon, uses Running Mean smoother to do the forecast on the spectral level
  # Defaults to 0 which is no forecast, i.e. same as lacv.leftwp function
  
  # this function is not exported so checks on inputs are done by parent functions.
  
  # rest of function
  TT <- length(x)
  filter<-wavethresh::filter.select(filter.number,family)
  
  Jp <-ceiling(logb(TT,2))
  add<-2^Jp-TT
  
  Nh<-length(filter$H!=0) # NvSK says Nh is # of non-zero elements in the filter.
  
  # The support of the discrete wavelets is something like this:
  # Lj<-(2^j -1)*(Nh-1)+1
  #
  # For Haar wavelets, this is 2^j, not really surprising.  The modified periodogram in Fryzlewicz 2003
  # extends the left hand wavelet value to the left hand edge, for locations 0,...,Lj-2.
  #
  
  xa<-c(rep(0,times=add),x)
  lxa<-length(xa)#should be 2^(J+1) if TT not equal to  2^J
  
  dsname = deparse(substitute(x))
  if (binwidth==0){
    binwidth=locits::AutoBestBW(x=xa,filter.number = filter.number, family = family, 
                                smooth.dev = smooth.dev, AutoReflect = AutoReflect,...)
  }
  if(binwidth>=TT){
    # i.e. binwidth is larger than data length
    binwidth=locits::AutoBestBW(x=x[(TT-2^{Jp-1}+1):TT],filter.number = filter.number, family = family, 
                                smooth.dev = smooth.dev, AutoReflect = AutoReflect,...)
  }
  
  #if(binwidth<=h){stop(paste('Automatic Bandwidth Selection not compatible with h-step ahead forecast.  Choose a forecast horizon less than',binwidth,'.'))}
  if(binwidth<=h){binwidth = h + 1}
  Pwd <- locits::ewspec3(x = xa, filter.number = filter.number, family = family, 
                         smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = WPsmooth.type, 
                         binwidth = binwidth, ...)$S
  
  P<-matrix(Pwd$D,ncol=lxa,byrow=TRUE)
  
  K<-(2^(1:Jp)-1)*(Nh-1)+1
  Jt<-max(which(K<=TT))
  
  # modify left periodogram edge:
  #   Smat <- matrix(NA, nrow = TT+h, ncol = Jt) # extra rows for the forecast points to be added below
  #   for (j in 1:Jt){
  #     Lj<-K[j]
  #     Smat[Lj:TT,j]<-P[j,(1+add):(lxa-Lj+1)]
  #     Smat[1:(Lj-1),j]<-rep(Smat[Lj,j],times=Lj-1)
  #   }
  # keep central periodogram but get rid of the artefacts introduced by the non 2^J and put some blank spaces ready for the extrapolation
  Smat <- matrix(NA, nrow = TT+h, ncol = Jt) # extra rows for the forecast points to be added below
  for (j in 1:Jt){
    Lj<-K[j]
    Smat[1:TT,j]<-P[j,(1+add):lxa]
  }
  
  # forecast the spectrum at the end of the data
  if(h>0){
    for(htmp in 1:h){
      Smat[TT+htmp,]=apply(Smat[(TT-binwidth+htmp):(TT+htmp-1),],2,mean) 
      # does an extrapolated average that uses previous extrapolations
    }    
  }
  
  #  S <- EWS$S
  #  J <- Pwd$nlevels
  Psi <- PsiJmat(-Jt, filter.number = filter.number, family = family)
  nc <- ncol(Psi)
  L <- (nc - 1)/2
  dimnames(Psi) <- list(NULL, c(-L:0, 1:L))
  if (is.null(lag.max)) 
    the.lacv <- Smat %*% Psi[, (L + 1):ncol(Psi)]
  else {
    if (L + 1 + lag.max > ncol(Psi)) {
      warning(paste("lag.max too high. Have reset it to ", 
                    ncol(Psi) - L - 1, ". Higher lags are zero"))
      lag.max <- ncol(Psi) - L - 1
    }
    the.lacv <- Smat %*% Psi[, (L + 1):(L + 1 + lag.max)]
  }
  the.lacor <- sweep(the.lacv, 1, the.lacv[, 1], FUN = "/")
  out=list(lacv = the.lacv, lacr = the.lacor,S=Smat,binwidth=binwidth,filter.number=filter.number,family=family)
  class(out)='lacv'
  return(out)
}
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}

sig.sq <-function(predeq, minimiser = solve(predeq$B, predeq$RHS))
{
  
  # input checks
  if(!is(predeq,'predeq')){stop('predeq must be of class "predeq" as output from predeq.est.')}
  if(!is.numeric(minimiser)){stop('minimiser must be numeric.')}
  
  # rest of function
  s <- length(predeq$RHS)
  B <- matrix(0, s + 1, s + 1)
  B[1:s, 1:s] <- predeq$B
  B[1:s, s + 1] <- predeq$RHS
  B[s + 1, 1:s] <- predeq$RHS
  B[s + 1, s + 1] <- predeq$extra.var
  b <- c(minimiser, -1)
  sqerr <- sum(b * (B %*% b))
  return(sqerr)
}

reg.xyr=function(pr)
{
  # Function to implement regularization from Xie, Yu, Ranneby 2007
  # pr is output from pred.est.gw
  
  # no input checks as this function is not exported and check are done in parent function.
  
  fn=function(lambda){
    is.positive.definite=function(lambda){
      val = try(chol(pr$B-diag(lambda,nrow=nrow(pr$B),ncol=ncol(pr$B))), silent = TRUE)
      #      if (class(val) == "try-error") {
      # MAN edited:
      if ("try-error" %in% class(val)) {
        return(FALSE)
      }
      else {
        return(TRUE)
      }      
    }
    check=is.positive.definite(lambda) # checks satisfies (18) from xyr paper
    if(check){
      inv.bracket=solve(pr$B-diag(lambda,nrow=nrow(pr$B),ncol=ncol(pr$B)))
      return((pr$RHS%*%inv.bracket%*%inv.bracket%*%pr$RHS -1)^2)
    }
    else{return(Inf)}
  }
  
  if(is.matrix(pr$B)){
    op.lambda=suppressWarnings(optimize(fn,interval=c(-1,1))$minimum)
    pr$B=pr$B-diag(op.lambda,nrow=nrow(pr$B),ncol=ncol(pr$B))    
  }# if not a matrix then there is no problems with regularity and thus nothing to be done.
  return(pr)
}

lpacf.leftwp.end=function (x, binwidth, lag.max=NULL, filter.number=1, family="DaubExPhase", smooth.dev=var,
                           AutoReflect=TRUE, tol=0.1, maxits=5, ABBverbose=0, lapplyfn=lapply) 
{
  
  # this is the lpacf function from GPN 
  # this is still central facing, code for modified to be left facing is commented out
  # also modified to do non 2^J lengths
  # it only returns the lpacf for the final time point
  # thus alot of the function has been removed as we only need to do it for one timepoint
  
  # below has been taken from lacv.fc.leftwp
  # rest of function
  TT <- length(x)
  filter<-wavethresh::filter.select(filter.number,family)
  
  Jp <-ceiling(logb(TT,2))
  add<-2^Jp-TT
  
  Nh<-length(filter$H!=0) # NvSK says Nh is # of non-zero elements in the filter.
  
  # The support of the discrete wavelets is something like this:
  # Lj<-(2^j -1)*(Nh-1)+1
  #
  # For Haar wavelets, this is 2^j, not really surprising.  The modified periodogram in Fryzlewicz 2003
  # extends the left hand wavelet value to the left hand edge, for locations 0,...,Lj-2.
  #
  
  xa<-c(rep(0,times=add),x)
  lxa<-length(xa)#should be 2^(J+1) if TT not equal to  2^J
  
  dsname = deparse(substitute(x))
  
  if (missing(binwidth) || binwidth == 0) 
    binwidth <- locits::AutoBestBW(x = xa, filter.number = filter.number, 
                                   family = family, smooth.dev = smooth.dev, AutoReflect = AutoReflect, 
                                   tol = tol, maxits = maxits, plot.it = FALSE, 
                                   
                                   verbose = ABBverbose,ReturnAll=FALSE)
  
  if(binwidth>=TT){
    # i.e. binwidth is larger than data length
    binwidth=locits::AutoBestBW(x=x[(TT-2^{Jp-1}+1):TT],filter.number = filter.number, family = family, 
                                smooth.dev = smooth.dev, AutoReflect = AutoReflect,tol = tol, maxits = maxits,
                                plot.it = FALSE, verbose = ABBverbose,ReturnAll=FALSE)
  }
  
  #if(binwidth<=lag.max){stop(paste('Automatic Bandwidth Selection not compatible with lag.max-step ahead forecast.  Choose a forecast horizon less than',binwidth,'.'))}
  if(binwidth<=lag.max){binwidth = lag.max + 1}
  # above has been taken from lacv.fc.leftwp  
  
  n <- length(x)
  
  if (is.null(lag.max)) 
    lag.max <- floor(10 * (log10(n)))
  
  start <- (n - round(binwidth/2)) # (n-binwidth+1) # commented is for left facing
  end <- n
  
  out=list(lpacf=pacf(x[seq(from=start,to=end)],lag.max=lag.max,plot=FALSE)$acf[,,1],binwidth=binwidth)
  class(out)='lpacf'
  return(out)
  
}

forecastlpacf <-function (x, h = 1, binwidth = 0, regularize = TRUE, lag.max = max(10, 2 * 
                                                         h), forecast.type = NULL, ...) 
{
    xdata = x
    
    if (sum(is.na(x)) > 0) {
      stop("Time series must not contain NA values.")
    }
    if ((h <= 0) | !is.wholenumber(h)) {
      stop("Forecast horizon, h, must be a positive integer.")
    }
    if (!is.logical(regularize)) {
      stop("regularize must be logical.")
    }
    if (lag.max <= h) {
      stop(paste("lag.max must be atleast:", h + 1))
    }
    if (lag.max > length(x)) {
      stop("lag.max must be within the range of the data")
    }
    if (is.null(forecast.type)) {
      stop("forecast.type must be either recursive, fixed or extend.")
    }
    else if (!(forecast.type %in% c("recursive", "fixed", "extend"))) {
      stop("forecast.type must be either recursive, fixed or extend.")
    }
    
    len = length(x)
    
    sp.lacv = lacv.fc.leftwp(x, h, binwidth = binwidth, ...)
    lpacf = lpacf.leftwp.end(x,binwidth, lag.max = lag.max, ...)$lpacf
    ci = qnorm((1 + 0.95)/2)/sqrt(sp.lacv$binwidth)
    p = ifelse(!length(which(((lpacf > ci) | (lpacf < -ci)) == 
                               FALSE)), lag.max, which(((lpacf > ci) | (lpacf < -ci)) == 
                                                         FALSE)[1] - 1)
    h = 1:h
    if (p == 0) {
      p = 1
    }
    if (p == lag.max) {
      warning("It is advised that you increase your maximum lag to aid your forecast estimation")
    }
    pmean = NULL
    std.err = NULL
    if (forecast.type != "recursive") {
      p = rep(p, length(h))
      if (any(p <= max(h))) {
        warning("p is less than the maximum forecast horizon, using forecast.type=\"extend\" or \"fixed\" matters!")
        if (forecast.type == "extend") {
          p[p < h] = h[p < h]
        }
      }
      B <- matrix(0, max(p), max(p))
      for (i in 1:max(p)) {
        B[i, i:max(p)] = sp.lacv$lacv[len - max(p) + i, 1:length(i:max(p))]
      }
      tmp = diag(B)
      B = t(B) + B
      diag(B) = tmp
      for (i in 1:length(p)) {
        if (p[i] == 1) {
          RHS = 0
        }
        else {
          RHS <- c(0, sp.lacv$lacv[len + h[i], 2:p[i]])
        }
        extra.var = sp.lacv$lacv[len + h[i], 1]
        pr = list(B = B[(max(p) - p[i] + 1):max(p), (max(p) - 
                                                       p[i] + 1):max(p)], RHS = RHS, extra.var = extra.var)
        class(pr) = "predeq"
        if (regularize == TRUE) {
          pr = reg.xyr(pr)
        }
        b <- solve(pr$B, pr$RHS)
        pmean[i] <- sum(x[(len - p[i] + 1):len] * b)
        std.err[i] <- sqrt(max(sig.sq(pr, b), 1e-10))
      }
    }
    else {
      xext = x
      for (H in h) {
        B <- matrix(0, p, p)
        for (i in 1:p[1]) {
          B[i, i:p] = sp.lacv$lacv[len - p + i, 1:length(i:p)]
        }
        tmp = diag(B)
        B = t(B) + B
        diag(B) = tmp
        if (p == 1) {
          RHS = 0
        }
        else {
          RHS <- c(0, sp.lacv$lacv[len + H, 2:p])
        }
        extra.var = sp.lacv$lacv[len + H, 1]
        pr = list(B = B, RHS = RHS, extra.var = extra.var)
        class(pr) = "predeq"
        if (regularize == TRUE) {
          pr = reg.xyr(pr)
        }
        b <- solve(pr$B, pr$RHS)
        pmean[H] = sum(xext[(length(xext) - p + 1):length(xext)] * 
                         b)
        xext = c(xext, pmean[H])
        std.err[H] <- sqrt(max(sig.sq(pr, b), 1e-10))
      }
    }
    out = list(mean = pmean, std.err = std.err, lpacf = lpacf, 
               ci = ci, binwidth = sp.lacv$binwidth, p = p, x = xdata, d=0)
    class(out) = "forecastlpacf"
    return(out)
  }

RunLSW <- function(name) {
  #Load data saved from Matlab
  #source('temp/MissingSamples.R')
  
  source(name)
  
  #Run one sample t test and save results
  
  sp_ts = ts(si,start = c(1))
  forecast <- forecastlpacf(sp_ts,h=Li,lag.max=min(Li+1,Ni), filter.number=8,family="DaubExPhase",
                            forecast.type='extend',regularize=TRUE)
  
  #Save results to Matlab data file
  
  imputed_values = forecast$mean
  
  #writeMat('temp/Imputevalues.mat',imp=imputed_values)
  writeMat(paste(name,'.mat',sep=''),imp=imputed_values)
}


args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
RunLSW(name)

###########################################################################################
# first generate some non-stationary data we want to forecast
#set.seed(1)
#x=tvar2sim()
#predict 1-step ahead using Daubechies wavelets with 2 vanishing moments, although 
#other choices for the wavelet family and filter are possible (including Haar)
#pred<-forecastlpacf(x,h=1,filter.number=8,family="DaubExPhase",forecast.type='recursive')

#pred$mean gives the predicted value, while pred$std.err gives the prediction erro
###########################################################################################