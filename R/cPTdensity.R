cPTdensity <- function(XY, tau = 0.95, raw = TRUE, structure = "min", prior, mcmc, state, status,
                       data=sys.frame(sys.parent()),na.action=na.fail)
  UseMethod("cPTdensity")

cPTdensity.default <- function(XY, tau = 0.95, raw = TRUE, structure = "min", prior, mcmc, state,status,
                               data=sys.frame(sys.parent()),na.action=na.fail) {
  if((is.matrix(XY)|is.data.frame(XY))==FALSE)
    stop('data needs to be a matrix or data frame')
  dim<-ncol(XY)
  if(dim==2){
    y<-XY[,2]
  }
  else{
      ## Convert data to common margins, if raw == TRUE
  if(raw == TRUE) {
    n <- dim(XY)[1]
    FZ<-apply(XY[,-1], 2, function(c) -1/log(length(c)/(length(c) + 1)*ecdf(c)(c)))
                  }
  if(raw == FALSE) {
    FZ<-XY[,-1]
                  }
  if(structure == "min"){
      y<- apply(FZ,1,min)
  } 
  if(structure == "max"){
      y<- apply(FZ,1,max)
  } 
  if(structure == "sum"){
      y<- apply(FZ,1,sum)
  } 
                               }
    threshold <-quantile(y, tau)
 
   ## Basic input validation
  if (as.numeric(threshold) >= max(y))
    stop('threshold cannot exceed max(y)')
  
  ## Initialize variables
  T <- length(y)
  w <- which(y > threshold) / T
  k <- length(w) 
  
  ##########CODE FOR PTDENSITY
  ############################
  PTdensitybeta <- function(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
    UseMethod("PTdensitybeta")
  
  PTdensitybeta.default <- function(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
  {
    #########################################################################################
    # call parameters
    #########################################################################################
    cl <- match.call()
    y <- na.action(as.matrix(y))	
    
    #########################################################################################
    # data structure
    #########################################################################################
    nrec <- nrow(y)
    nvar <- ncol(y)
    
    if(nvar>1)
    {
      left <- rep(0,2)
      right <- rep(0,2)
      
      left[1] <- min(y[,1])-0.5*sqrt(var(y[,1]))
      right[1] <- max(y[,1])+0.5*sqrt(var(y[,1]))
      
      left[2] <- min(y[,2])-0.5*sqrt(var(y[,2]))
      right[2] <- max(y[,2])+0.5*sqrt(var(y[,2]))
    } 
    else
    {
      left <- 0
      right <-1
    }    
    
    #########################################################################################
    # prior information
    #########################################################################################
    
    jfr <- c(0,0)
    if(nvar==1)
    {
      tau <- c(-1,1)
      m0 <- 0
      S0 <- 1
    } 
    else
    {
      nu0 <- -1
      tinv <- diag(1,nvar)
      m0 <- rep(0,nvar)
      S0 <- diag(1,nvar)
    }
    
    if(is.null(prior$al))
    {
      arand <- 1
      if(is.null(prior$m0))
      {
        jfr[1] <- 1
      }
      else
      {
        m0 <- prior$m0
        S0 <- prior$S0
      }
    }
    else
    {  
      arand <- 0
      al <- prior$al
    }
    
    if(is.null(prior$be))
    {
      brand <- 1
      if(nvar==1)
      { 
        if(is.null(prior$tau1))
        {
          jfr[2] <- 1
        }
        else
        {
          tau <- c(prior$tau1,prior$tau2)
        }
      }
      else
      {
        if(is.null(prior$nu0))
        {
          jfr[2] <- 1
        }
        else
        {
          nu0 <- prior$nu0
          tinv <- prior$tinv
        }   
      }
    }
    else
    {  
      brand <- 0
      be <- prior$be
    }
    
    if(is.null(prior$a0))
    {
      ca <- -1
      cb <- -1 
      cpar <- prior$alpha
      crand <- 0
    }
    else
    {
      ca <- prior$a0
      cb <- prior$b0
      cpar <- 10
      crand <- 1
    }
    ab <- c(ca,cb)
    
    #########################################################################################
    # mcmc specification
    #########################################################################################
    mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
    nsave <- mcmc$nsave
    
    if(is.null(mcmc$tune1))         
    {
      tune1 <- 1.1
    }
    else
    {
      tune1 <- mcmc$tune1
    }
    
    if(is.null(mcmc$tune2))         
    {
      tune2 <- 1.1
    }
    else
    {
      tune2 <- mcmc$tune2
    }
    
    if(is.null(mcmc$tune3))         
    {
      tune3 <- 1.1
    }
    else
    {
      tune3 <- mcmc$tune3
    }
    
    
    #########################################################################################
    # output
    #########################################################################################
    acrate <- rep(0,3)
    
    if(nvar==1)
    {
      f <- rep(0,ngrid)
    }
    else
    {
      ngrid <- as.integer(sqrt(ngrid))
      f <- matrix(0,nrow=ngrid,ncol=ngrid)
      fun1 <- rep(0,ngrid)
      fun2 <- rep(0,ngrid)
    }
    
    thetasave <- matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+1)
    randsave <- matrix(0,nrow=nsave,ncol=nvar)
    estimasave<-matrix(0,nrow=nsave,ncol=ngrid)
    
    #########################################################################################
    # parameters depending on status
    #########################################################################################
    
    if(status==TRUE)
    {
      if(nvar==1)
      {
        if(arand==1)
        {
          ymed <- mean(y)
          al<-ymed*(((ymed*(1-ymed))/var(y))-1)
        }
        if(brand==1)
        {
          ymed <- mean(y)
          be <- (1-ymed)*(((ymed*(1-ymed))/var(y))-1)
        }
        
      }}
    
    if(status==FALSE)
    {
      cpar <- state$alpha
      al <- state$al 
      be <- state$be
    }    
    
    #########################################################################################
    # working space
    #########################################################################################
    
    acrate <- rep(0,3)
    cpo <- rep(0,nrec)
    if(nvar==1)
    {
      if(is.null(grid))
      {
        grid <- seq(left,right,length=ngrid)
      }
      else
      {
        grid <- as.matrix(grid)
        ngrid <- nrow(grid)
        grid <- as.vector(grid)
      }
    }        
    seed <- c(sample(1:29000,1),sample(1:29000,1))
    
    #########################################################################################
    # calling the fortran code
    #########################################################################################
    
    if(nvar==1)
    {
      #if(is.null(prior$M))
      # {
      #  whicho <- rep(0,nrec)
      #  whichn <- rep(0,nrec)
      #  foo <- .Fortran("ptdensitybetau",
      #  ngrid      =as.integer(ngrid),
      #   nrec       =as.integer(nrec),
      #    y          =as.double(y),
      #   ab         =as.double(ab),
      #   arand     =as.integer(arand),
      #   brand  =as.integer(brand),
      #   jfr        =as.integer(jfr),
      #   m0         =as.double(m0),
      #   s0         =as.double(S0),  
      #   tau        =as.double(tau),
      #   mcmcvec    =as.integer(mcmcvec),
      #   nsave      =as.integer(nsave),
      #   tune1      =as.double(tune1),
      #   tune2      =as.double(tune2),
      #  tune3      =as.double(tune3),
      #   acrate     =as.double(acrate),
      #  f          =as.double(f),
      #   thetasave  =as.double(thetasave),		
      #   cpo        =as.double(cpo),		
      #  cpar       =as.double(cpar),		
      #   a         =as.double(a),		
      #   b      =as.double(b),		
      #   grid       =as.double(grid),		
      #  seed       =as.integer(seed),
      #  whicho     =as.integer(whicho),
      #   whichn     =as.integer(whichn),
      #   PACKAGE    ="extremis")
      # }	
      # else
      # {
      nlevel <- prior$M
      ninter <- 2**nlevel
      assign <- matrix(0,nrow=nrec,ncol=nlevel)
      accums <- matrix(0,nrow=nlevel,ncol=ninter)
      counter <- matrix(0,nrow=nlevel,ncol=ninter)
      endp <- rep(0,ninter-1)
      intpn <- rep(0,nrec)
      intpo <- rep(0,nrec)
      prob <- rep(0,ninter)
      rvecs <- matrix(0,nrow=nlevel,ncol=ninter)
      
      #foo <- .Fortran("ptdensityubp",
      foo <- .Fortran("ptdensitybetaupmh",
                      ngrid      =as.integer(ngrid),
                      nrec       =as.integer(nrec),
                      y          =as.double(y),
                      ab         =as.double(ab),
                      arand     =as.integer(arand),
                      brand  =as.integer(brand),
                      jfr        =as.integer(jfr),
                      m0         =as.double(m0),
                      s0         =as.double(S0),  
                      tau        =as.double(tau),
                      nlevel     =as.integer(nlevel),
                      ninter     =as.integer(ninter),
                      mcmcvec    =as.integer(mcmcvec),
                      nsave      =as.integer(nsave),
                      tune1      =as.double(tune1),
                      tune2      =as.double(tune2),
                      tune3      =as.double(tune3),
                      acrate     =as.double(acrate),
                      f          =as.double(f),
                      thetasave  =as.double(thetasave),		
                      estimasave  =as.double(estimasave),
                      cpo        =as.double(cpo),		
                      cpar       =as.double(cpar),		
                      al         =as.double(al),		
                      be         =as.double(be),		
                      grid       =as.double(grid),		
                      intpn     =as.integer(intpn),		
                      intpo     =as.integer(intpo),		
                      accums    =as.double(accums),
                      assign    =as.integer(assign),
                      counter   =as.integer(counter),
                      endp      =as.double(endp),
                      prob      =as.double(prob),
                      rvecs     =as.double(rvecs),
                      seed      =as.integer(seed),
                      PACKAGE    ="extremis")
    }
    # }   
    
    
    #########################################################################################
    # save state
    #########################################################################################
    model.name <- "Bayesian Density Estimation Using MPT"		
    
    varnames<-colnames(y)
    if(is.null(varnames))
    {
      varnames<-all.vars(cl)[1:nvar]
    }
    
    state <- list(alpha=foo$cpar,
                  al=foo$al,
                  be=matrix(foo$be,nrow=nvar,ncol=nvar)
    )
    
    thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=(nvar+nvar*(nvar+1)/2+1))
    estimasave<-matrix(foo$estimasave,nrow=nsave,ncol=ngrid)
   
    if(nvar>1)
    {
      randsave <- matrix(foo$randsave,nrow=nsave,ncol=nvar)
      colnames(randsave) <- varnames
    }   
    
    coeff<-apply(thetasave,2,mean) 
    
    pnames1<-NULL
    for(i in 1:nvar)
    {
      pnames1<-c(pnames1,paste("al",varnames[i],sep=":"))
    }
    pnames2<-NULL
    for(i in 1:nvar)
    {
      for(j in i:nvar)
      {
        if(i==j)
        {
          tmp <- varnames[i]
        }
        else
        {
          tmp <- paste(varnames[i],varnames[j],sep="-")
        }   
        pnames2 <- c(pnames2,paste("be",tmp,sep=":"))
      }	
    }
    
    names(coeff) <- c(pnames1,pnames2,"alpha")
    colnames(thetasave) <- c(pnames1,pnames2,"alpha")
    save.state <- list(thetasave=thetasave,randsave=randsave)
    
    if(crand==0)
    {
      acrate<-foo$acrate[1:2]
    }
    else
    {
      acrate<-foo$acrate
    }
    
    x1<-NULL
    x2<-NULL
    dens<-NULL
    
    if(nvar==1)
    {
      x1<-foo$grid
      dens<-foo$f
      f<-foo$f
      grid1<-foo$grid
      grid2<-NULL
      fun1<-foo$f
      fun2<-NULL
    }
    else
    {
      x1 <- grid1
      x2 <- grid2
      dens <- matrix(foo$f,nrow=ngrid,ncol=ngrid)
      f <- matrix(foo$f,nrow=ngrid,ncol=ngrid)
      
      dist1 <- grid2[2]-grid2[1] 
      dist2 <- grid1[2]-grid1[1] 
      fun1 <- (dist1/2)*(dens[,1]+dens[,ngrid]+2*apply(dens[,2:(ngrid-1)],1,sum))
      fun2 <- (dist2/2)*(dens[1,]+dens[ngrid,]+2*apply(dens[2:(ngrid-1),],2,sum))
    }   
    
    z<-list(call=cl,y=y,varnames=varnames,modelname=model.name,cpo=foo$cpo,
            prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
            nvar=nvar,crand=crand,coefficients=coeff,f=f,grid1=grid1,grid2=grid2,
            fun1=fun1,fun2=fun2,x1=x1,x2=x2,dens=dens,acrate=acrate,estimasave=estimasave)
    
    cat("\n\n")
    class(z)<-"PTdensitybeta"
    return(z)
  }
  
  
  ###                    
  ### Tools
  ###
  ### Copyright: Alejandro Jara Vallejos, 2006
  ### Last modification: 28-11-2006.
  ###
  
  
  "print.PTdensitybeta"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
  {
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
                  quote = FALSE) 
    
    cat("\nPosterior Inference of Parameters:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                  quote = FALSE)
    
    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Variables:",x$nvar,"\n")        
    cat("\n\n")
    invisible(x)
  }
  
  
  "summary.PTdensitybeta"<-function(object, hpd=TRUE, ...) 
  {
    stde<-function(x)
    {
      n<-length(x)
      return(sd(x)/sqrt(n))
    }
    
    hpdf<-function(x)
    {
      alpha<-0.05
      vec<-x
      n<-length(x)         
      alow<-rep(0,2)
      aupp<-rep(0,2)
      a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                  alow=as.double(alow),aupp=as.double(aupp),PACKAGE="extremis")
      return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
      alpha<-0.05
      vec<-x
      n<-length(x)         
      alow<-rep(0,2)
      aupp<-rep(0,2)
      a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                  alow=as.double(alow),aupp=as.double(aupp),PACKAGE="extremis")
      return(c(a$alow[2],a$aupp[2]))
    }
    
    thetasave<-object$save.state$thetasave
    nvar<-object$nvar
    
    if(object$crand==0)
    {
      dimen<-(nvar+nvar*(nvar+1)/2)
      mat<-matrix(thetasave[,1:dimen],ncol=2) 
      
    }
    else
    {
      dimen<-(nvar+nvar*(nvar+1)/2+1)
      mat<-thetasave[,1:dimen]
    }
    
    coef.p<-object$coefficients[1:dimen]
    coef.m <-apply(mat, 2, median)    
    coef.sd<-apply(mat, 2, sd)
    coef.se<-apply(mat, 2, stde)
    
    if(hpd){             
      limm<-apply(mat, 2, hpdf)
      coef.l<-limm[1,]
      coef.u<-limm[2,]
    }
    else
    {
      limm<-apply(mat, 2, pdf)
      coef.l<-limm[1,]
      coef.u<-limm[2,]
    }
    
    names(coef.m)<-names(object$coefficients[1:dimen])
    names(coef.sd)<-names(object$coefficients[1:dimen])
    names(coef.se)<-names(object$coefficients[1:dimen])
    names(coef.l)<-names(object$coefficients[1:dimen])
    names(coef.u)<-names(object$coefficients[1:dimen])
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    
    if(hpd)
    {
      dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                                                    "95%HPD-Low","95%HPD-Upp"))
    }
    else
    {
      dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                                                    "95%CI-Low","95%CI-Upp"))
    }
    
    ans <- c(object[c("call", "modelname")])
    
    ans$coefficients<-coef.table
    
    
    ### CPO
    ans$cpo<-object$cpo
    
    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec
    ans$nvar<-object$nvar
    
    class(ans) <- "summaryPTdensity"
    return(ans)
  }
  
  
  "print.summaryPTdensitybeta"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
  {
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
                  quote = FALSE) 
    
    if (length(x$coefficients)) {
      cat("\nBaseline parameters:\n")
      print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
    else cat("No coefficients\n")
    
    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Variables:",x$nvar,"\n")            
    cat("\n\n")
    invisible(x)
  }
  
  ##############################
  ##############################
  c <- PTdensitybeta(w, ngrid  = T,prior = prior,mcmc = mcmc,state,status = status,
                     data=sys.frame(sys.parent()),na.action=na.fail)
  ###############
  ## Organize and return outputs    
  outputs <- list(c = c, w = w, k = k, T = T, XY = XY)
  class(outputs) <- "cPTdensity"
  return(outputs)
}

plot.cPTdensity <- function(x, rugrep = TRUE,
                            original = TRUE, main = "", 
                            CI= FALSE, rquantiles = FALSE, ...) {
  if ( rquantiles == TRUE){
    salida<-cbind(x$c$x1,t(x$c$estimasave))
    delta<-salida[,1][2]-salida[,1][1]
    acumulado<-apply(salida[,-1],2,cumsum)
    mean<-apply(acumulado*delta,1,mean)
    intV<-t(apply(acumulado*delta, 1, quantile, probs = c(0.05, 0.95)))
    indexw<-numeric()
    for (j in 1:length(x$w)){
      indexw[j]<-sum(x$c$x1<=x$w[j])}
    varaux<-mean[indexw]
    min<-intV[,1][indexw]
    max<-intV[,2][indexw]
    #===
    varaux<-c(varaux,tail(mean,1))
    min<-c(min,tail(intV[,1],1))
    max<-c(max,tail(intV[,2],1))
    #===
    ylim <- range(varaux)
    n<-length(varaux)
    x <- qnorm(ppoints(n))[order(order(varaux))]
    y<-qnorm(varaux)
    ymin<-qnorm(min)[order(x)]
    ymax<-qnorm(max)[order(x)]
    ymax[which(is.na(ymax))]<-max(na.omit(ymax))
    ymin[which(is.na(ymin))]<-max(na.omit(ymin))
    xymin<- qnorm(ppoints(n))[order(order(ymin))]
    xymax<-qnorm(ppoints(n))[order(order(ymax))]
    interval<-data.frame(x=x[order(x)],mean=y[order(x)],
                         min=ymin[order(xymin)],max=ymax[order(xymax)])
    par(pty="s")
    plot(interval$x,interval$x,ylab="Sample quantiles",xlab="Theoretical quantiles",
         xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),type="l",...)
    polygon(c(rev(interval$x), interval$x), c(rev(interval$min), 
                                              interval$max), col = 'lightgrey', border = NA)
    lines(interval$x,interval$mean,type="S",lty=2)
  }
  if (rquantiles == FALSE){
  aux<-data.frame(x$XY[, 1],x$c$x1,x$c$dens)
  aux<-aux[which(aux[,2]>0.01 & aux[,2]<.99),]
  dim<-ncol(x$XY)
  if(dim==2){
    labs<-"Scedasis Density"
  }
  else{labs<-"Structure Scedasis Density"}
  if(CI==TRUE){
    ci<-data.frame(x$XY[,1],x$c$x1,t(x$c$estimasave))
    ci<-ci[which(ci[,2]>0.01 & ci[,2]<.99),]
    intV <- t(apply(ci[,-c(1,2)], 1, quantile, probs = c(0.1, 0.85)));
    min<-intV[,1];max<-intV[,2];
    yl<-max(min,max)
    if(original == TRUE) { 
      plot(aux[,1],aux[,3], xlab = "Time", ylab = labs, 
           main = "",
           type = "S", ylim=c(0,yl),...)
      polygon(c(rev(ci[,1]), ci[,1]), c(rev(min), max), col = 'lightgrey', border = NA)
      lines(aux[,1],aux[,3],type="S")
      if(rugrep == TRUE)
        rug(x$XY[x$w * x$T, 1])
    }
    if(original == FALSE) {
      par(pty = "s")
      plot(aux[,2],aux[,3], xlab = "w", ylab = labs, 
           main = "",ylim=c(0,yl), type = "S", ...)
      polygon(c(rev(ci[,2]), ci[,2]), c(rev(min), max), col = 'lightgrey', border = NA)
      lines(aux[,2],aux[,3],type="S")
      if (rugrep == TRUE)
        rug(x$w)
    }
  }
  else{
    if(original == TRUE) {
    plot(aux[,1],aux[,3], xlab = "Time", ylab = labs, 
         main = "",
         type = "S", ...)
    if(rugrep == TRUE)
      rug(x$XY[x$w * x$T, 1])
        }
  if(original == FALSE) {
    par(pty = "s")
    plot(aux[,2],aux[,3], xlab = "w", ylab = labs, 
         main = "",
         type = "S", ...)
    if (rugrep == TRUE)
      rug(x$w)
  }}}
  
 }




