##  ========================================================================  ##
##  Miguel de Carvalho                                                        ##
##  Copyright (C) 2018                                                        ##
##  ------------------------------------------------------------------------  ##
##  This program is free software; you can redistribute it and/or modify      ##
##  it under the terms of the GNU General Public License as published by      ##
##  the Free Software Foundation; either version 2 of the License, or         ##
##  (at your option) any later version.                                       ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful,           ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of            ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             ##
##  GNU General Public License for more details.                              ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License         ##
##  along with this program; if not, a copy is available at                   ##
##  http://www.r-project.org/Licenses/                                        ##
##  ========================================================================  ##

cdensity <- function(XY, tau= 0.95,raw=TRUE, ...)
    UseMethod("cdensity")

cdensity.default <- function(XY, tau=0.95,raw=TRUE, ...) {
    if((is.matrix(XY)|is.data.frame(XY))==FALSE)
        stop('data needs to be a matrix or data frame')
    dim<-ncol(XY)
    if(dim==2){
        y <- XY[,2]
        raw<-TRUE
    }
    else{
        ## Convert data to common margins, if raw == TRUE
        if(raw == TRUE) {
            n <- dim(XY)[1]
            FX <- ecdf(XY[, 2])
            FY <- ecdf(XY[, 3])
            x1 <- -1 / log(n / (n + 1) * FX(XY[, 2]))
            x2 <- -1 / log(n / (n + 1) * FY(XY[, 3]))
        }
        if(raw == FALSE) {
            x1 <- XY[, 2]
            x2 <- XY[, 3]
        }
        y<- apply(cbind(x1,x2),1,min)}
    threshold <-quantile(y, tau)
    ## Basic input validation
    if (as.numeric(threshold) >= max(y))
        stop('threshold cannot exceed max(y)')

    ## Initialize variables
    T <- length(y)
    w <- which(y > threshold) / T
    k <- length(w) 
    c <- density(w, n = T, from = 1 / T, to = 1, ...)
    
    ## Organize and return outputs    
    outputs <- list(c = c, w = w, k = k, T = T, XY = XY)
    class(outputs) <- "cdensity"
    return(outputs)
}

    plot.cdensity <- function(x, rugrep = TRUE,
                          original = TRUE, main = "", ...) {
    if(original == TRUE) {
        plot(x$XY[, 1], x$c$y, xlab = "Time", ylab = "Scedasis Density", 
             main = "", type = "S", ...)
        if(rugrep == TRUE)
            rug(x$XY[x$w * x$T, 1])
    }
    if(original == FALSE) {
        par(pty = "s")
        plot(x$c, xlab = "w", ylab = "Scedasis Density", 
             main = "", ...)
        if (rugrep == TRUE)
          rug(x$w)
    }
}

    
