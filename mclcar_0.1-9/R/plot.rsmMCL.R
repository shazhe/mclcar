### plot the rsm results
plot.rsmMCL <- function(x, family, trace.all = TRUE, plain=TRUE, ...){
    plot.exact <- family == "gauss"
    ## extract all the elemnts in the list object of the rsm results
    ##for(i in 1:length(x)){
    ##    tempobj=x[[i]]
    ##    eval(parse(text=paste(names(x)[[i]],"= tempobj")))
###}
	RSM.1 <- x$RSM.1
	Psi <- x$Psi
	psi0 <- x$psi0
	FVbox <- x$FVbox
	DAbox <- x$DAbox
	ESAs <- x$ESAs
	Exacts <- x$Exacts
	N.iter <- x$N.iter
	RSM.2 <- x$RSM.2
      N.SA <- length(ESAs)

    if(trace.all){
        Psis <- do.call(rbind, Psi)
        Psis <- rbind(psi0, Psis)
        for (i in 1:N.iter) {
            if (i <= N.SA){
                SAn <- !is.null(ESAs[[i]])
            }else{
                SAn <- FALSE
            }
            if (SAn){
                par(mfrow = c(1,3))
                if (is.null(RSM.2[[i]])){
                    RSM.x <- RSM.1[[i]]
                }else{
                    RSM.x <-RSM.2[[i]]
                }
                plot.RSM(RSM.x, psi = Psis[i, ], bounds = list(x1=c(-1.1,2), x2 = c(-1.5,1.5)),
                         fvbox = FVbox[[i]], DAbox = DAbox[[i]], A.path = "steepest",
                         distance = ESAs[[i]]$ds, exact = Exacts[[i]], exact.vals = plot.exact)
                title(main = list(paste("Iteration ", i)))
                SA <- ESAs[[i]][[1]]
                plot(mc.lr ~ dist, data = SA$SA.fit, xlab = "distance", ylab = "Monte Carlo likelihood")
                lines(SA$pred.val ~ SA$pred.dist)
                plot(log(mc.var) ~ dist, data = SA$SA.fit, xlab = "distance",
                     ylab = "Monte Carlo variance")
                lines(SA$pred.var ~ SA$pred.dist)
            }else{
                par(mfrow = c(1,1))
                plot.RSM(RSM.2[[i]], psi = Psis[i, ], bounds = list(x1=c(-1.1,2), x2 = c(-1.5,1.5)),
                         fvbox = FVbox[[i]], DAbox = DAbox[[i]], A.path = "canonical",
                         exact = Exacts[[i]], exact.vals = plot.exact)
                title(main = list(paste("Iteration ", i), cex = 1.1))
            }
        }
        }else{
            ## plot the RSM of the final iteration
            if(N.iter <=N.SA){
                SAn <- is.null(ESAs[[N.iter]])
            }else{
                SAn <- FALSE
            }
            if (SAn){
                par(mfrow = c(1,3))
                if (is.null(RSM.2[[N.iter]])){
                    RSM.x <- RSM.1[[N.iter]]
                }else{
                    RSM.x <-RSM.2[[N.iter]]
                }
                plot.RSM(RSM.x, psi = Psi[[N.iter-1] ], bounds = list(x1=c(-1.1,2), x2 = c(-1.5,1.5)),
                         fvbox = FVbox[[N.iter]], DAbox = DAbox[[N.iter]], A.path = "steepest",
                         distance = ESAs[[N.iter]]$ds, exact = Exacts[[N.iter]], exact.vals = plot.exact)
                title(main = "Final Iteration")
                SA <- ESAs[[N.iter]][[1]]
                plot(mc.lr ~ dist, data = SA$SA.fit, xlab = "distance", ylab = "Monte Carlo likelihood")
                lines(SA$pred.val ~ SA$pred.dist)
                plot(log(mc.var) ~ dist, data = SA$SA.fit, xlab = "distance",
                     ylab = "Monte Carlo variance")
                lines(SA$pred.var ~ SA$pred.dist)
            }
            else{
                par(mfrow = c(1,1))
                plot.RSM(RSM.2[[N.iter]], psi = Psi[[N.iter-1]], bounds = list(x1=c(-1.1,2), x2 = c(-1.5,1.5)),
                         fvbox = FVbox[[N.iter]], DAbox = DAbox[[N.iter]], A.path = "canonical",
                         exact = Exacts[[N.iter]], exact.vals = plot.exact)
                title(main = "Final Iteration")
            }

        }
}
