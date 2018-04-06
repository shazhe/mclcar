#### Sub functions for RSM
#### Estimate the suitable range of the design area
range.tt <- function(psi, N, s, rho.r, k = 10, torus = TRUE, torus.a = 0) {
    ## n -- no obs
    ## s -- no MC samples
    ## k -- furthurest level set
    lambda1 <- 1/rho.r[1]
    lambdan <- 1/rho.r[2]
    rho.psi <- psi[1]
    ##sigma.psi <- psi[2]

    r.ts <- sqrt( 1 - (k^2*s/4 + 1)^(-2/N))
    if(rho.psi > (rho.r[1] + rho.r[2])/2){
        alpha <-  ifelse(torus, lambdan*(1 - rho.psi*lambdan)/(2-rho.psi*lambdan)/(lambdan-lambda1),0.5)
        b <-  sqrt(1 - (k^2*s/4 + 1)^(-2/(alpha*N)))
        r1.tr <- -b*(1/lambdan - rho.psi)/(1-b)
        r2.tr <- b*(1/lambdan - rho.psi)/(1+b)
    }else{
        alpha <- ifelse(torus.a == 0, -lambda1*(1 - rho.psi*lambda1)/(2-rho.psi*lambda1)/(lambdan-lambda1), torus.a)
        b <-  sqrt(1 - (k^2*s/4 + 1)^(-2/(alpha*N)))
        r1.tr <- b*(1/lambda1 - rho.psi)/(1+b)
        r2.tr <- -b*(1/lambda1 - rho.psi)/(1-b)
    }

    ## control the boundary
    ds1 <- psi[1] - rho.r[1]
    ds2 <- rho.r[2] - psi[1]
    if(abs(ds1/rho.r[1]) < 0.1){
        r1.tr <-  -(ds1 - 1e-6)
        r2.tr <- ds1*2
    }else if(abs(ds2/rho.r[2]) < 0.1){
         r2.tr <- ds2 - 1e-6
         r1.tr <- -ds2*2
     }

    return(c(max(rho.r[1] + 1e-6,r1.tr), min(rho.r[2]-1e-6,r2.tr), r.ts))
}


#### CAR models(rho + sigma):
#### Generate all the design points for evaluation
############################################################
RSM.design.CAR <- function(design.tt, psi, n01=4, n02=2, rho.rr=c(-0.25, 0.25)){
    rho.psi <- psi[1]
    sigma.psi <- psi[2]
    sigma.r <- sigma.psi*(-exp(-design.tt[3]) + exp(design.tt[3]))/2
    rho.r <- (design.tt[2] - design.tt[1])/2
    sigma.c <- sigma.psi*(exp(-design.tt[3]) + exp(design.tt[3]))/2
    rho.c <- rho.psi + (design.tt[2] + design.tt[1])/2

    f1 <- as.formula(paste0("x1 ~ (rho - ", eval(rho.c), ")/", eval(rho.r)))
    f2 <- as.formula(paste0("x2 ~ (sigma - ", eval(sigma.c), ")/", eval(sigma.r)))
    expt1 <- cube(~ x1 + x2, n0=n01, coding = c(f1, f2))
    expt2 <- star(expt1, n0 = n02, alpha = "orthogonal")
    expt2.dd <- any(decode.data(expt2)$rho > rho.rr[2] | decode.data(expt2)$rho < rho.rr[1])
    if(expt2.dd){
        expt2 <- star(expt1, n0 = n02, alpha = "face")
    }
    return(list(expt1, expt2))
}

#### Direct CAR models (rho + sigma):
#### Evaluate at the design points
############################################################
Exp.dCAR <- function(exp.design, data, psi, n0, mc.data){
    if(is.null(mc.data)){
         stop("Please provide the Monte Carlo samples by specifying mc.data!")
    }
    n.samples <- ncol(mc.data$simY)
    nb.sim <- floor(n.samples/n0)

    mc.lr <- mc.var <- rep(0, 4+n0)
    exp.design$mc.lr <- mc.lr
    exp.design$mc.var <- mc.var
    beta0 <- as.numeric(psi[-c(1,2)])
    decode.design <- decode.data(exp.design)

    for (k in 1:nrow(exp.design)){
        std.o <- exp.design$std.order[k]
        pars.i <-  c(as.numeric(decode.design[k, c("rho", "sigma")]), beta0)
        if(std.o > 4){
            indx <- sample(n.samples, nb.sim)
            mcdatai <- list(simY = mc.data$simY[,indx], t10 = mc.data$t10,
                            t20 = mc.data$t20[indx])

            mc.res <- mcl.dCAR(pars = pars.i, data=data, simdata=mcdatai, rho.cons= NULL, Evar=TRUE)
            exp.design[k,]$mc.lr <- mc.res[1]
            exp.design[k,]$mc.var <- mc.res[2]
        }
        else{
            mc.res <- mcl.dCAR(pars = pars.i, data=data, simdata=mc.data, rho.cons=NULL, Evar=TRUE)
            exp.design[k,]$mc.lr <- mc.res[1]
            exp.design[k,]$mc.var <- mc.res[2]
        }
    }
    return(exp.design)
}

#### Direct CAR models (rho + sigma):
#### Evaluate along the steepest/canonical path
############################################################
ExpPath.dCAR <- function(exp.design, data, psi, mc.data){
    if(is.null(mc.data)){
       stop("Please provide the Monte Carlo samples by specifying mc.data!")
   }
    kk <- nrow(exp.design)
    mc.lr <- mc.var <- rep(0, kk)
    exp.design$mc.lr <- mc.lr
    exp.design$mc.var <- mc.var
    mc.res <- apply(exp.design[, c("rho", "sigma")], 1, function(x) mcl.dCAR(c(x, psi[-c(1,2)], rho.cons = NULL),
                    data=data, simdata=mc.data, Evar=TRUE))
    exp.design$mc.lr <- mc.res[1,]
    exp.design$mc.var <- mc.res[2,]
    return(exp.design)
}


ESA.dCAR <- function(exp.design, data, psi, mc.data, rho.r = c(-0.25, 0.25),
                      plot = FALSE, ds.only = TRUE){
    steep <- steepest(exp.design, dist = seq(0, 3, 0.5))
    steep <- steep[which((steep$rho > rho.r[1]) & (steep$rho < rho.r[2]) &
                         (steep$sigma > 0.01)), ]
    max.step <- max(steep$dist)
    expt.steep <- ExpPath.dCAR(steep, data = data, psi = psi, mc.data = mc.data)
    sa.val <- lm(mc.lr ~ poly(dist,2, raw = TRUE),  data = expt.steep)
    sa.var <- lm(log(mc.var) ~ poly(dist,2, raw = TRUE), data = expt.steep)
    coefs.val <- coef(sa.val)
    coefs.var <- coef(sa.var)
    opt.val <- as.numeric(-coefs.val[2]/coefs.val[3]/2)
    opt.var <- as.numeric(ifelse(coefs.var[3] == 0, max.step, -coefs.var[2]/coefs.var[3]/2))
    ds.val <- ifelse(opt.val > 0, min(max.step, round(opt.val, digits = 1)), max.step)
    ds.var <- min(max.step, floor(opt.var))
    ds <- ifelse(ds.var < 1, ds.val, min(ds.val, ds.var))
    pred.dist <-  seq(0, max.step, by = 0.1)
    pred.val <-  predict(sa.val, newdata = data.frame(dist =pred.dist))
    pred.var <-  predict(sa.var, newdata = data.frame(dist =pred.dist))

    if(plot){
        par(mfrow = c(1,2))
        plot(mc.lr ~ dist, data = expt.steep)
        lines(pred.dist, pred.val)
        plot(log(mc.var) ~ dist, data = expt.steep)
        lines(pred.dist, pred.var)
    }

    sa.res <- list(SA.fit = expt.steep, pred.dist = pred.dist,
                   pred.val=pred.val, pred.var = pred.var)
    if(ds.only){
        return(ds)
    }
    else{
        return(list(sa.res, ds = ds))
    }
}

#### latent CAR models (rho + sigma): -- Binomial or Poisson
#### Evaluate at the design points
############################################################
Exp.CAR.glm <- function(exp.design, data, family, psi, n0, mc.data){
    if(is.null(mc.data)){
       stop("Please provide the Monte Carlo samples by specifying mc.data!")
   }
    n.samples <- nrow(mc.data$Zy)
    nb.sim <- floor(n.samples/n0)

    mc.lr <- mc.var <- rep(0, 4+n0)
    exp.design$mc.lr <- mc.lr
    exp.design$mc.var <- mc.var
    beta0 <- as.numeric(psi[-c(1,2)])
    decode.design <- decode.data(exp.design)

    for (k in 1:nrow(exp.design)){
        std.o <- exp.design$std.order[k]
        pars.i <-  c(as.numeric(decode.design[k, c("rho", "sigma")]), beta0)
        if(std.o > 4){
            indx <- sample(n.samples, nb.sim)
            mcdatai <- mc.data
            mcdatai$Zy <- mcdatai$Zy[indx,]
            mcdatai$lHZy.psi <- mcdatai$lHZy.psi[indx]

            mc.res <- mcl.glm(pars = pars.i, mcdata = mcdatai, family = family, Evar = TRUE)
            exp.design[k,]$mc.lr <- mc.res[1]
            exp.design[k,]$mc.var <- ifelse(mc.res[2]==0, 1e-8, min(1e8, mc.res[2]))
        }
        else{
            mc.res <- mcl.glm(pars = pars.i, mcdata = mc.data, family = family, Evar = TRUE)
            exp.design[k,]$mc.lr <- mc.res[1]
            exp.design[k,]$mc.var <- ifelse(mc.res[2]==0, 1e-8, min(1e8, mc.res[2]))
        }
    }
    return(exp.design)
}

#### latent CAR models (rho + sigma): -- Binomial or Poisson
#### Evaluate points on steepest/canonical path
############################################################
ExpPath.CAR.glm<- function(exp.design, data, family, psi, mc.data){
   if(is.null(mc.data)){
       stop("Please provide the Monte Carlo samples by specifying mc.data!")
   }
    kk <- nrow(exp.design)
    mc.lr <- mc.var <- rep(0, kk)
    exp.design$mc.lr <- mc.lr
    exp.design$mc.var <- mc.var
    beta0 <- as.numeric(psi[-c(1,2)])
    pars <- exp.design[, c("rho", "sigma")]
    mc.res <- apply(pars, 1, function(x) mcl.glm(pars = c(x, beta0),
                                                     mcdata = mc.data, family = family, Evar = TRUE))
    exp.design$mc.lr <- mc.res[1,]
    exp.design$mc.var <-  ifelse(mc.res[2]==0, 1e-8, min(1e8, mc.res[2]))
    return(exp.design)
}

ESA.CAR.glm <- function(exp.design, data, family, psi, mc.data, rho.r = c(-0.25, 0.25),
                              plot = FALSE, ds.only = TRUE){
    steep <- steepest(exp.design, dist = seq(0, 3, 0.5))
    steep <- steep[which((steep$rho > rho.r[1]) & (steep$rho < rho.r[2]) &
                         (steep$sigma > 0.01)), ]
    max.step <- max(steep$dist)
    expt.steep <- ExpPath.CAR.glm(steep, data = data, family = family, psi = psi, mc.data = mc.data)
    sa.val <- lm(mc.lr ~ poly(dist,2, raw = TRUE),  data = expt.steep)
    sa.var <- lm(log(mc.var) ~ poly(dist,2, raw = TRUE), data = expt.steep)
    coefs.val <- coef(sa.val)
    coefs.var <- coef(sa.var)
    opt.val <- as.numeric(-coefs.val[2]/coefs.val[3]/2)
    opt.var <- as.numeric(ifelse(coefs.var[3] == 0, max.step, -coefs.var[2]/coefs.var[3]/2))
    ds.val <- ifelse(opt.val > 0, min(max.step, round(opt.val, digits = 1)), max.step)
    ds.var <- min(max.step, floor(opt.var))
    ds <- ifelse(ds.var < 1, ds.val, min(ds.val, ds.var))
    pred.dist <-  seq(0, max.step, by = 0.1)
    pred.val <-  predict(sa.val, newdata = data.frame(dist =pred.dist))
    pred.var <-  predict(sa.var, newdata = data.frame(dist =pred.dist))

    if(plot){
        par(mfrow = c(1,2))
        plot(mc.lr ~ dist, data = expt.steep)
        lines(pred.dist, pred.val)
        plot(log(mc.var) ~ dist, data = expt.steep)
        lines(pred.dist, pred.var)
    }

    sa.res <- list(SA.fit = expt.steep, pred.dist = pred.dist,
                   pred.val=pred.val, pred.var = pred.var)
    if(ds.only){
        return(ds)
    }
    else{
        return(list(sa.res, ds = ds))
    }
}

#### update psi according to ESA
Psi.ESA <- function(rsm.fitted, ds, psi.lim){ ## not to exceed the limit on the pars
    r.lim <- psi.lim$r.lim
    s.lim <- psi.lim$s.lim
    rs.new <- as.numeric(steepest(rsm.fitted, dist = ds)[,c("rho", "sigma")])
    rs.new[1] <- max(min(rs.new[1], r.lim[2]), r.lim[1])
    rs.new[2] <- max(min(rs.new[2], s.lim[2]), s.lim[1])
    rs.new
}

#### wrap the direct and latent CAR into one function
###########################################################
Exp.CAR <- function(exp.design, data, psi, n0, mc.data, family){
    if(family == "gauss"){
        Exp.dCAR(exp.design, data, psi, n0, mc.data)
    }else{
        Exp.CAR.glm(exp.design, data, family, psi, n0, mc.data)
    }
}

ExpPath.CAR <- function(exp.design, data, psi, mc.data, family){
     if(family == "gauss"){
        ExpPath.dCAR(exp.design, data, psi, mc.data)
    }else{
        ExpPath.CAR.glm(exp.design, data, family, psi, mc.data)
    }
}

ESA.CAR <- function(exp.design, data, psi, mc.data, family, rho.r = c(-0.25, 0.25),
                    plot = FALSE, ds.only = TRUE){
     if(family == "gauss"){
        ESA.dCAR(exp.design, data, psi, mc.data, rho.r, plot, ds.only)
    }else{
        ESA.CAR.glm(exp.design, data, family, psi, mc.data, rho.r, plot, ds.only)
    }
}

#### Plot the fitted response
#### exact values can be plotted for direct CAR models but
#### not for latent!!!
############################################################
plot.RSM <- function(rsm.fitted, psi, bounds, fvbox = NULL, DAbox = NULL,
                     rho.r = c(-0.25, 0.25),
                     A.path = c("steepest", "canonical"), distance,
                     exact = NULL, exact.vals = TRUE){
    rsm.res <- contour(rsm.fitted, ~ x1 + x2, bounds = bounds, plot.it = FALSE)
    x.r <- range(rsm.res[[1]]$x)
    x.lim <- c(max(rho.r[1], x.r[1]), min(x.r[2], rho.r[2]))
    y.r <- range(rsm.res[[1]]$y)
    y.lim <- c(max(1e-4, y.r[1]), min(y.r[2], 2*psi[2]))
    levs <- seq(-15, 3, 3)

    ## When rho.r > 0.25, z = -Inf and sigma.r > 2*psi[[2]]
    ind.r <- (rsm.res[[1]]$x <= rho.r[2]) & (rsm.res[[1]]$x >= rho.r[1])
    ind.s <- rsm.res[[1]]$y >= 2*psi[2]
    rsm.res[[1]]$z[!ind.r,] <- -Inf
    rsm.res[[1]]$z[,ind.s] <- -Inf


    if(exact.vals){
        if(is.null(exact)){
            stop("Exact values can not be NULL.")
        }
        rhos <- exact$r.vec
        sigmas <- exact$s.vec
        rho.range <- range(rhos)
        sigma.range <- range(sigmas)
        x.lim.new <- c(max(rho.range[1], x.lim[1]), min(rho.range[2], x.lim[2]))
        y.lim.new <- c(max(sigma.range[1], y.lim[1]), min(sigma.range[2], y.lim[2]))
        lrs <- exact$lrs
        levs.new <- seq(-15, 0, 3)
        ## image of the exact values
        image.plot(rhos, sigmas, lrs, xlim = x.lim.new, ylim = y.lim.new,
                   breaks = c(levs.new, max(lrs)), nlevel = length(levs.new)+1,
                   col = terrain.colors(length(levs.new), alpha = 0.7),
                   xlab = expression(rho), ylab = expression(sigma^2))
        ## Contour of the fitted RSM
        contour(rsm.res[[1]], levels = levs, add = TRUE)
    }
    else{
        ## Contour of the fitted RSM
        contour(rsm.res[[1]], levels = levs, xlim = x.lim, ylim = y.lim,
                xlab = expression(rho), ylab = expression(sigma^2))
    }

    ## finite variance area
    if(!is.null(fvbox)){
        polygon(fvbox[,1],fvbox[,2], border = "red")
    }
    ## Plot the design region
    if(!is.null(DAbox)){
        polygon(DAbox[,1], DAbox[,2], border = "red")
        p.centre <- c(mean(DAbox[,1]), mean(DAbox[,2]))
        points(p.centre[1],  p.centre[2])
    }
    ## Plot the importance sampler
    points(psi[1], psi[2], col = "blue", cex = 1, pch = 15) # design centre

    ## Plot the canonical path
    if ("canonical" %in% A.path){
        ## plot the stationary point
        xss <- code2val(xs(rsm.fitted), codings(rsm.fitted))
        points(xss[1],  xss[2], col = 2, cex = 1, pch = 16)

        ## Plot the canonical path
        canonical.data <- canonical.path(rsm.fitted, dist = seq(-3,3, by = 0.5))
        rs.path <- canonical.data[, c("rho", "sigma")]
        n.points <- nrow(rs.path)
        nss <- seq(n.points, 2, by = -2)
        arrows(x0 = rs.path[nss,1], y0 = rs.path[nss,2],
               x1 = rs.path[nss-1,1], y1 = rs.path[nss-1,2],
               length = 0.07, col = "red")
    }
    ## Plot the Steepest Ascent path
    if("steepest" %in% A.path){
        steepest.data <- steepest(rsm.fitted)
        rs.path <- steepest.data[, c("rho", "sigma")]
        n.points <- nrow(rs.path)
        nss <- seq(2, n.points, by = 2)
        arrows(x0 = rs.path[nss-1,1], y0 = rs.path[nss-1,2],
               x1 = rs.path[nss,1], y1 = rs.path[nss,2],
               length = 0.07, col = "red")
        ss.next <- as.numeric(steepest(rsm.fitted, dist = distance)[, c("rho", "sigma")])
        points(ss.next[1], ss.next[2], col = 2, pch = 17) #next centre
    }
}
