################################################################
####     RSM Optimisation using an iterative procedure      ####
################################################################
rsmMCL <- function(data, psi0, family, exact0=NULL, control = list(), mc.control = list()){
    N.data <- nrow(data$W)
    if(family == "gauss"){
        ## general algorithm default control
        con <- list(mc.samples = c(500, 2, 1e4), n.iter = 20, time.max = 2,
                    K = c(6, 3), mc.var = FALSE,
                    psi.lim = list(r.lim = c(-0.2499, 0.2499), s.lim = c(0.1, 8)),
                    exacts = list(eval = TRUE, rho = c(-0.25, 0.25), sigma = c(0.5, 2), length = 100),
                    rsm.fit = list(n01 = 4, n02= 2, Rsq = 0.9, lof = 5e-2,
                    st.diff = 50, mcl.diff = 10), plotrsm = TRUE, trace.all = TRUE, verbose = TRUE)
        con[names(control)] <- control

        n.s <- n.s0 <- con$mc.samples[1]
        n.increase <- con$mc.samples[2]
        n.max <- con$mc.samples[3]
        rho.lim <- con$psi.lim$r.lim
        sigma.lim <- con$psi.lim$s.lim
        rho.r0 <- con$exacts$rho
        sigma.r0 <- con$exacts$sigma
        vec.length <- con$exacts$length
		mc.var <- con$mc.var
		plotrsm=con$plotrsm
        r.vec <- seq(rho.r0[1] + 0.01, rho.r0[2] - 0.01, length.out = vec.length)
        s.vec <- seq(sigma.r0[1], sigma.r0[2], length.out = vec.length)
        rs.mat <- expand.grid(r.vec, s.vec)
        Exacts <- list()
    }else{
        ## MCMC default control
        mc.con <- list(method = "mala", N.Zy = 1e4, Scale = 1.65/(N.data^(2/6)), thin = 5,
                       burns = 5e3, scale.fixed = TRUE, c.ad = c(1, 0.7))
        mc.con[names(mc.control)] <- mc.control

        ## general algorithm default control
        con <- list(mc.samples = c(2, 1e4), n.iter = 20, time.max = 2, K = c(6, 3), mc.var = FALSE,
                    psi.lim = list(r.lim = c(-0.2499, 0.2499), s.lim = c(0.1, 8)),
                    exacts = list(eval = FALSE),
                    rsm.fit = list(n01 = 4, n02= 2, Rsq = 0.9, lof = 5e-2, st.diff = 50, mcl.diff = 10),
                    plotrsm = TRUE, trace.all = TRUE, verbose = TRUE)
        con[names(control)] <- control

        n.s <- n.s0 <- (mc.con$N.Zy - mc.con$burns)/mc.con$thin
        n.increase <- con$mc.samples[1]
        n.max <- con$mc.samples[2]
        rho.lim <- con$psi.lim$r.lim
        sigma.lim <- con$psi.lim$s.lim
        
        plotrsm=con$plotrsm
    }

    ## Initializing
    kk <- con$K[1]
    kk.decrease <- con$K[2]
    mcl.diff <- con$rsm.fit$mcl.diff
    n.iter <- con$n.iter
    time.max <- con$time.max
    psi <- psi0

    Psi  <- list()                      # optimal points in each iteration
    FVbox <- list()                     # Finite variance boxes
    DAs <- list()                       # transformed Design areas
    DAbox <- list()                     # Design area box
    Ttime <- list()                     # Time for one iteration

    Expts <- list()            # generated design points
    Expt.Val1 <- list()        # evaluation at 1st order points
    Expt.Val2 <- list()        # evaluation at 2nd order points
    RSM.1 <- list()            # 1st order rsm
    RSM.2 <- list()            # 2nd order rsm
    ESAs <- list()             # Steepest ascent analysis
    Sim.data <- list()         # Monte Carlo samples simulated in each iteration


    flag.converge <- FALSE
    flag.stationary <- FALSE
    sDA <- FALSE
    n.fst = 0
    i = 1
    t.total <- 0
    par(mfrow = c(4, 3))

    ## Start rsm iteration
    while (i <= n.iter & (!flag.stationary | !flag.converge) & !sDA){

        if(n.s > n.max){
            warning("Iteration terminated: Evaluations can take too long with the maximum Monte Carlo sample size.")
            break
        }

        if(t.total/3600 > time.max){
            warning("Iteration terminated: Maximum evaluation time reached.")
            break
        }

        if(con$exacts$eval){
            if(family == "gauss"){
                betai <- psi[-c(1,2)]
                lls <- matrix(apply(rs.mat, 1, function(x) loglik.dCAR(c(x, betai), data = data)),
                              vec.length, vec.length)
                exactsi <- list(r.vec = r.vec, s.vec = s.vec, lrs = lls)
                exactsi$lrs <- exactsi$lrs - loglik.dCAR(psi, data = data)
                Exacts[[i]] <- exactsi
            }else{
                stop("No exact value for non-Guassin data!")
            }
        }else{
            exactsi <- exact0
        }


        t.time <- system.time({
            ## Design area
            fvbox <- data.frame(rhos = c(rho.lim[1], rho.lim[2], psi[1]),
                                sigmas = c(0, 0, psi[2]*2))
            da <- range.tt(psi, N = N.data, s = n.s, rho.r = rho.lim, k = kk)
            sigma.r <- psi[2]*c(exp(-da[3]), exp(da[3]))
            rho.r <- psi[1] +(da[1:2])
            dabox <- data.frame(rhos = c(rho.r, rho.r[2], rho.r[1]),
                                sigmas = c(sigma.r[1], sigma.r[1], sigma.r[2], sigma.r[2]))

            ## Generate the design points
            expts <- RSM.design.CAR(design.tt = da, psi = psi)

            ## Generate MC samples given current psi
            if(family == "gauss"){
                sim.data <- mcl.prep.dCAR(psi, n.s, data)
            }else{
                sim.data <- mcl.prep.glm(data, family, psi, mcmc.control = mc.con)
            }

            ## fit 1st order model
            expt.1 <- Exp.CAR(expts[[1]], data = data, psi = psi, n0 = con$rsm.fit$n01,
                              mc.data = sim.data, family)

            rsm.1 <- rsm(mc.lr ~ FO(x1, x2), data = expt.1, weight = 1/mc.var)

            ## Save the results 1
            FVbox[[i]] <- fvbox
            DAs[[i]] <- da
            DAbox[[i]] <- dabox
            Expts[[i]] <- expts
            Expt.Val1[[i]] <- expt.1
            RSM.1[[i]] <- rsm.1
            Sim.data[[i]] <- sim.data

            ## goodness of fit of the 1st order model
            Rsq1 <- summary(rsm.1)$r.squared
            lof.p <- summary(rsm.1)$lof[[5]][3]
            max.mcl <- max(abs(expt.1$mc.lr))

            ## fit 2nd order
            if(Rsq1 < con$rsm.fit$Rsq & lof.p < con$rsm.fit$lof){
                expt.2 <- Exp.CAR(expts[[2]], data = data, psi = psi, n0 = con$rsm.fit$n02,
                                  mc.data = sim.data, family)
                rsm.2 <- rsm(mc.lr ~ SO(x1, x2), data = djoin(expt.1, expt.2),
                             weights = 1/mc.var)
                max.mcl <- max(c(abs(expt.2$mc.lr), max.mcl))
                st.point.c <- xs(rsm.2)
                flag.in <- all(abs(st.point.c) < 1)
                flag.stationary <- flag.in & mean(abs(c(expt.1$mc.lr, expt.2$mc.lr))) <
                    con$rsm.fit$st.diff

                if(flag.stationary){
                    st.point <- as.numeric(code2val(st.point.c,
                                                    codings(expts[[2]]))[c("rho", "sigma")])
                    if(plotrsm){
                    plot.RSM(rsm.2, psi = psi, bounds = list(x1=c(-2,2), x2 = c(-2,2)),
                             fvbox = fvbox, DAbox = dabox, A.path = "canonical", distance = 0,
                             exact = exactsi, exact.vals = con$exacts$eval)
                }
                    rs.new <- st.point
                }else{
                    esa.res <- ESA.CAR(rsm.2, data = data, psi = psi, mc.data = sim.data, family,
                                       plot = plotrsm, ds.only = FALSE)
                    ds <- esa.res$ds
                    ESAs[[i]] <- esa.res
                    rs.newESA <- Psi.ESA(rsm.2, ds, con$psi.lim)
                    ## Control the error at the boundary of the parameter
                    flag.boundary <- min(abs((rs.newESA[1]-rho.lim)/rho.lim)) < 5e-3 | min((abs(rs.newESA[2] - sigma.lim)/sigma.lim)) < 5e-3
                    if(flag.boundary){
                        rs.new <- (ds/3.5)*rs.newESA + (1-ds/3.5)*psi[1:2]
                    }else{
                        rs.new <- rs.newESA
                    }
                    if(plotrsm){
                    plot.RSM(rsm.2, psi = psi, bounds = list(x1=c(-2,2), x2 = c(-2,2)),
                             fvbox = fvbox, DAbox = dabox, A.path = "steepest", distance = ds,
                             exact = exactsi, exact.vals = con$exacts$eval)
                    }
                }
                Expt.Val2[[i]] <- expt.2
                RSM.2[[i]] <- rsm.2
            }else{
                esa.res <- ESA.CAR(rsm.1, data = data, psi = psi, mc.data = sim.data, family,
                                   plot = plotrsm, ds.only = FALSE)
                ds <- esa.res$ds
                ESAs[[i]] <- esa.res
                rs.newESA <- Psi.ESA(rsm.1, ds, con$psi.lim)
                ## Control the error at the bounary of the parameter
                flag.boundary <- min(abs((rs.newESA[1]-rho.lim)/rho.lim)) < 1e-4 | min((abs(rs.newESA[2] - sigma.lim)/sigma.lim)) < 1e-4
                if(flag.boundary){
                        rs.new <- (ds/3.1)*rs.newESA + (1-ds/3.1)*psi[1:2]
                    }else{
                        rs.new <- rs.newESA
                    }
                if(plotrsm){
                plot.RSM(rsm.1, psi = psi, bounds = list(x1=c(-2,2), x2 = c(-2,2)),
                         fvbox = fvbox, DAbox = dabox, A.path = "steepest", distance = ds,
                         exact = exactsi, exact.vals = con$exacts$eval)
            }
            }
            ## Update the psi and simulate new MC samples
            if(family == "gauss"){
                beta.new <- get.beta.lm(rs.new, data)
            }else{
                beta.new <- get.beta.glm(psi[-c(1,2)], rs.new, sim.data, family)$x
            }
            psi.new <- c(rs.new, beta.new)
            Psi[[i]] <- psi.new

            diff.psi <- abs(psi.new - psi)
            flag.diff <- all(diff.psi < 1/n.s)
            flag.mcl <- round(max.mcl) <= mcl.diff
            psi <- psi.new
        })


        sDA <- kk < con$K[1] & flag.converge
        flag.converge <- (flag.mcl | flag.diff) & (n.s > n.s0) & flag.stationary
        n.fst <- ifelse(flag.stationary, n.fst + 1, n.fst)
        kk <- ifelse(n.fst == 1 | flag.converge, kk/kk.decrease, kk)
        Ttime[[i]] <- t.time
        t.total <- t.total + t.time[3]
        if(con$verbose){
            print(paste(c("iteration =", "MC.samples =", "stationary =",
                          "converged =", "time.elapsed ="),
                        c(i, n.s, flag.stationary, flag.converge, t.total)))
        }
        if((flag.mcl | flag.diff) & n.s == n.s0){
            n.s <- n.s * n.increase
            if(family != "gauss"){
                mc.con$N.Zy <- n.s * mc.con$thin + mc.con$burns
            }
        }
        i = i + 1
    }

    if(family != "gauss"){
        Exacts <- rep(NA, i - 1)
    }

    if(con$trace.all){
        ans <- list(psi0= psi0, Psi=Psi, FVbox = FVbox, DAs=DAs, DAbox=DAbox, Ttime=Ttime,
                    Expts = Expts, Expt.Val1=Expt.Val1, Expt.Val2=Expt.Val2,
                    RSM.1 = RSM.1, RSM.2 = RSM.2, ESAs = ESAs, Sim.data = Sim.data, Exacts = Exacts,
                    data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "rsmMCL"
        return(ans)
    }else{
        ans <- list(psi0= psi0, Psi = Psi[(i-2):(i-1)], FVbox = FVbox[i-1], DAs = DAs[i-1],
                    DAbox = DAbox[i-1], Ttime = Ttime[i-1], Expts = Expts[i-1],
                    Expt.Val1 = Expt.Val1[i-1], Expt.Val2 = Expt.Val2[length(Expt.Val2)],
                    RSM.1 = RSM.1[i-1], RSM.2 = RSM.2[length(RSM.2)],
                    ESAs = ESAs[length(ESAs)], Sim.data = Sim.data[i-1], Exacts = Exacts[i-1],
                    data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "rsmMCL"
        return(ans)
    }
}
