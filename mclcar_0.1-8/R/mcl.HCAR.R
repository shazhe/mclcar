#### Functions for evaluating the Monte Carlo likelihood for HCAR models

## The Monte Carlo likelihood
mcl.HCAR <- function(pars, mcdata, data, Evar = FALSE){
  if (is.environment(data)) {
    y <- get("y", envir=data)
    X <- get("X", envir=data)
    W <- get("W", envir=data)
    M <- get("M", envir=data)
    Z <- get("Z", envir=data)
    n <- length(y)
    K <- nrow(M)
    In <- get("In", envir=data) #diag(1, n)
    Ik <- get("Ik", envir=data) #diag(1, K)
   } else {
    y <- data$y
    X <- data$X
    W <- data$W
    M <- data$M
    Z <- data$Z

    n <- length(y)
    K <- nrow(M)
    In <- data$In #diag(1, n)
    Ik <- data$Ik #diag(1, K)
  }
    lpsi <- mcdata$lpsi
    const.psi <- mcdata$const.psi
    u.y <- mcdata$u.y
    n.samples <- mcdata$n.samples


    if(is.null(W)){ # W NULL
        lambda <- pars[1]
        sigma.e <- pars[2]
        sigma.u <- pars[3]
        beta <- pars[-c(1:3)]

        Qe.s <- In/sigma.e #as.spam(In/sigma.e)
        Qu.s <- (Ik - lambda*M)/sigma.u #as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res)/sigma.e + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpars <- apply(u.y, 2, log.uy)
        logdetQes.half <- sum(log(diag(Qe.s)))/2
        Lu.s <-  chol(Qu.s, pivot = "MMD", memory = list(nnzcolindices = 6.25 * K))
        logdetQus.half <- sum(log(diag(Lu.s)))
        const.pars <- logdetQes.half + logdetQus.half
    }else{ # W NULL
        rho <- pars[1]
        lambda <- pars[2]
        sigma.e <- pars[3]
        sigma.u <- pars[4]
        beta <- pars[-c(1:4)]

        Qe.s <- as.spam((In - rho*W)/sigma.e)
        Qu.s <- as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res, Qe.s) %*% res + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpars <- apply(u.y, 2, log.uy)

        Le.s <-  chol(Qe.s, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n))
        logdetQes.half <- sum(log(diag(Le.s)))
        Lu.s <-  chol(Qu.s, pivot = "MMD", memory = list(nnzcolindices = 6.25 * K))
        logdetQus.half <- sum(log(diag(Lu.s)))
        const.pars <- logdetQes.half + logdetQus.half
    } #W

    Lrs <- exp(lpars + const.pars - lpsi - const.psi)
    mc.lr <- log(mean(Lrs))

    if (abs(mc.lr) == Inf){
        warning("Design area too large: Inf produced! Choose a parameter closer to psi.")
        mc.lr <- sign(mc.lr)*1e5
    }

    if (Evar){
        if(abs(mc.lr) == Inf){
            v.lr <- 100
        }else{
            v.lr <- var(Lrs/exp(mc.lr))/n.samples
            return(c(mc.lr, v.lr))
        }
    }
    else{
        return(mc.lr)
    }
}

## const.pars <- sum(log(1-rho*aa)) + sum(log(1-lambda*bb))-n/2*log(sigma.e) - K/2*log(sigma.u)
##  the gradient of the Monte Carlo likelihood
mcl.grad.HCAR <- function(pars, mcdata, data, Evar = FALSE){
  if (is.environment(data)) {
    y <- get("y", envir=data)
    X <- get("X", envir=data)
    W <- get("W", envir=data)
    M <- get("M", envir=data)
    Z <- get("Z", envir=data)
    n <- length(y)
    K <- nrow(M)
    In <- get("In", envir=data) #diag(1, n)
    Ik <- get("Ik", envir=data) #diag(1, K)
    bb <- get("bb", envir=data)
   } else {
    y <- data$y
    X <- data$X
    W <- data$W
    M <- data$M
    Z <- data$Z

    n <- length(y)
    K <- nrow(M)
    In <- data$In #diag(1, n)
    Ik <- data$Ik #diag(1, K)
    bb <- data$bb
  }
    lpsi <- mcdata$lpsi
    u.y <- mcdata$u.y
    n.samples <- mcdata$n.samples

    if(is.null(W)){
        lambda <- pars[1]
        sigma.e <- pars[2]
        sigma.u <- pars[3]
        beta <- pars[-c(1:3)]
        if(is.null(bb)){
            bb <- eigen(M)$values
        }
        Qe.s <- In/sigma.e #as.spam(In/sigma.e)
        Qu.s <- (Ik - lambda*M)/sigma.u #as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res)/sigma.e + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpars <- apply(u.y, 2, log.uy)

        grad.uy <- function(u.y){
            Z.uy <- Z %*% u.y
            ee <- y -Xb - Z.uy
            g.beta <- crossprod(X, ee)/sigma.e
            g.sigma.e <- crossprod(ee)/(2*sigma.e^2) - n/(2*sigma.e)
            g.lambda <- crossprod(u.y, M) %*% u.y/(2*sigma.u) - sum(bb/(1-lambda*bb))/2
            g.sigma.u <- crossprod(u.y, Qu.s) %*% u.y/(2*sigma.u) - K/(2*sigma.u)
            c(g.lambda, g.sigma.e, g.sigma.u, g.beta)
        }
    }else{
         rho <- pars[1]
        lambda <- pars[2]
        sigma.e <- pars[3]
        sigma.u <- pars[4]
        beta <- pars[-c(1:4)]


        if(is.null(data$aa)){
            aa <- eigen(W)$values
        }else{
            aa <- data$aa}
        if(is.null(data$bb)){
            bb <- eigen(M)$values
        }else{
            bb <-data$bb}

        Qe.s <- as.spam((In - rho*W)/sigma.e)
        Qu.s <- as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res, Qe.s) %*% res + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpars <- apply(u.y, 2, log.uy)


        grad.uy <- function(u.y){
            Z.uy <- Z %*% u.y
            ee <- y -Xb - Z.uy
            g.beta <- crossprod(X, Qe.s) %*% ee
            g.rho <- crossprod(ee, W) %*% ee/(2*sigma.e) - sum(aa/(1-rho*aa))/2
            g.sigma.e <- crossprod(ee, Qe.s) %*% ee /(2*sigma.e) - n/(2*sigma.e)
            g.lambda <- crossprod(u.y, M) %*% u.y/(2*sigma.u) - sum(bb/(1-lambda*bb))/2
            g.sigma.u <- crossprod(u.y, Qu.s) %*% u.y/(2*sigma.u) - K/(2*sigma.u)
            c(g.rho, g.lambda, g.sigma.e, g.sigma.u, g.beta)
        }
     }

    grad.coef <- apply(u.y, 2, grad.uy)
    Lrs <- exp(lpars - lpsi)
    grad.w <- Lrs/mean(Lrs)
    grad.lrs <- grad.coef * grad.w
    grad.lr <- rowMeans(grad.lrs)

    if (Evar){
        v.lr <- var(t(grad.lrs - grad.lr))/n.samples
        return(list(grad.lr = grad.lr, grad.var = v.lr))
    }
    else{
        return(grad.lr)
    }

}

##  the Hessian of the Monte Carlo likelihood
mcl.Hessian.HCAR <- function(pars, mcdata, data){
  if (is.environment(data)) {
    y <- get("y", envir=data)
    X <- get("X", envir=data)
    W <- get("W", envir=data)
    M <- get("M", envir=data)
    Z <- get("Z", envir=data)
    n <- length(y)
    K <- nrow(M)
    In <- get("In", envir=data) #diag(1, n)
    Ik <- get("Ik", envir=data) #diag(1, K)
    bb <- get("bb", envir=data)
   } else {
    y <- data$y
    X <- data$X
    W <- data$W
    M <- data$M
    Z <- data$Z

    n <- length(y)
    K <- nrow(M)
    In <- data$In #diag(1, n)
    Ik <- data$Ik #diag(1, K)
    bb <- data$bb
}
    lpsi <- mcdata$lpsi
    u.y <- mcdata$u.y
    n.samples <- mcdata$n.samples
    n.pars <- length(pars)

    if(is.null(W)){
        lambda <- pars[1]
        sigma.e <- pars[2]
        sigma.u <- pars[3]
        beta <- pars[-c(1:3)]

        if(is.null(bb)){
            bb <- eigen(M)$values
        }

        Qe.s <- In/sigma.e #as.spam(In/sigma.e)
        Qu.s <- (Ik - lambda*M)/sigma.u #as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res)/sigma.e + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpars <- apply(u.y, 2, log.uy)


        grad.uy <- function(u.y){
            Z.uy <- Z %*% u.y
            ee <- y -Xb - Z.uy
            g.beta <- crossprod(X,ee)/sigma.e
            g.sigma.e <- crossprod(ee) /(2*sigma.e^2) - n/(2*sigma.e)
            g.lambda <- crossprod(u.y, M) %*% u.y/(2*sigma.u) - sum(bb/(1-lambda*bb))/2
            g.sigma.u <- crossprod(u.y, Qu.s) %*% u.y/(2*sigma.u) - K/(2*sigma.u)
            c(g.lambda, g.sigma.e, g.sigma.u, g.beta)
        }

        grad.coef <- apply(u.y, 2, grad.uy)

        H.beta2 <- - crossprod(X,Qe.s) %*% X
        H.lambda2 <- - sum((bb/(1-lambda*bb))^2)/2

        Hessian.uy <- function(u.y){
            Z.uy <- Z %*% u.y
            ee <- y -Xb - Z.uy
            H.mat <- matrix(0, n.pars, n.pars)

            H.sige2 <- - crossprod(ee)/(sigma.e^3) + n/(2*sigma.e^2)
            H.sigu2 <- -crossprod(u.y, Qu.s) %*% u.y/(sigma.u^2) + K/(2*sigma.u^2)
            H.lambda.sigu <- - crossprod(u.y, M) %*% u.y/(2*sigma.u^2)
            H.sig.beta <- - crossprod(X, ee)/sigma.e^2

            H.mat[1,1] <- H.lambda2
            H.mat[1,3] <- H.mat[3,1] <- H.lambda.sigu
            H.mat[2,2] <- H.sige2
            H.mat[2,4:n.pars] <- H.mat[4:n.pars,2] <- H.sig.beta
            H.mat[3,3] <- H.sigu2
            H.mat[4:n.pars, 4:n.pars] <- H.beta2
            H.mat
        }
    }else{
        rho <- pars[1]
        lambda <- pars[2]
        sigma.e <- pars[3]
        sigma.u <- pars[4]
        beta <- pars[-c(1:4)]


        if(is.null(data$aa)){
            aa <- eigen(W)$values
        }else{
            aa <- data$aa}
        if(is.null(data$bb)){
            bb <- eigen(M)$values
        }else{
            bb <-data$bb}

        Qe.s <- as.spam((In - rho*W)/sigma.e)
        Qu.s <- as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res, Qe.s) %*% res + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpars <- apply(u.y, 2, log.uy)


        grad.uy <- function(u.y){
            Z.uy <- Z %*% u.y
            ee <- y -Xb - Z.uy
            g.beta <- crossprod(X, Qe.s) %*% ee
            g.rho <- crossprod(ee, W) %*% ee/(2*sigma.e) - sum(aa/(1-rho*aa))/2
            g.sigma.e <- crossprod(ee, Qe.s) %*% ee /(2*sigma.e) - n/(2*sigma.e)
            g.lambda <- crossprod(u.y, M) %*% u.y/(2*sigma.u) - sum(bb/(1-lambda*bb))/2
            g.sigma.u <- crossprod(u.y, Qu.s) %*% u.y/(2*sigma.u) - K/(2*sigma.u)
            c(g.rho, g.lambda, g.sigma.e, g.sigma.u, g.beta)
        }

        grad.coef <- apply(u.y, 2, grad.uy)

        H.beta2 <- - crossprod(X,Qe.s) %*% X
        H.lambda2 <- - sum((bb/(1-lambda*bb))^2)/2
        H.rho2 <- - sum((aa/(1-rho*aa))^2)/2

        Hessian.uy <- function(u.y){
            Z.uy <- Z %*% u.y
            ee <- y -Xb - Z.uy
            H.mat <- matrix(0, n.pars, n.pars)

            H.sige2 <- - crossprod(ee, Qe.s) %*% ee/(sigma.e^2) + n/(2*sigma.e^2)
            H.sigu2 <- -crossprod(u.y, Qu.s) %*% u.y/(sigma.u^2) + K/(2*sigma.u^2)
            H.rho.sige <- - crossprod(ee, W) %*% ee/(2*sigma.e^2)
            H.lambda.sigu <- - crossprod(u.y, M) %*% u.y/(2*sigma.u^2)
            H.rho.beta <- - crossprod(X, Qe.s) %*% ee
            H.sig.beta <- - crossprod(X, Qe.s) %*% ee/sigma.e

            H.mat[1,1] <- H.rho2
            H.mat[1,3] <- H.mat[3,1] <- H.rho.sige
            H.mat[1,5:n.pars] <- H.mat[5:n.pars, 1] <- H.rho.beta
            H.mat[2,2] <- H.lambda2
            H.mat[2,4] <- H.mat[4,2] <- H.lambda.sigu
            H.mat[3,3] <- H.sige2
            H.mat[3,5:n.pars] <- H.mat[5:n.pars, 3] <- H.sig.beta
            H.mat[4,4] <- H.sigu2
            H.mat[5:n.pars, 5:n.pars] <- H.beta2
            H.mat
        }
    }

    Lrs <- exp(lpars - lpsi)
    W.Lrs <- Lrs/mean(Lrs)

    ## -(grad(Z)*WZ)^2
    grad.wi <- sapply(1:n.samples, function(x) grad.uy(u.y[,x]) * W.Lrs[x])
    grad.w <- rowMeans(grad.wi)
    H.part1 <- -crossprod(t(grad.w))
    ## mean(grad(Z)^2 * WZ)
    grad.2wi <- sapply(1:n.samples, function(x) crossprod(t(grad.uy(u.y[,x]))) * W.Lrs[x])
    grad.2w <- rowMeans(grad.2wi)
    H.part2 <- matrix(grad.2w, n.pars, n.pars)
    ## mean(Hess(Z)*WZ)
    Hess.wi <- sapply(1:n.samples, function(x) Hessian.uy(u.y[,x]) * W.Lrs[x])
    Hess.w <- rowMeans(Hess.wi)
    H.part3 <-  matrix(Hess.w, n.pars, n.pars)

    H <- H.part1 + H.part2 + H.part3
    H

}
## Given rho and sigma, find beta
get.beta.HCAR <- function(beta0, rho.sig, mcdata, data){
    grad.beta <- function(beta){
        pars.b <- c(rho.sig, beta)
        n.rs <- length(rho.sig)
        grad.b <- mcl.grad.HCAR(pars.b, mcdata, data)[-c(1:n.rs)]
        grad.b <- ifelse(abs(grad.b) == Inf, sign(grad.b) * 1e5, grad.b)
        if (any(abs(beta0 - beta) > 5)){ # penalise when beta is too far away from the beta0
            return(sign(grad.b)*1e5)
        }else{
            return(grad.b)
        }
    }
    nleqslv(beta0, grad.beta, method = "Newton", control= list(ftol = 1, maxit = 3 ))$x
}


##  Variance of the Monte Carlo MLE
vmle.HCAR <- function(MLE, mcdata, data){
    par <- MLE$estimate
    A <- mcl.grad.HCAR(par, mcdata, data, Evar = TRUE)$grad.var
    B <- MLE$hessian
    Binv <- solve(B)
    Binv %*% A %*% Binv # The matrix are usually small so just compute directly
}
