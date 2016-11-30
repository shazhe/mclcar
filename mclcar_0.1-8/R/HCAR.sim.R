#### Functions for simulating samples from the importance distribution of a HCAR
#### sampling eta | y
sim.HCAR <- function(psi, data, n.samples){
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

    if(is.null(W)){
        lambda <- psi[1]
        sigma.e <- psi[2]
        sigma.u <- psi[3]
        beta <- psi[-c(1:3)]

        Qe.s <- In/sigma.e #as.spam(In/sigma.e)
        Qu.s <- (Ik - lambda*M)/sigma.u # as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        Q.new <- crossprod(Z)/sigma.e  + Qu.s #as.spam(crossprod(Z)/sigma.e  + Qu.s)
        b <-  solve(Q.new, crossprod(Z, as.spam(Qe.s %*% (y-Xb))))
        CholQ <- chol(Q.new,  pivot = "MMD", memory = list(nnzcolindices = 6.25 * K))
        z.err <- matrix(rnorm(K*n.samples), nrow = K, ncol = n.samples)
        u.y <- apply(z.err, 2, function(x) b + backsolve(CholQ, x))

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res)/sigma.e + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpsi <- apply(u.y, 2, log.uy)
        logdetQes.half <- sum(log(diag(Qe.s)))/2
        Lu.s <-  chol(Qu.s, pivot = "MMD", memory = list(nnzcolindices = 6.25 * K))
        logdetQus.half <- sum(log(diag(Lu.s)))
        const.psi <- logdetQes.half + logdetQus.half


    }else{
        rho <- psi[1]
        lambda <- psi[2]
        sigma.e <- psi[3]
        sigma.u <- psi[4]
        beta <- psi[-c(1:4)]

        Qe.s <- as.spam((In - rho*W)/sigma.e)
        Qu.s <- as.spam((Ik - lambda*M)/sigma.u)
        Xb <- X %*% beta

        Q.new <- as.spam(crossprod(Z, Qe.s) %*% Z  + Qu.s)
        b <-  solve(Q.new, crossprod(Z, Qe.s %*% (y-Xb)))
        CholQ <- chol(Q.new,  pivot = "MMD", memory = list(nnzcolindices = 6.25 * K))
        z.err <- matrix(rnorm(K*n.samples), nrow = K, ncol = n.samples)
        u.y <- apply(z.err, 2, function(x) b + backsolve(CholQ, x))

        log.uy <- function(u.y){
            res <- y - Xb - Z%*%u.y
            as.numeric(-0.5 *(crossprod(res, Qe.s) %*% res + crossprod(u.y, Qu.s) %*% u.y))
        }

        lpsi <- apply(u.y, 2, log.uy)

        Le.s <-  chol(Qe.s, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n))
        logdetQes.half <- sum(log(diag(Le.s)))
        Lu.s <-  chol(Qu.s, pivot = "MMD", memory = list(nnzcolindices = 6.25 * K))
        logdetQus.half <- sum(log(diag(Lu.s)))
        const.psi <- logdetQes.half + logdetQus.half
    } # W

    return(list(psi = psi, n.samples = n.samples, u.y = u.y, lpsi = lpsi, const.psi = const.psi))

}
