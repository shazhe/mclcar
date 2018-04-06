#### Simulate CAR data on n1 * n2 torus with rook style neighbours and given precision
CAR.simTorus <- function(n1, n2, rho, prec){
    ## build the weight matrix -- block circulant matrix
    x <- list(c(2, n2), rep(1,2))
    Wx <- circulant.spam(x,n = n2)
    Ix <- diag.spam(1,n1)
    W <- kronecker(Ix,Wx) + kronecker(Wx,Ix)
    ew <- eigen(W,only.values = TRUE)$values
    min <-1/min(ew)
    max <- 1/max(ew)
    if (rho < min | rho > max) stop(paste("rho should be within", min, max)) # positive-definite
    ## build the precision matrix Q
    I <- diag.spam(1, n1*n2)
    Q <- as.spam(prec * (I - rho*W))
    ## Simulate from N(0, Q)
    cholR <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n1 * n2)) ## upper triagle
    X <- backsolve(cholR, rnorm(n1*n2))
    ## return results
    W <- as.matrix(W)
    result <- list(W = W, X = X)
    return(result)
}

#### Simulate CAR data with given rho, precision and weight matrix
CAR.simWmat <- function(rho, prec, W){
    n <- nrow(W)
    I <- diag.spam(1, n)
    Q <- as.spam(prec * (I - rho*W))
    ## Simulate from N(0, Q)
    cholR <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * n)) ## upper triagle
    X <- backsolve(cholR, rnorm(n))
    ## return results
    X
}

#### Simulate samples from a linear model with CAR error
CAR.simLM <- function(pars, data){
    data.vec <- data$data.vec
    rho <- pars[1]
    sigma <- pars[2]
    W <- data$W
    if(length(pars) > 2){
        beta <- pars[-c(1,2)]
        X <-  model.matrix(y~.,data=data.vec)
        Xb <- X %*% beta
        as.numeric(CAR.simWmat(rho, 1/sigma, W) + Xb)
    }else{
        CAR.simWmat(rho, 1/sigma, W)
    }
}

#### Simulate binomial or poisson data with CAR latent varianbles with given spatial weight matrix
CAR.simGLM <- function(method = c("binom", "poisson"), W, n = NULL, pars, Xs = NULL, n.trial = 1){
    rho <- pars[1]
    sigma <- pars[2]
    if(is.null(n)){
        X <- CAR.simWmat(rho, 1/sigma, W)
        Z.car <- list(X=X, W=W)
    }else{
        Z.car <- CAR.simTorus(n[1],n[2],rho,1/sigma)
    }
    if(length(pars) > 2){
        beta <- pars[-c(1,2)]
        eta <- Xs %*% beta + Z.car$X
    }else{
        eta <- Z.car$X
    }
    if(method == "binom"){
        ps <- exp(eta)/(1+exp(eta))
        Emean <- ps
        N <- length(ps)
        Y <- rbinom(N, n.trial, prob = ps)} else{
            lambda <- exp(eta)
            Emean <- lambda
            N <- length(lambda)
            Y <- rpois(N,lambda)
        }
    list(rho = rho, sigma = sigma, beta = beta,
         y = Y, covX = Xs, W = Z.car$W, Z.true = Z.car$X,
         eta=eta, Emean = Emean, n.trial = n.trial)
}
