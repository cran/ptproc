linear.cond.int <- function(params, eval.pts, pts = NA, data = NULL, TT = NULL) {
    mu <- params[1]
    beta <- params[-1]

    if(is.null(TT)) {
        ## Evaluate
        ci <- mu + eval.pts %*% beta
        ci <- as.vector(ci)
    }
    else {
        ## Integrate
        total.vol <- prod(apply(TT, 2, diff))
        m.vol <- sapply(1:ncol(TT), function(i)
                    {
                        z <- TT[,-i,drop=FALSE]
                        prod(apply(z, 2, diff))
                    })
        d <- apply(TT^2 / 2, 2, diff)
        ci <- mu * total.vol + (beta * d) %*% m.vol
    }
    ci
}

hawkes.cond.int <- function(params, eval.pts, pts = NA, data = NULL, TT = NULL) {
    mu <- params[1]
    C <- params[2]    
    ak <- params[-(1:2)]
    K <- length(ak)

    if(K < 1) 
        stop("K must be >= 1")
    if(is.null(TT)) {
        S <- sapply(as.vector(eval.pts), function(x, times, ak, C)
                {
                    use <- times < x
                    if(!is.na(use) && any(use)) {
                        d <- x - times[use]
                        k <- 0:(length(ak)-1)
                        lxk <- outer(log(d), k) + (-C * d)
                        sum(exp(lxk) %*% ak)
                    }
                    else 0
                }, times = as.vector(pts), ak = ak, C = C)
        ci <- mu + S
    }
    else {
        Rfunc <- function(x, L, c) {
            k <- 0:(L-1)
            g <- gamma(k+1) / c^(k+1)
            
            ## Evaluate the (incomplete) Gamma function for all x and k
            o <- outer(x, k+1, pgamma, scale = 1/c)
            r <- t(t(o) * g)
            r
        }
        times <- as.vector(pts)
        ci <- mu*(TT[2,1]-TT[1,1]) 
        S <- double(2)
        
        for(i in 1:2) {
            use <- times < TT[i,1]

            if(any(use)) {
                r <- Rfunc(TT[i,1] - times[use], K, C
                           )
                S[i] <- sum(r %*% ak)
            }
        }
        ci <- ci + (S[2] - S[1])
    }
    ci
}

hPois.cond.int <- function(params, eval.pts, pts = NA, data = NULL, TT = NULL) {
    mu <- params[1]

    if(is.null(TT))
        rep(mu, nrow(eval.pts))
    else {
        vol <- prod(apply(TT, 2, diff))
        mu * vol
    }
}
