ptproc <- function(pts = NA, cond.int, params, fixed.params = rep(NA, length(params)),
                   condition = NULL, ranges = NULL, data = NULL, ndim = NULL, is.pois = FALSE) {
    if(is.na(pts)) {
        cat("Point process model for simulation only\n")

        if(is.null(ranges))
            stop("`ranges' needs to be specified")
        if(is.null(ndim))
            stop("`ndim' needs to be specified")
    }
    else {
        if(is.vector(pts) || is.data.frame(pts))
            pts <- as.matrix(pts)
        ndim <- NCOL(pts)
        
        if(is.null(ranges)) 
            ranges <- apply(pts, 2, range)        
    }
    if(is.character(cond.int)) {
        stop("No built-in CIF's right now")
    }
    cif <- match.fun(cond.int)
    fun.name <- toupper(sub(".cond.int", "", deparse(substitute(cond.int))))

    if(length(params) != length(fixed.params))
        stop("length(params) != length(fixed.params)")
    if(is.null(names(params)))
        names(params) <- paste("p", 1:length(params), sep ="")
    names(fixed.params) <- names(params)
    ppm <- list(pts = pts, cond.int = cif, params = params,
                fixed.params = fixed.params, ranges = ranges,
                model = ifelse(is.character(cond.int), toupper(cond.int), fun.name),
                condition = condition, data = data, ndim = ndim, is.pois = is.pois)
    class(ppm) <- "ptproc"
    invisible(ppm)
}


residuals.ptproc <- function(object, type = c("ordinary", "approx"),
                             m = NULL, K = NULL, ...) {
    ppobj <- object
    type <- match.arg(type)
    pts <- ppobj$pts
    
    if(is.na(pts))
        stop("Cannot construct residual process when `pts' is `NA'")
    ci <- evalCIF(ppobj)

    if(type == "ordinary") {
        if(is.null(m))
            stop("`m' must be specified for ordinary thinned residuals")
        if(any(ci < m))
            stop("Some conditional intensity values < m")
        p <- m / ci
        coin <- rbinom(length(ci), 1, p)
        keep <- pts[coin == 1, , drop = FALSE]
        attr(keep, "rate") <- m
    }
    else {
        if(is.null(K))
            stop("`K' must be specified for approximate thinned residuals")
        p <- (1 / ci) / sum(1 / ci)
        idx <- 1:nrow(pts)
        r.idx <- sample(idx, K, rep = FALSE, prob = p)
        r.idx <- sort(r.idx)
        keep <- pts[r.idx, , drop = FALSE]
        attr(keep, "rate") <- K / prod(apply(object$ranges, 2, diff))
    }
    attr(keep, "type") <- type
    colnames(keep) <- colnames(ppobj$pts)
    keep
}

make.box <- function(ppobj, idx = 1:2) {
    ranges <- ppobj$ranges
    a <- c(ranges[,idx[1]], rev(ranges[,idx[1]]))
    b <- c(rep(ranges[1,idx[2]],2),rep(ranges[2,idx[2]],2))
    cbind(a,b)
}

evalCIF <- function(ppobj, xpts = ppobj$pts, ...) {
    if(!inherits(ppobj, "ptproc")) 
        stop("Only use with ptproc objects!\n")
    if(!is.matrix(xpts))
        xpts <- as.matrix(xpts)
    cif <- ppobj$cond.int    
    ci <- cif(params = params(ppobj), eval.pts = xpts, pts = ppobj$pts,
              data = ppobj$data, ...)
    ci
}

integrateCIF <- function(ppobj, TT = ppobj$ranges, ...) {
    if(!inherits(ppobj, "ptproc")) 
        stop("Only use with ptproc objects!\n")
    if(missing(TT))
        TT <- ppobj$ranges
    cif <- ppobj$cond.int
    int <- cif(params = params(ppobj), eval.pts = NULL, pts = ppobj$pts,
               data = ppobj$data, TT = TT, ...)
    as.numeric(int)
}

logLik.ptproc <- function(object, negative = FALSE, ...) {
    ci <- evalCIF(object)

    if(any(ci < 0)) {
        stop("Some conditional intensity values < 0")
    }
    L1 <- sum(log(ci))
    L2 <- integrateCIF(object)
    LL <- L1 - L2
    if(negative)
        LL <- -LL    
    attr(LL, "df") <- (length(params(object)) - length(na.omit(fixed(object))))
    class(LL) <- "logLik"
    LL
}

ptproc.fit <- function(ppobj, optim.control = list(), method = "Nelder-Mead",
                       alpha = 0, ...) {
    if(!inherits(ppobj, "ptproc")) 
        stop("Only use with `ptproc' objects!")
    ppobj$initial.params <- ppobj$params
    names(ppobj$initial.params) <- names(ppobj$params)
    fixed <- ppobj$fixed.params
    mask <- is.na(fixed)
    pstart <- ppobj$params[mask]
    optim.logLik <- make.optim.logLik(ppobj, alpha)

    fit <- optim(pstart, optim.logLik, method = method, control = optim.control, ...)

    ##if(fit$convergence > 0) 
    ##    warning(paste("optim returned code =", fit$convergence))
    ppobj$params[mask] <- fit$par
    ppobj$hessian <- fit$hessian
    ppobj$convergence <- fit$convergence
    ppobj
}

make.optim.logLik <- function(ppobj, alpha) {
    optim.logLik <- function(params) {
        fixed <- ppobj$fixed.params
        mask <- is.na(fixed)
        ppobj$params[mask] <- params
        condition <- ppobj$condition
        returnflag <- NULL

        if(!is.null(condition)) {
            if(!is.expression(condition))
                warning("condition must be an expression")
            else
                eval(condition)    
        }
        if(!is.null(returnflag))
            return(returnflag)
        LL <- logLik(ppobj, negative = TRUE)
    }
}


ptproc.sim <- function(ppobj, M) {
    if(!inherits(ppobj, "ptproc")) 
        stop("Only use with ptproc objects!")
    ranges <- ppobj$ranges
    n <- ceiling(M * prod(apply(ranges, 2, diff)))

    ## Simulate homogeneous Poisson with rate M
    u <- sapply(1:ppobj$ndim, function(i) runif(n,ranges[1,i],ranges[2,i]))

    if(is.null(dim(u)))
        u <- matrix(u, byrow = TRUE, ncol = ppobj$ndim)

    ## Order by the first column (in case they are times)
    u <- u[order(u[,1]), , drop = FALSE]

    if(ppobj$is.pois) {
        ci <- evalCIF(ppobj, xpts = u)

        if(any(ci > M))
            stop("Some conditional intensity values > M")    
        coin <- rbinom(n, 1, ci / M)
        u.keep <- u[coin > 0, , drop = FALSE]
    }
    else {
        u.keep <- NULL
        sim.obj <- ppobj
        sim.obj$pts <- NA
        
        for(i in 1:n) {
            ci <- evalCIF(sim.obj, xpts = u[i, , drop = FALSE])
            
            if(ci > M) {
                cat("Conditional intensity value > M;  Stopping.\n")
                return(u.keep)
            }
            coin <- rbinom(1, 1, ci / M)
            
            if(coin > 0) {
                u.keep <- rbind(u.keep, u[i, ])
                
                if(is.na(sim.obj$pts))
                    sim.obj$pts <- u[i, , drop = FALSE]
                else
                    sim.obj$pts <- rbind(sim.obj$pts, u[i, ])
            }    
        }
    }
    u.keep
}


penalty <- function(code = NULL, condition = FALSE) {
    if(!is.character(condition))
        stop("`condition' should be a character string")
    if(!is.null(code) && !is.character(code))
        stop("`code' should be a character string")
    condition <- paste("if(", condition, ") returnflag <- alpha", sep = "")

    if(!is.null(code)) 
        condition <- paste(code, condition, sep = ";")
    parse(text = condition)
}



## Some methods

fixed <- function(x) {
    if(!inherits(x, "ptproc"))
        stop("Only use with `ptproc' objects!")
    x$fixed.params
}

"fixed<-" <- function(x, value) {
    if(!inherits(x, "ptproc"))
        stop("Only use with `ptproc' objects!")
    if(length(value) == 1)
        x$fixed.params <- rep(value, length(coef(x)))
    else        
        x$fixed.params <- value
    x
}

params <- function(x) {
    if(!inherits(x, "ptproc"))
        stop("Only use with `ptproc' objects!")
    x$params
}

"params<-" <- function(x, value) {
    if(!inherits(x, "ptproc"))
        stop("Only use with `ptproc' objects!")
    x$params <- value
    x
}

condition <- function(x) {
    if(!inherits(x, "ptproc"))
        stop("Only use with `ptproc' objects!")
    x$condition
}

"condition<-" <- function(x, value) {
    if(!inherits(x, "ptproc"))
        stop("Only use with `ptproc' objects!")
    if(!is.expression(value))
        stop("Right hand side needs to be an R expression")
    x$condition <- value
    x
}

print.ptproc <- function(x, digits = getOption("digits") - 3, ...) {
    cat("Model type:", x$model, "\n\n")
    cat("Parameter Values:\n")
    print.default(format(params(x), digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\nInitial Values:\n")
    print.default(format(x$initial.params, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\nFixed Parameters:\n")
    fp <- as.numeric(x$fixed.params)
    names(fp) <- names(x$fixed.params)
    print.default(fp)

    if(!is.null(x$condition)) {
        e <- paste(deparse(x$condition, width = 20)[1], "...", sep = "")
        cat("\nCondition:", e, "\n")
    }
}

summary.ptproc <- function(object, ...) {
    npts <- NROW(object$pts)
    vol <- prod(apply(object$ranges, 2, diff))
    poisson.AIC <- -2 * (npts * log(npts / vol) - npts) + 2
    model.AIC <- AIC(object)

    s <- list(poisson.AIC = poisson.AIC, model.AIC = model.AIC, ppobj = object)
    class(s) <- "summary.ptproc"
    s
}

print.summary.ptproc <- function(x, digits = getOption("digits") - 3, ...) {
    cat("Model type:", x$ppobj$model, "\n\n")
    cat("Parameter Values:\n")
    print.default(format(params(x$ppobj), digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\nModel AIC:\t", x$model.AIC)
    cat("\nH. Pois. AIC:\t", x$poisson.AIC)
    cat("\n")
}



log.surv <- function(res, theoretical = TRUE, xlab = "Interevent time",
                     ylab = "Cumulative number", ...) {
    u <- diff(res)
    plot(sort(u), length(u):1, log="y", xlab = xlab, ylab = ylab,
         main="Log Survivor Plot\nof Interevent Times", ...)
    
    if(theoretical) 
        lines(sort(u), length(u) * (1 - pexp(sort(u), attr(res, "rate"))))
}


stationarity <- function(res, h, std = TRUE, xlab = "Transformed time",
                         ylab="Standardized # of events per interval",
                         type = "b", ...) {
    nbins <- ceiling((max(res) - min(res)) / h)
    breaks <- seq(min(res), max(res), len = nbins)
    b <- hist(res, breaks, plot = FALSE)$counts
    std.hi <- 1:3
    std.lo <- -std.hi
    m <- h * attr(res, "rate")
    scaled.freq <- (b - m) / sqrt(m)
    plot(breaks[-length(breaks)], scaled.freq,
         ylim = range(c(std.hi, std.lo, scaled.freq)),
         xlab = xlab, ylab = ylab,
         main = paste("Stationarity Plot\nh = ", h, sep=""),
         type = type, ...)
    if(std)
        abline(h=c(std.hi,std.lo), lty=3)
}


        
