maxlike <- function(formula, raster, points, starts) {

    if(identical(formula, ~1))
        stop("At least one predictor variable must be specified in the formula")
    call <- match.call()

    npts <- nrow(points)
    npix <- prod(dim(raster)[1:2])
    varnames <- all.vars(formula)
    layernames <- layerNames(raster)
    if(!all(varnames %in% layernames))
        stop("at least 1 variable in the formula is not in layerNames(raster).")

    x <- as.data.frame(matrix(extract(raster, points), npts))
    z <- as.data.frame(matrix(getValues(raster), npix))
    names(x) <- names(z) <- layernames
    X <- model.matrix(formula, x)
    Z <- model.matrix(formula, z)

    npars <- ncol(X)
    if(missing(starts)) {
        starts <- rep(0, npars)
        names(starts) <- colnames(X)
        }
    else
       names(starts) <- colnames(X)

    nll <- function(pars) {
        psix <- plogis(X %*% pars)
        psiz <- plogis(Z %*% pars)
        -1*sum(log(psix/sum(psiz)))
        }

    fm <- optim(starts, nll, hessian=TRUE)
    par <- fm$par
    vc <- try(solve(fm$hessian))
    if(identical(class(vc), "matrix"))
        se <- sqrt(diag(vc))
    else {
        vc <- matrix(NA, npars, npars)
        se <- rep(NA, npars)
        }
    aic <- 2*fm$value + 2*npars
    out <- list(Est=cbind(Est=par, SE=se), vcov=vc, AIC=aic, call=call)
    class(out) <- c("maxlikeFit", "list")
    return(out)
    }


print.maxlikeFit <- function(fm, ...) {
    cat("\nCall:", paste(deparse(fm$call)), "\n\n")
    cat("Coefficients:\n")
    print(fm$Est, ...)
    cat("\nAIC:", fm$AIC, "\n\n")
    }



coef.maxlikeFit <- function(fm) fm$Est[,"Est"]


vcov.maxlikeFit <- function(fm) fm$vcov


