maxlike <- function(formula, raster, points, starts, hessian=TRUE, ...) {

    if(identical(formula, ~1))
        stop("At least one predictor variable must be specified in the formula")
    call <- match.call()

    npts <- nrow(points)
    npix <- prod(dim(raster)[1:2])
    varnames <- all.vars(formula)
    layernames <- layerNames(raster)
    if(!all(varnames %in% layernames))
        stop("at least 1 variable in the formula is not in layerNames(raster).")

    cellID <- cellFromXY(raster, points)
    duplicates <- duplicated(cellID)
    uniqueCells <- unique(cellID)
    nptsU <- length(uniqueCells)
    if(nptsU < npts)
        warning(paste("Some cells contain multiple points. \n\tDuplicate points have been discarded.", nptsU, "points were retained"))
    x <- as.data.frame(matrix(extract(raster, uniqueCells), nptsU))
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

    fm <- optim(starts, nll, hessian=hessian, ...)
    par <- fm$par
    vc <- try(solve(fm$hessian))
    if(identical(class(vc), "matrix"))
        se <- sqrt(diag(vc))
    else {
        vc <- matrix(NA, npars, npars)
        se <- rep(NA, npars)
        }
    aic <- 2*fm$value + 2*npars
    out <- list(Est=cbind(Est=par, SE=se), vcov=vc, AIC=aic, call=call,
                retained=points[!duplicates,])
    class(out) <- c("maxlikeFit", "list")
    return(out)
    }




# S3 methods

print.maxlikeFit <- function(x, ...) {
    cat("\nCall:", paste(deparse(x$call)), "\n\n")
    cat("Coefficients:\n")
    print(x$Est, ...)
    cat("\nAIC:", x$AIC, "\n\n")
    }



coef.maxlikeFit <- function(object, ...) object$Est[,"Est"]


vcov.maxlikeFit <- function(object, ...) object$vcov



predict.maxlikeFit <- function(object, newdata, ...) {
    form <- as.formula(object$call$formula)
    if(missing(newdata))
        stop("newdata must be supplied")
    if(!identical(class(newdata)[1], "data.frame"))
        stop("newdata must be a data.frame")
    X <- model.matrix(form, newdata)
    E <- drop(plogis(X %*% coef(object)))
    return(E)
    }



