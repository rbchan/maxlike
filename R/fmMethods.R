
# S3 methods

print.maxlikeFit <- function(x, ...) {
    converge <- x$optim$converge
    cat("\nCall:", paste(deparse(x$call)), "\n\n")
    cat("Coefficients:\n")
    print(x$Est, ...)
    cat("\nAIC:", x$AIC, "\n\n")
    if(converge != 0)
        warning("Model did not converge")
    }



coef.maxlikeFit <- function(object, ...) object$Est[,"Est"]


vcov.maxlikeFit <- function(object, ...) object$vcov





summary.maxlikeFit <- function(object, digits=3, ...) {
    s <- object$Est
    z <- s[,"Est"] / s[,"SE"]
    p <- 2*pnorm(abs(z), lower.tail = FALSE)
    s <- cbind(s, z=z, "P(>|z|)"=p)
    converge <- object$optim$converge
    cat("\nCall:", paste(deparse(object$call)), "\n", fill=TRUE)
    cat("Coefficients:\n")
    print(s, digits=digits, ...)
    cat("\noptim convergence code:", converge, "\n")
    cat("\nAIC:", object$AIC, "\n\n")
    if(converge != 0)
        warning("Model did not converge")
    invisible(s)
    }




logLik.maxlikeFit <- function(object, ...) {
    -object$optim$value
}


AIC.maxlikeFit <- function(object, ..., k=2) {
    2*object$optim$value + k*length(coef(object)[object$not.fixed])
}




predict.maxlikeFit <- function(object, rasters, ...) {
    e <- coef(object)
    if(missing(rasters)) {
        rasters <- try(get(as.character(object$call$rasters)))
        if(identical(class(rasters)[1],  "try-error"))
            stop("could not find the raster data")
        warning("The raster data was not supplied. Using the data found in the workspace.")
    }
    link <- object$link
    cd.names <- layerNames(rasters)
    npix <- prod(dim(rasters)[1:2])
    z <- as.data.frame(matrix(getValues(rasters), npix))
    names(z) <- cd.names
    formula <- object$call$formula
    varnames <- all.vars(formula)
    if(!all(varnames %in% cd.names))
        stop("at least 1 covariate in the formula is not in rasters.")
    Z.mf <- model.frame(formula, z, na.action="na.pass")
    Z <- model.matrix(terms(Z.mf), Z.mf) # Requires R>2.13.0
    eta <- drop(Z %*% coef(object))
    if(identical(link, "logit"))
        psi.hat <- .Call("logit_linkinv", eta, PACKAGE="stats")
    else if(identical(link, "cloglog"))
        psi.hat <- 1-exp(-exp(eta))
    else
        stop("link function should be either 'logit' or 'cloglog'")
    psi.mat <- matrix(psi.hat, dim(rasters)[1], dim(rasters)[2],
                      byrow=TRUE)
    psi.raster <- raster(psi.mat)
    extent(psi.raster) <- extent(rasters)
    psi.raster
}



