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


predict.maxlikeFit <- function(object, newdata=NULL, ...) {
    e <- coef(object)
 
# allowing for data.frames (Roeland Kindt)
    if (is.null(newdata) == F) {
# allowing for new RasterStack (Roeland Kindt, Nov 2019)
        if ("RasterStack" %in% class(newdata)) {
            rasters <- newdata
            cd.names <- names(newdata)
            npix <- prod(dim(newdata)[1:2])
            z <- as.data.frame(matrix(getValues(newdata), npix))
            names(z) <- cd.names
        }else{
            z <- newdata
        }
        cd.names <- names(z)
    }else{
        rasters <- object$rasters
        if(is.null(rasters)) {
            rasters <- try(get(as.character(object$call$rasters)))
            if(identical(class(rasters)[1],  "try-error"))
                stop("could not find the raster data")
            warning("raster data were not saved with object, using the data found in the workspace instead.")
        }
#       link <- object$link
        cd.names <- names(rasters)
        npix <- prod(dim(rasters)[1:2])
        z <- as.data.frame(matrix(getValues(rasters), npix))
        names(z) <- cd.names
    }
    formula <- object$call$formula
    varnames <- all.vars(formula)
# reversed cd.names and varnames to allow for polynomials etc via I(), Roeland Kindt
# re-reversed...
#    if(!all(cd.names %in% varnames)) {
    if(!all(varnames %in% cd.names)) {
        if (is.null(newdata) == F) {
            stop("at least 1 covariate in the formula is not in new data")
        }else{
            stop("at least 1 covariate in the formula is not in rasters")
        }
    }
    Z.mf <- model.frame(formula, z, na.action="na.pass")
    Z.terms <- attr(Z.mf, "terms")
    Z <- model.matrix(Z.terms, Z.mf)
    eta <- drop(Z %*% coef(object))
    link <- object$link
    if(identical(link, "logit"))
        psi.hat <- plogis(eta)
    else if(identical(link, "cloglog"))
        psi.hat <- 1-exp(-exp(eta))
    else
        stop("link function should be either 'logit' or 'cloglog'")
# Corrected Roeland Kindt, Nov 2019
    if (is.null(newdata) == F  && ("RasterStack" %in% class(newdata)) == F) {    
        return(psi.hat)
    }else{
        psi.mat <- matrix(psi.hat, dim(rasters)[1], dim(rasters)[2], byrow=TRUE)
        psi.raster <- raster(psi.mat)
        extent(psi.raster) <- extent(rasters)
        return(psi.raster)
    }
}
 

