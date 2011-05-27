maxlike <- function(formula, covData, ptData, removeDuplicates=FALSE,
                    starts, hessian=TRUE, na.action="na.omit", ...)
{
    if(identical(formula, ~1))
        stop("At least one covariate must be specified in the formula")
    varnames <- all.vars(formula)
    call <- match.call()
    npts <- nrow(ptData)
    cd.class <- class(covData)[1]
    if(!cd.class %in% c("RasterStack", "data.frame"))
        stop("covData must be a raster stack or a data.frame")
    pt.class <- class(ptData)[1]
    if(!pt.class %in% c("matrix", "data.frame"))
        stop("ptData must be a matrix or a data.frame")
    pt.names <- colnames(ptData)
    if(identical(cd.class, "data.frame")) {
        cd.names <- colnames(covData)
        npix <- nrow(covData)
        if(!all(cd.names %in% pt.names))
            stop("When covData is a data.frame, ptData must be a data.frame with the same column names")
        if(removeDuplicates)
            warning("removeDuplicates ignored when covData is a data.frame")
        duplicates <- rep(FALSE, npts)
        x <- as.data.frame(ptData)
        z <- covData
        }
    if(identical(cd.class, "RasterStack")) {
        cd.names <- layerNames(covData)
        npix <- prod(dim(covData)[1:2])
        cellID <- cellFromXY(covData, ptData)
        duplicates <- duplicated(cellID)
        if(removeDuplicates) {
            cellID <- unique(cellID)
            npts <- length(cellID)
            }
        x <- as.data.frame(matrix(extract(covData, cellID), npts))
        z <- as.data.frame(matrix(getValues(covData), npix))
        names(x) <- names(z) <- cd.names
        }
    if(!all(varnames %in% cd.names))
        stop("at least 1 covariate in the formula is not in covData.")
    X.mf <- model.frame(formula, x, na.action=na.action)
    X.mf.a <- attributes(X.mf)
    pts.removed <- integer(0)
    if("na.action" %in% names(X.mf.a)) {
        pts.removed <- X.mf.a$na.action
        npts.removed <- length(pts.removed)
        if(npts.removed > 0)
            warning(paste(npts.removed, "points removed due to missing values"))
        }
    X <- model.matrix(formula, X.mf)
    Z.mf <- model.frame(formula, z, na.action=na.action)
    Z.mf.a <- attributes(Z.mf)
    pix.removed <- integer(0)
    if("na.action" %in% names(Z.mf.a)) {
        pix.removed <- Z.mf.a$na.action
        npix.removed <- length(pix.removed)
        if(npix.removed > 0)
            warning(paste(npix.removed, "pixels removed due to missing values"))
        }
    Z <- model.matrix(formula, Z.mf)
    npars <- ncol(X)
    parnames <- colnames(X)
    if(!"(Intercept)" %in% parnames)
        stop("The intercept must be estimated")
    if(missing(starts)) {
        starts <- rep(0, npars)
        names(starts) <- parnames
        }
    else
       names(starts) <- parnames

    nll <- function(pars) {
        psix <- plogis(X %*% pars)
        psiz <- plogis(Z %*% pars)
        -1*sum(log(psix/sum(psiz)))
        }

    fm <- optim(starts, nll, hessian=hessian, ...)
    par <- fm$par
    if(hessian) {
        vc <- try(solve(fm$hessian))
        if(identical(class(vc), "matrix"))
           se <- sqrt(diag(vc))
        }
    else {
        vc <- matrix(NA, npars, npars)
        se <- rep(NA, npars)
        }
    aic <- 2*fm$value + 2*npars
    out <- list(Est=cbind(Est=par, SE=se), vcov=vc, AIC=aic, call=call,
                retained=ptData[!duplicates,], pts.removed=pts.removed,
                pix.removed=pix.removed)
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



