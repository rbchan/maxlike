maxlike <- function(formula, rasters, points, starts, hessian=TRUE,
                    fixed, removeDuplicates=FALSE, na.action="na.omit", ...)
{
    if(identical(formula, ~1))
        stop("At least one continuous covariate must be specified in the formula")
    varnames <- all.vars(formula)
    call <- match.call()
    npts <- nrow(points)
    cd.class <- class(rasters)[1]
    if(cd.class != "RasterStack")
        stop("rasters must be a raster stack")
    pt.class <- class(points)[1]
    if(!pt.class %in% c("matrix", "data.frame"))
        stop("points must be a matrix or a data.frame")
    if(ncol(points) != 2)
        stop("points must have 2 columns containing the x- and y- coordinates")
    pt.names <- colnames(points)
    if(identical(cd.class, "RasterStack")) {
        cd.names <- layerNames(rasters)
        npix <- prod(dim(rasters)[1:2])
        cellID <- cellFromXY(rasters, points)
        duplicates <- duplicated(cellID)
        if(removeDuplicates) {
            cellID <- unique(cellID)
            npts <- length(cellID)
            }
        x <- as.data.frame(matrix(extract(rasters, cellID), npts))
        z <- as.data.frame(matrix(getValues(rasters), npix))
        names(x) <- names(z) <- cd.names
        }
    if(!all(varnames %in% cd.names))
        stop("at least 1 covariate in the formula is not in rasters.")
    X.mf <- model.frame(formula, x, na.action=na.action)
    X.mf.a <- attributes(X.mf)
    pts.removed <- integer(0)
    if("na.action" %in% names(X.mf.a)) {
        pts.removed <- X.mf.a$na.action
        npts.removed <- length(pts.removed)
        if(npts.removed > 0)
            warning(paste(npts.removed,
                          "points removed due to missing values"))
        }
    X <- model.matrix(formula, X.mf)
    Z.mf <- model.frame(formula, z, na.action=na.action)
    Z.mf.a <- attributes(Z.mf)
    pix.removed <- integer(0)
    if("na.action" %in% names(Z.mf.a)) {
        pix.removed <- Z.mf.a$na.action
        npix.removed <- length(pix.removed)
        if(npix.removed > 0)
            warning(paste(npix.removed,
                          "pixels removed due to missing values"))
        }
    Z <- model.matrix(formula, Z.mf)
    npars <- ncol(X)
    parnames <- colnames(X)
    if(!"(Intercept)" %in% parnames)
        stop("The intercept must be estimated or fixed")
    if(missing(starts)) {
        starts <- rep(0, npars)
        names(starts) <- parnames
        }
    else
       names(starts) <- parnames

    nll <- function(pars) {
        psix <- plogis(X %*% pars)
        psiz <- plogis(Z %*% pars)
        -1*sum(log(psix/sum(psiz) + .Machine$double.xmin))
        }

    is.fixed <- rep(FALSE, npars)
    if(!missing(fixed)) {
        if(length(fixed) != length(starts))
            stop("fixed should be a vector with the same length as the number of parameters to be estimated")
        if(sum(is.real(fixed)) < 1)
            stop("fixed must contain at least one real value")
        is.fixed <- !is.na(fixed)
        if(sum(!is.fixed) < 1)
            stop("you cannot fix all parameters in the model")
        npars <- sum(!is.fixed)
        nll.fix <- function(p) {
            p[is.fixed] <- fixed[is.fixed]
            do.call("nll", list(pars=p))
        }
        fm <- optim(starts, nll.fix, hessian=hessian, ...)
        fm$par[is.fixed] <- fixed[is.fixed]
    } else {
        fm <- optim(starts, nll, hessian=hessian, ...)
    }
    not.fixed <- !is.fixed

    par <- fm$par
    if(hessian) {
        vcTry <- try(solve(fm$hessian[not.fixed, not.fixed]))
        if(identical(class(vcTry), "matrix")) {
            vc <- matrix(0, length(par), length(par))
            vc[not.fixed, not.fixed] <- vcTry
            se <- sqrt(diag(vc))
        }
        else {
            vc <- matrix(NA, npars, npars)
            se <- rep(NA, npars)
        }
    } else {
        vc <- matrix(NA, npars, npars)
        se <- rep(NA, npars)
    }
    dimnames(vc) <- list(parnames, parnames)
    aic <- 2*fm$value + 2*npars

    fitted <- plogis(Z %*% par)
    out <- list(Est=cbind(Est=par, SE=se), vcov=vc, AIC=aic, call=call,
                pts.removed=pts.removed, pix.removed=pix.removed,
#                duplicates=duplicates,
                optim=fm, not.fixed=not.fixed)
    class(out) <- c("maxlikeFit", "list")
    return(out)
    }



