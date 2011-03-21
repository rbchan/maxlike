maxlike <- function(formula, raster, points, starts) {

    npts <- nrow(points)
    npix <- prod(dim(raster)[1:2])
    varnames <- all.vars(formula)
    layernames <- layerNames(raster)
    if(!all(varnames %in% layernames))
        stop("names in formula are not in layerNames(raster).")

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
    vc <- solve(fm$hessian)
    se <- sqrt(diag(vc))
    aic <- 2*fm$value + 2*npars
    return(list(Est=cbind(Est=par, SE=se), vcov=vc, AIC=aic))
    }

