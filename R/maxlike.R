maxlike <- function(formula, raster, points, starts) {

    x <- as.data.frame(extract(raster, points))
    z <- as.data.frame(getValues(raster))
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
        psix <- plogis(X %*% pars) # observed locations
        psiz <- plogis(Z %*% pars) # all locations
        -1*sum(log(psix/sum(psiz)))
        }

    fm <- optim(starts, nll, hessian=TRUE)
    par <- fm$par
    vc <- solve(fm$hessian)
    se <- sqrt(diag(vc))
    aic <- 2*fm$value + 2*npars
    return(list(Est=cbind(Est=par, SE=se), vcov=vc, AIC=aic))
    }

