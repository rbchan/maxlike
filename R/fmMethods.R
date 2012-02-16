
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






