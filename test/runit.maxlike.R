test.maxlike.fit.simple.1 <- function() {

    data(MaungaWhau)
    elev <- raster(MaungaWhau$elev, 0, 61, 0, 87)
    precip <- raster(MaungaWhau$precip, 0, 61, 0, 87)
    xy <- MaungaWhau$xy

    # Stack them and make sure they are named
    ep <- stack(elev, precip)
    layerNames(ep) <- c("elev", "precip")

    # Fit a model
    fm <- maxlike(~elev + I(elev^2) + precip, ep, xy)

    # Check estimates
    checkEqualsNumeric(coef(fm),
                       c(0.5366934, 2.4465578, -2.3575862, 2.1310296),
                       tol=1e-6)

    # Check variance-covariance matrix
    checkEqualsNumeric(vcov(fm), matrix(c(
                       0.05204765,  0.03300724, -0.03589617,  0.03561091,
                       0.03300724,  0.03921600, -0.03373490,  0.03307982,
                       -0.03589617, -0.03373490, 0.03618861, -0.03398646,
                       0.03561091,  0.03307982, -0.03398646,  0.04569138),
                       4, 4, byrow=TRUE), tol=1e-6)

    # Add missing values and refit
    elev2 <- elev
    elev2[c(1,5)] <- NA
    xy2 <- xy
    xy2[2,] <- NA
    ep2 <- stack(elev2, precip)
    layerNames(ep2) <- c("elev", "precip")
    fm2 <- maxlike(~elev + I(elev^2) + precip, ep2, xy2)
    checkEqualsNumeric(fm2$pix.removed, c(1,5))
    checkEqualsNumeric(fm2$pts.removed, 2)

}




