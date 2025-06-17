
#' @useDynLib myglmnet, .registration = TRUE
#' @export
myglmnet <- function(x, y, family = c("gaussian", "logistic", "binomial", "cox", "quantile"), weights = NULL, offset = NULL,
                     alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs <= nvars, 0.01, 1e-04),
                     lambda = NULL, standardize = TRUE, intercept = TRUE, thresh = 1e-07, dfmax = nvars +
                     1, pmax = min(dfmax * 2 + 20, nvars), exclude = NULL, penalty.factor = rep(1,
                     nvars), lower.limits = -Inf, upper.limits = Inf, maxit = 1e+04, beta0=NULL, h=0.25, tau=0.5) {

    this.call <- match.call()
    ### Need to do this first so defaults in call can be satisfied
    np <- dim(x)
    ## check dims
    if (is.null(np) | (np[2] <= 1))
        stop("x should be a matrix with 2 or more columns")
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    family <- match.arg(family)
    if(family == "binomial"){
        family <- "logistic"
    }
    if (alpha > 1) {
        warning("alpha >1; set to 1")
        alpha <- 1
    }
    if (alpha < 0) {
        warning("alpha<0; set to 0")
        alpha <- 0
    }
    alpha <- as.double(alpha)
    nlam <- as.integer(nlambda)
    y <- drop(y)  # we dont like matrix responses unless we need them
    if (is.null(weights))
        weights <- rep(1/nobs, nobs) else if (length(weights) != nobs)
        stop(paste("number of elements in weights (", length(weights), ") not equal to the number of rows of x (",
            nobs, ")", sep = "")) else weights <- weights/(sum(weights))
    dimy <- dim(y)
    nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
    if (nrowy != nobs)
        stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (",
            nobs, ")", sep = ""))
    if(inherits(x, "PlinkMatrix")){
        vnames <- x@colname
    } else{
        vnames <- colnames(x)
    }
    if (is.null(vnames))
        vnames <- paste("V", seq(nvars), sep = "")
    ne <- as.integer(dfmax)
    nx <- as.integer(pmax)
    if (is.null(exclude))
        exclude <- integer(0)
    if (any(penalty.factor == Inf)) {
        exclude <- c(exclude, seq(nvars)[penalty.factor == Inf])
        exclude <- sort(unique(exclude))
    }
    ju <- rep(1L, nvars)
    if (length(exclude) > 0) {
        jd <- match(exclude, seq(nvars), 0)
        if (!all(jd > 0))
            stop("Some excluded variables out of range")
        penalty.factor[jd] <- 1  #ow can change lambda sequence
        ju[jd] <- 0L
    }
    vp <- as.double(penalty.factor)
    internal.parms <- glmnet.control()

    # if (internal.parms$itrace) trace.it <- 1 else { if (trace.it) {
    # glmnet.control(itrace = 1) on.exit(glmnet.control(itrace = 0)) } } check on
    # limits
    if (any(lower.limits > 0)) {
        stop("Lower limits should be non-positive")
    }
    if (any(upper.limits < 0)) {
        stop("Upper limits should be non-negative")
    }
    lower.limits[lower.limits == -Inf] <- -internal.parms$big
    upper.limits[upper.limits == Inf] <- internal.parms$big
    if (length(lower.limits) < nvars) {
        if (length(lower.limits) == 1)
            lower.limits <- rep(lower.limits, nvars) else stop("Require length 1 or nvars lower.limits")
    } else lower.limits <- lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars) {
        if (length(upper.limits) == 1)
            upper.limits <- rep(upper.limits, nvars) else stop("Require length 1 or nvars upper.limits")
    } else upper.limits <- upper.limits[seq(nvars)]
    cl <- rbind(lower.limits, upper.limits)
    if (any(cl == 0)) {
        ### Bounds of zero can mess with the lambda sequence and fdev; ie nothing happens
        ### and if fdev is not zero, the path can stop fdev <- glmnet.control()$fdev if
        ### (fdev != 0) { glmnet.control(fdev = 0) on.exit(glmnet.control(fdev = fdev)) }
        stop("This case has not been implemented")
    }
    storage.mode(cl) <- "double"
    ### end check on limits

    isd <- as.integer(standardize)

    if(family == "cox") {
        if(!missing(intercept)){
            warning("Cox model has no intercept")
        }
        intercept = FALSE
        if(!survival::is.Surv(y)){
            stop("y has to be a Surv object for Cox model")
        }
        # y[,1] are the time and y[,2] are the status
        ytime = y[,1]
        yorder = order(ytime)
        ytime = ytime[yorder]
        rankmin = rank(ytime, ties.method="min") - 1L
        rankmax = rank(ytime, ties.method="max") - 1L
        yorder = yorder - 1L

        weights = weights * y[, 2]
        y = list(yorder, rankmin, rankmax)
    }

    intr <- as.integer(intercept)

    # Don't have this yet jsd <- as.integer(standardize.response)
    thresh <- as.double(thresh)
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1)
            stop("lambda.min.ratio should be less than 1")
        flmin <- as.double(lambda.min.ratio)
        ulam <- double(1)
    } else {
        flmin <- as.double(1)
        if (any(lambda < 0))
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    is.sparse <- FALSE
    ix <- jx <- NULL
    if (inherits(x, "sparseMatrix")) {
        stop("sparse matrices not implemented yet")
        ## Sparse case
        is.sparse <- TRUE
        x <- as(x, "CsparseMatrix")
        x <- as(x, "dgCMatrix")
        ix <- as.integer(x@p + 1)
        jx <- as.integer(x@i + 1)
        x <- as.double(x@x)
    }

    maxit <- as.integer(maxit)
    lmu <- integer(1)
    a0 <- double(nlam)
    ca <- double(nx * nlam)
    ia <- integer(nx)
    nin <- integer(nlam)
    devratio <- double(nlam)
    alm <- double(nlam)
    nlp <- integer(1)
    jerr <- integer(1)
    nulldev <- double(1)
    # compute weighted residuals for each lambda
    # save it in this matrix
    residuals <- double(nlam * nobs)

    if (is.null(offset)) {
        offset <- double(0)
        has_offset <- 0L
    } else if (length(offset) != nobs) {
        stop("The length of offset must be the same as the number of observations")
    } else {
        has_offset <- 1L
    }
    mxitnr <- internal.parms$mxitnr
    alpha_mod <- alpha
    if (family == "gaussian") {
        mxitnr <- 1
        sdy <- sd(y)
        y <- y/sdy
        cl <- cl/sdy
        ulam <- ulam * (1 - alpha + alpha/sdy)
        alpha_mod <- alpha / (sdy * (1- alpha)+ alpha) #modified alpha, keep original
        if(!is.null(beta0)){
            beta0 <- beta0/sdy
        }
    }
    mxitnr <- as.integer(mxitnr)

    if(inherits(x, "PlinkMatrix")){
        if(x@ncov > 0){
            covmat <- x@covs
            if(intercept) {
                # Maybe the mean and sd should be computed using the
                # observation weights?
                covmean <- apply(covmat, 2, mean)
                covmat <- sweep(covmat, 2, covmean)
            }
            if(standardize) {
                covsd <- apply(covmat, 2, sd)
                covmat <- sweep(covmat, 2, covsd, "/")
                if(!is.null(beta0)){
                    beta0[1:x@ncov] <- beta0[1:x@ncov] * covsd
                }
            }
            x_list = list("Plink", np[1], np[2], x@ptr, x@xim, x@ncov, covmat)
        } else {
             x_list = list("Plink", np[1], np[2], x@ptr, x@xim, 0L)
        }

    } else {
        x_list = list("Dense", np[1], np[2], x)
    }

    # warm start support
    iy = integer(np[2])
    mm = integer(np[2])
    nino = integer(1)
    warm=0L # warm start flag
    if(!is.null(beta0)){
        if(length(beta0) != np[2]) {
            stop("warm start dimension does not match matrix dimenstion")
        }
        if(is.null(lambda)) {
            stop("warm start must be given regularization parameters")
        }
        nnz_ind = which(beta0 != 0)
        nino[1] = length(nnz_ind)
        if(nino > nx) {
            stop("Already have too many variables in the model")
        }
        if(nino > 0){
            iy[nnz_ind] = 1L
            mm[nnz_ind] = 1:nino
            ia[1:nino] = nnz_ind - 1L # 0 based index
        }
        beta0 = as.double(beta0)
        warm = 1L
    } else {
        beta0 = double(np[2])
    }
    .Call("solve", alpha_mod, x_list, y, weights, ju, vp, cl, nx, nlam, flmin, ulam, thresh,
        isd, intr, maxit, lmu, a0, ca, ia, nin, devratio, alm, nlp, family, offset,
        has_offset, mxitnr, nulldev, jerr, beta0, iy, mm, nino, warm, residuals, h, tau)

    # if (trace.it) { if (relax) cat('Training Fit\n') pb <- createPB(min = 0, max =
    # nlam, initial = 0, style = 3) }

    # adjust the result back to original scale
    if (family == "gaussian") {
        ca <- ca * sdy
        a0 <- a0 * sdy
        alm <- alm/(1 - alpha + alpha/sdy)
        residuals <- residuals * sdy
    }
    fit <- list(jerr = jerr, a0 = a0, nin = nin, lmu = lmu, alm = alm, ca = ca, ia = ia)
    outlist <- getcoef(fit, nvars, nx, vnames)
    dev <- devratio[seq(lmu)]
    outlist <- c(outlist, list(dev.ratio = dev, nulldev = nulldev, npasses = nlp,
        jerr = jerr, offset = has_offset))

    if(inherits(x, "PlinkMatrix") && x@ncov > 0) {
        beta_tmp = outlist$beta[1:x@ncov, , drop=F]
        if(standardize){
            beta_tmp = Matrix::Diagonal(x@ncov,1/covsd) %*% beta_tmp
            outlist$beta[1:x@ncov, ] = beta_tmp
        }
        if(intercept){
            outlist$a0 = outlist$a0 - Matrix::t(covmean %*% beta_tmp)
        }
    }

    residuals <- residuals[seq(lmu*nobs)]
    residuals <- matrix(residuals, nrow = nobs, ncol = lmu)
    outlist$residuals <- residuals
    outlist$call <- this.call
    outlist$nobs <- nobs
    class(outlist) <- c(class(outlist), "glmnet")
    outlist
}

#' @export
PlinkPredict = function(myfit, x){
    p = nrow(myfit$beta)
    dense_beta = matrix(myfit$beta, nrow=p)
    nlam = ncol(dense_beta)
    n = nrow(x)
    if(ncol(x) != p) {
        stop("incorrect x dimension")
    }
    result = double(n*nlam)
    intercepts = as.numeric(myfit$a0)
    .Call("PlinkMultvC", x@ptr, x@xim, x@covs, n, p , x@ncov, dense_beta, intercepts, result)
    result = matrix(result, nrow=n)
    return(result)
}
