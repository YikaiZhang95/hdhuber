#' Cross-validation for Selecting the Tuning Parameter in Penalized Huber Regression
#'
#' Performs k-fold cross-validation for \code{\link{hdhuber}}.
#'
#' @param x A numerical matrix with \eqn{n} rows (observations) and \eqn{p} columns (variables).
#' @param y Response variable.
#' @param lambda Optional; a user-supplied sequence of \code{lambda} values. If \code{NULL}, 
#'   \code{\link{hdhuber}} selects its own sequence.
#' @param nfolds Number of folds for cross-validation. Defaults to 5.
#' @param foldid Optional vector specifying the indices of observations in each fold.
#'   If provided, it overrides \code{nfolds}.
#'@param delta (\code{delta}) used in the loss function.
#' @param ... Additional arguments passed to \code{\link{hdhuber}}.
#'
#' @details
#' This function computes the average cross-validation error and provides the standard error.
#'
#' @return
#' An object with S3 class \code{cv.hdhuber} consisting of
#'   \item{lambda}{Candidate \code{lambda} values.}
#'   \item{cvm}{Mean cross-validation error.}
#'   \item{cvsd}{Standard error of the mean cross-validation error.}
#'   \item{cvup}{Upper confidence curve: \code{cvm} + \code{cvsd}.}
#'   \item{cvlo}{Lower confidence curve: \code{cvm} - \code{cvsd}.}
#'   \item{lambda.min}{\code{lambda} achieving the minimum cross-validation error.}
#'   \item{lambda.1se}{Largest \code{lambda} within one standard error of the minimum error.}
#'   \item{cv.min}{Cross-validation error at \code{lambda.min}.}
#'   \item{cv.1se}{Cross-validation error at \code{lambda.1se}.}
#'   \item{hdhuber.fit}{a fitted \code{\link{hdhuber}} object for the full data.}
#'   \item{nzero}{Number of non-zero coefficients at each \code{lambda}.}
#' @keywords models regression
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' delta = 0.5
#' cv.fit <- cv.hdhuber(x = x, y = y, delta = delta)
#' @export
#'


cv.hdhuber <- function(x, y, lambda=NULL, delta, nfolds=5L, foldid, ...) {
  ####################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row) 
    stop("x and y have different number of observations.")  
  ####################################################################
  hdhuber.object <- hdhuber(x, y, lambda=lambda, delta=delta, ...)
  lambda <- hdhuber.object$lambda
  nz <- sapply(coef(hdhuber.object, type="nonzero"), length) 
  if (missing(foldid)) 
    foldid <- sample(rep(seq(nfolds), 
                         length=x.row)) else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=5 recommended.")
  outlist <- as.list(seq(nfolds))
  ## fit the model nfold times and save them
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- hdhuber(x=x[!which, , drop=FALSE], 
                         y=y[!which], delta=delta, lambda=lambda, ...)
  }
  ## select the lambda according to predmat
  cvstuff <- cvpath.hdhuber(outlist, x, y, delta, lambda, 
                         foldid, x.row, ...)
  
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlower = cvm - cvsd,
              nzero = nz, name = cvname,
              hdhuber.fit = hdhuber.object)
  obj = c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) = "cv.hdhuber"
  obj
} 

cvpath.hdhuber <- function(outlist, x, y, delta, lambda, foldid, x.row, ...) {
  nfolds <- max(foldid)
  predmat <- matrix(NA, x.row, length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[whichfold, , drop = FALSE])
    nlami <- length(fitobj$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- huber_loss(y-predmat, delta)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- colMeans(cvraw, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(colMeans(scaled, na.rm = TRUE) / (N - 1))
  out <- list(cvm=cvm, cvsd=cvsd, cvraw=cvraw)
  out
} 