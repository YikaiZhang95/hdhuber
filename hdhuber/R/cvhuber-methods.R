#' Extract Coefficients from a `cv.hdhuber` Object
#'
#' Retrieves coefficients from a cross-validated `hdhuber()` model, using the
#' stored `"hdhuber.fit"` object and the optimal `lambda` value determined during
#' cross-validation.
#'
#' @param object A fitted [cv.hdhuber()] object from which coefficients are to be extracted.
#' @param s Specifies the value(s) of the penalty parameter `lambda` for which coefficients are desired.
#'   The default is `s = "lambda.1se"`, which corresponds to the largest value of `lambda` such that the
#'   cross-validation error estimate is within one standard error of the minimum. Alternatively,
#'   `s = "lambda.min"` can be used, corresponding to the minimum of the cross-validation error estimate.
#'   If `s` is numeric, these are taken as the actual values of `lambda` to use.
#' @param ... Not used.
#'
#' @return Returns the coefficients at the specified `lambda` values.
#' @seealso [cv.hdhuber()], [predict.cv.hdhuber()]
#' @method coef cv.hdhuber
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' delta <- 0.5
#' lam2 <- 0.01
#' cv.fit <- cv.hdhuber(x = x, y = y, delta = delta, lam2 = lam2)
#' coef(cv.fit, s = c(0.02, 0.03))

coef.cv.hdhuber <- function(object,
                         s = c("lambda.1se", "lambda.min"),
                         ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else stop("Invalid form for s")
  coef(object$hdhuber.fit, s = lambda, ...)
}

#' Make Predictions from a `cv.hdhuber` Object
#'
#' Generates predictions using a fitted `cv.hdhuber()` object. This function utilizes the
#' stored `hdhuber.fit` object and an optimal value of `lambda` determined during the
#' cross-validation process.
#'
#' @param newx Matrix of new predictor values for which predictions are desired.
#'   This must be a matrix and is a required argument.
#' @param object A fitted `cv.hdhuber()` object from which predictions are to be made.
#' @param s Specifies the value(s) of the penalty parameter `lambda` at which predictions
#'   are desired. The default is `s = "lambda.1se"`, representing the largest value of `lambda`
#'   such that the cross-validation error estimate is within one standard error of the minimum.
#'   Alternatively, `s = "lambda.min"` can be used, corresponding to the minimum of the
#'   cross-validation error estimate. If `s` is numeric, these are taken as the actual values
#'   of `lambda` to use for predictions.
#' @param ... Not used.
#' @return Returns a matrix or vector of predicted values corresponding to the specified
#'   `lambda` values.
#' @seealso \code{\link{cv.hdhuber}}, \code{\link{coef.cv.hdhuber}}
#' @method predict cv.hdhuber
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' delta <- 0.5
#' lam2 <- 0.01
#' cv.fit <- cv.hdhuber(x = x, y = y, delta = delta, lam2 = lam2)
#' predict(cv.fit, newx = x[50:60, ], s = "lambda.min")

predict.cv.hdhuber <- function(object, newx,
                            s = c("lambda.1se", "lambda.min"),
                            ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else stop("Invalid form for s")
  predict(object$hdhuber.fit, newx, s = lambda, ...)
}