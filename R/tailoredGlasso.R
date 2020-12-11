#' Perform the tailored graphical lasso
#'
#' Implements the tailored graphical lasso for data integration in Gaussian graphical models. The tailored graphical lasso utilizes useful prior information more effectively than the weighted graphical lasso without involving any risk of loss of accuracy should the prior information be misleading.
#'
#' The tailored graphical lasso is an extension of the weighted graphical lasso for graph reconstruction. The objective is to get better utilisation of the available prior information, while ensuring that the introduction of prior information may not decrease the accuracy of the resulting inferred graph. The method takes a data matrix (or a covariance matrix) for which a weighted graphical lasso graph is to be inferred, and a prior weight matrix, and transforms the prior weights to more appropriate values with the logistic function. Unless provided, the steepness parameter \eqn{k} is chosen by the extended BIC with additional edge penalty parameter \eqn{\gamma} and the sigmoid midpoint \eqn{w_0} is chosen as the lower \eqn{\beta}-quantile of the non-zero prior weights. \eqn{\beta} is the variability threshold used in the StARS selection of the penalty parameter \eqn{\lambda} for the unweighted graphical lasso graph of the data. If \code{x} is a covariance matrix, StARS selection of \eqn{\lambda} cannot be performed and the parameter must be provided by the user.
#'
#' @param x Either an \eqn{n} by \eqn{p} data matrix or a \eqn{p} by \eqn{p} empirical covariance matrix. If \code{x} is a covariance matrix, tailoredGlasso cannot perform StARS selection of \eqn{\lambda} for the unweighted grapihcal lasso graphs and \code{lambda.opt} must be provided. This is only recommended if another selection routine for \eqn{\lambda} than StARS is desired. The function automatically identifies the input matrix by checking the symmetry. \eqn{n} is the sample size and \eqn{p} is the dimension.
#'
#' @param prior.matrix The \eqn{p} by \eqn{p} prior weight matrix to use. Must be symmetric and have elements in \eqn{[0,1]}.
#'
#' @param n The number of observations used in the computation of the empirical covariance matrix. This quantity is used to compute the value of log-likelihood. Only required if \code{x} is an empirical covariance matrix.
#'
#' @param ebic.gamma The value of \eqn{\gamma} to use in the extended BIC (eBIC) selection criterion. Negative values are not valid. The default value is \eqn{0}.
#'
#' @param k.max The maximum value to consider for the steepness parameter \eqn{k} of the logistic function.
#'
#' @param w0 If provided, the sigmoid midpoint to use for the logistic function. If \code{NULL}, it will be set to be the lower \code{stars.thresh}-quantile of the non-zero prior weights. Must be in \eqn{[0,1]}.
#'
#' @param k If provided, no selection of \eqn{k} will be performed. Instead, the method will perform the weighted graphical lasso value with the prior weights transformed by the logistic function with logistic function with steepness parameter \eqn{k}. \code{NULL} by default.
#'
#' @param lambda.opt If provided, the optimal \eqn{\lambda} for the unweighted graph. If not provided, it will be selected by StARS.
#'
#' @param scale If \code{scale=TRUE}, all variables will be scaled. Default value is \code{FALSE}.
#'
#' @param stars.thresh The variability threshold to use in the StARS selection of \eqn{\lambda} for the unweighted graph if \code{lambda.opt=NULL}.The default value is \eqn{0.05}. Will also be used to determine the sigmoid midpoint if \code{w0=NULL}.
#'
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#'
#'
#' @return Object of class \code{"list"}. Contains the following items:
#' \describe{
#'   \item{theta.opt}{The tailored graphical lasso precision matrix estimate.}
#'   \item{k.opt}{The selected value of \eqn{k} that minimizes the eBIC.}
#'   \item{opt.ebic}{The eBIC value of the optimal model.}
#'   \item{opt.sparsity}{The sparsity of the optimal model.}
#'   \item{lambda.common}{The value of \eqn{\lambda} used in the optimal tailored graphical lasso graph.}
#'   \item{k_vals}{All values of the steepness parameter \eqn{k} in the logistic function considered in the selection.}
#'   \item{gamma}{The value of \eqn{\gamma} that was used in the eBIC selection.}
#'   \item{eBIC}{The eBIC scores of the models corresponding to the different values \eqn{k}.}
#'   \item{sparsity}{The sparsities of the models corresponding to different \eqn{k}.}
#'   \item{loglikelihoods}{The log likelihood values of the models corresponding to different \eqn{k}.}
#'   \item{thetas}{The estimated precision matrices corresponding to different \eqn{k}.}
#'   \item{w0}{The sigmoid midpoint used.}
#'   \item{opt.loglikelihood}{The log likelihood value of the selected precision matrix.}
#'   \item{lambda.glasso}{The value of \eqn{\lambda} selected for the unweighted graphical lasso graph.}

#' }
#'
#'
#' @seealso \code{\link[huge]{huge}}, \code{\link[huge]{huge.select}}, \code{\link[glasso]{glasso}}
#'
#' @export
#'
#' @author Camilla Lingjaerde
#'
#'
#'
#' @examples
#'
#' # example 1: simple example where the prior weight matrix is
#' # unrelated to the data of interest.
#' # generate data and prior matrix
#' set.seed(123)
#' x <- matrix(rnorm(40 * 20), ncol = 20)
#' x2 <- matrix(rnorm(40 * 20), ncol = 20)
#' prior.mat <- abs(cov2cor(var(x2)))
#' res <- tailoredGlasso(x, prior.mat)
#' res$theta.opt # the estimated tailoredGlasso precision matrix for x.
#' res$k.opt # the optimal selected value of k.
#' res$opt.sparsity # the sparsity of the estimated precision matrix.
#'
#' # example 2: scaling the data
#' set.seed(123)
#' res <- tailoredGlasso(x, prior.mat, scale = TRUE)
#'
#' # example 3: scale-free data where prior weight matrix is
#' # highly informative for the data of interest.
#' set.seed(123)
#' n <- 80
#' p <- 100
#' dat <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' prec.mat <- dat$omega # true precision matrix
#' prior.mat <- abs(cov2cor(prec.mat))
#' res <- tailoredGlasso(dat$data, prior.mat, scale = TRUE)
#' res$k.opt
#' precision(abs(prec.mat) > 1e-7, res$theta.opt != 0)
#' # k is chosen very large, high precision.
#'
#' # example 4: scale-free data where prior weight matrix is
#' #  completely uninformative for the data of interest.
#' # Create a completely unrealted prior data set
#' set.seed(123)
#' n <- 80
#' p <- 100
#' dat <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' dat.prior <- huge::huge.generator(n = n, d = p, graph = "scale-free")
#' prec.mat.prior <- dat.prior$omega # true precision matrix
#' prior.mat <- abs(cov2cor(prec.mat.prior))
#' res <- tailoredGlasso(dat$data, prior.mat, scale = TRUE)
#' res$k.opt
#' precision(abs(dat$omega) > 1e-7, res$theta.opt != 0)
#' # k is chosen to be very small, lower precision as less useful prior information.
tailoredGlasso <- function(x, prior.matrix, ebic.gamma = 0, k.max = 100, stars.thresh = 0.05, n = NULL, w0 = NULL, k = NULL, lambda.opt = NULL, scale = FALSE, verbose = TRUE) {

  # Check that the prior matrix has elements in [0,1].
  if (max(prior.matrix) > 1 | min(prior.matrix) < 0) {
    stop("the elements in prior.matrix must be between 0 and 1. \n")
  }
  # Check that w0 has a valid value if provided.
  if (!is.null(w0)) {
    if (w0 < 0 | w0 > 1) {
      stop("w0 must be between 0 and 1. \n")
    }
    else {
      if (verbose) cat("w0 has been provided, so it will not be selected by tailoredGlasso. \n")
    }
  }
  # Check that k has a valid value if provided.
  if (!is.null(k)) {
    if (k < 0) {
      stop("k cannot be negative. \n")
    }
  }
  # Check that lambda.opt has a valid value if provided.
  if (!is.null(lambda.opt)) {
    if (lambda.opt < 0) {
      stop("lambda.opt cannot be negative. \n")
    }
    else {
      if (verbose) cat("lambda.opt has been provided, so it will not be selected by tailoredGlasso. \n")
    }
  }
  # Check that ebic.gamma has a valid value
  if (ebic.gamma < 0) {
    stop("ebic.gamma cannot have a negative value. \n")
  }
  # Check that k.max has a valid value
  if (k.max < 5) {
    stop("larger k.max should be considered for appropriate selection. Try k.max = 50. \n")
  }
  # Check if the provided matrix is a data matrix or covariance matrix.
  if (!isSymmetric(x)) {
    n <- nrow(x)
    cov.mat <- stats::cov(x)
  }
  # if x is a covariance matrix, lambda.opt and n must be provided.
  else {
    cov.mat <- x
    if (is.null(lambda.opt)) {
      stop("penalty parameter selection for the unweighted graphical lasso graph cannot be performed when x is a covariance matrix. \n")
    }
    if (is.null(n)) {
      stop("n must be provided when x is a covariance matrix. \n")
    }
  }
  p <- nrow(cov.mat)

  # Scale the data if asked to.
  if (scale) {
    if (verbose) cat("Scaling data...\n")
    cov.mat <- stats::cov2cor(cov.mat)
    x <- scale(x)
  }

  # Unless sigmoid midpoint is explicitly provided, select it as the lower stars.thresh-quantile of the non-zero prior weights.
  if (is.null(w0)) {
    w0 <- stats::quantile(prior.matrix[which(prior.matrix != 0 & prior.matrix != 1, arr.ind = T)], stars.thresh)
  }

  # If the optimal penalty parameter lambda.opt for the unweighted graph is not provided, select it by StARS.
  if (is.null(lambda.opt)) {
    if (verbose) cat("Selecting lambda for the unweighted graph...\n")
    fit.huge <- huge::huge(x, method = "glasso", nlambda = 35, verbose = F)
    fit.stars <- huge::huge.select(fit.huge, criterion = "stars", stars.thresh, verbose = F)
    lambda.opt <- fit.stars$opt.lambda
  }

  # If k is provided, no selection routine is required.
  if (!is.null(k)) {
    if (verbose) cat("Since k is provided, its value will not be selected by tailoredGlasso. \n")
    ans <- list()
    fit.w <- tailoredGlasso_given_k(lambda.opt, cov.mat, prior.matrix, k, p, w0)
    theta.est <- fit.w[[1]]$wi
    # Make sure to avoid rounding errors.
    theta.est[which(abs(theta.est) < 1e-7)] <- 0
    ans$gamma <- ebic.gamma
    ans$ebic <- ans$opt.ebic <- tailoredGlasso::eBIC(theta.est, sample.cov = cov.mat, n = n, gamma = ebic.gamma)
    ans$k <- ans$k.opt <- k
    ans$sparsity <- ans$opt.sparsity <- tailoredGlasso::sparsity(theta.est != 0)
    ans$lambda.common <- fit.w$lambda
    ans$loglikelihoods <- ans$opt.loglikelihood <- tailoredGlasso::gaussianloglik(cov.mat, theta.est, n)
    ans$thetas <- ans$theta.opt <- theta.est
    ans$w0 <- w0
    ans$lambda.glasso <- lambda.opt
    return(ans)
  }

  # If k is not provided, consider a grid of values. The grid is finer in [0,5), as the results are sensitive to values in this area.
  k_vals <- c(seq(0, 4.99, by = 0.01), seq(5, k.max, by = 2))
  likelihoods <- rep(0, length(k_vals))
  theta.hats <- array(dim = c(p, p, length(k_vals)))
  lambdas <- rep(0, length(k_vals))

  # Fit models, one for each value of k:
  if (verbose) cat("Selecting k...\n")
  for (i in 1:length(k_vals)) {
    # Logistic transformation of the prior weights.
    weights <- matrix(1 - stats::plogis(prior.matrix, location = w0, scale = 1 / k_vals[i]), nrow = p, byrow = T)
    # New common penalty parameter found by penalty preservation rule.
    lambdas[i] <- (lambda.opt * p^2) / (sum(weights))
    # Fit weighted graphical lasso model with new penalty parameter and transformed weights.
    fit.w <- glasso::glasso(cov.mat, rho = lambdas[i] * weights, penalize.diagonal = F)
    # Save estimated precision matrix.
    theta.hats[, , i] <- fit.w$wi
    # Save the Multivariate Gaussian log likelihood.
    likelihoods[i] <- tailoredGlasso::gaussianloglik(cov.mat, fit.w$wi, n)
    # Print how far we have come last time the percentage done is dividable by 10.
    if (verbose) {
      done <- round(100 * i / length(k_vals))
      done.next <- round(100 * (i + 1) / length(k_vals))
      if (i == length(k_vals) | (done %% 10) == 0 & (done.next %% 10) != 0) cat(done, " % done \n")
    }
  }
  # The optimal model is the one minimizing the eBIC score.
  ans <- list()
  ebic <- apply(theta.hats, 3, tailoredGlasso::eBIC, sample.cov = cov.mat, n = n, gamma = ebic.gamma)
  ans$gamma <- ebic.gamma
  ans$ebic <- ebic
  ans$k <- k_vals
  ans$sparsity <- apply(theta.hats, 3, sparsity)
  ans$opt.sparsity <- tailoredGlasso::sparsity(theta.hats[, , which.min(ebic)])
  ans$k.opt <- k_vals[which.min(ebic)]
  ans$opt.ebic <- min(ebic)
  ans$lambda.common <- lambdas[which.min(ebic)]
  ans$loglikelihoods <- likelihoods
  ans$thetas <- theta.hats
  ans$theta.opt <- theta.hats[, , which.min(ebic)]
  ans$w0 <- w0
  ans$opt.loglikelihood <- likelihoods[which.min(ebic)]
  ans$lambda.glasso <- lambda.opt
  return(ans)
}

#' @keywords internal
tailoredGlasso_given_k <- function(lambda.opt, cov.mat, prior.matrix, k, p, w0 = NULL) {
  # Logistic transformation of prior weights.
  weights <- as.matrix(1 - stats::plogis(prior.matrix, location = w0, scale = 1 / k), nrow = p, byrow = T)
  # New common penalty parameter found by penalty preservation rule.
  lambda <- (lambda.opt * p^2) / (sum(weights))
  # Fit weighted graphical lasso model with new penalty parameter and transformed weights.
  fit.w <- glasso::glasso(cov.mat, rho = lambda * weights, penalize.diagonal = F)
  return(list(fit.w, lambda))
}
