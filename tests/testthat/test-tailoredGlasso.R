library(tailoredGlasso)

context("test-tailoredGlasso.R")

test_that("Test tailoredGlasso", {

  # Example: Highly informative prior.
  set.seed(123)
  n <- 80
  p <- 100
  dat <- huge::huge.generator(n = n, d = p, graph = "scale-free", verbose = F)
  prec.mat <- dat$omega # true precision matrix
  prior.mat <- abs(cov2cor(prec.mat))
  res <- tailoredGlasso(dat$data, prior.mat, verbose = F)

  # Tests -----------
  # Test that the classes of results are correct.
  list.names <- c("gamma", "ebic", "k", "sparsity", "opt.sparsity", "k.opt", "opt.ebic", "lambda.common", "loglikelihoods", "thetas", "theta.opt", "w0", "opt.loglikelihood", "lambda.glasso")
  expect_equal(class(res), "list")
  expect_true(mean(names(res) %in% list.names) == 1)
  expect_equal(class(res$k.opt), "numeric")
  expect_equal(class(res$theta.opt), "matrix")
  expect_equal(class(res$thetas), "array")
  expect_equal(class(res$opt.ebic), "numeric")
  expect_equal(class(res$loglikelihoods), "numeric")
  expect_equal(class(res$opt.sparsity), "numeric")
  expect_equal(class(res$sparsity), "numeric")
  expect_equal(class(res$k), "numeric")
  expect_equal(class(res$ebic), "numeric")
  expect_equal(class(res$gamma), "numeric")
  expect_equal(class(res$lambda.common), "numeric")

  # Check that results have correct lengths.
  expect_equal(max(res$k), 99) # Check that max value of k is 99.
  expect_equal(length(res$sparsity), length(res$k)) # Check length of results
  expect_equal(length(res$loglikelihoods), length(res$ebic)) # Check length of results
  expect_equal(dim(res$thetas), c(p, p, length(res$k))) # Check dimension of array.
  expect_equal(res$sparsity[which.min(res$ebic)], res$opt.sparsity) # Check that results correspond to minimized eBIC.


  # Test that the results are valid
  expect_true(det(res$theta.opt) > 0) # Positive definite.
  expect_true(res$opt.sparsity >= 0 & res$opt.sparsity <= 1) # Valid optimal sparsity.
  expect_true(mean(res$sparsity >= 0 & res$sparsity <= 1) == 1) # Valid sparsity for all k.
  expect_equal(res$opt.ebic, min(res$ebic)) # Check that eBIC is minimized.
  expect_equal(res$opt.loglikelihood, res$loglikelihoods[which.min(res$ebic)]) # Check that loglikelihood is the correct one.
  expect_true(res$lambda.common >= 0) # Positive common lambda.
  expect_true(res$gamma >= 0) # Positive gamma.
  expect_true(res$k.opt >= 0) # Positive k.opt
  expect_true(min(res$k) == 0) # All k positive.
  expect_true(res$w0 >= 0 & res$w0 <= 1) # Valid w0.

  # Test that arguments control what they are supposed to.
  expect_true(res$opt.sparsity >= tailoredGlasso(dat$data, prior.mat, ebic.gamma = 1, verbose = F)$opt.sparsity) # Test argument ebic.gamma controls sparsity
  expect_true(res$opt.sparsity >= tailoredGlasso(dat$data, prior.mat, stars.thresh = 0.01, verbose = F)$opt.sparsity) # Test argument stars.thresh controls sparsity
  expect_output(tailoredGlasso(dat$data, prior.mat, verbose = T)) # Test that verbose=T gives output
  expect_true(res$opt.loglikelihood != tailoredGlasso(dat$data, prior.mat, verbose = T, scale = T)$opt.loglikelihood) # Test that scaling gives different results.

  # Test that providing k, w0 and lambda.opt works, and gives same results.
  expect_equal(res$opt.sparsity, tailoredGlasso(dat$data, prior.mat, lambda.opt = res$lambda.glasso, k = res$k.opt, w0 = res$w0, verbose = T)$opt.sparsity) # Verbose=T to check that this works too.


  # Test defalt arguments
  expect_equal(res$k.opt, tailoredGlasso(dat$data, prior.mat, ebic.gamma = 0, verbose = F)$k.opt) # ebic.gamma
  expect_equal(res$k.opt, tailoredGlasso(dat$data, prior.mat, stars.thresh = 0.05, verbose = F)$k.opt) # stars.thresh
  expect_equal(res$k.opt, tailoredGlasso(dat$data, prior.mat, k.max = 100, verbose = F)$k.opt) # k.max
  expect_equal(res$k.opt, tailoredGlasso(dat$data, prior.mat, scale = F, verbose = F)$k.opt) # scaling

  # Test errors
  expect_error(tailoredGlasso(dat$data, prior.mat * 10)) # Too large prior weights.
  expect_error(tailoredGlasso(dat$data, -1 * prior.mat)) # Negative prior weights.
  expect_error(tailoredGlasso(dat$data, prior.mat, k.max = 1)) # Too small k.max
  expect_error(tailoredGlasso(cov(dat$data), prior.mat)) # Covariance matrix provided, but no lambda.opt
  expect_error(tailoredGlasso(cov(dat$data), prior.mat, lambda.opt = 0.2, verbose = F)) # Covariance matrix and lambda.opt provided, but no n.
  expect_error(tailoredGlasso(dat$data, prior.mat, ebic.gamma = -1, verbose = F)) # Negative ebic.gamma
  expect_error(tailoredGlasso(dat$data, prior.mat, lambda.opt = -1)) # Negative lambda.opt
  expect_error(tailoredGlasso(dat$data, prior.mat, w0 = -1)) # Negative w0
  expect_error(tailoredGlasso(dat$data, prior.mat, w0 = 2)) # Too large w0.
  expect_error(tailoredGlasso(dat$data, prior.mat, k = -1)) # Negative k
})
