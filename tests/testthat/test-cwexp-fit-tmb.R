# testing fitting the exp tmb lognormal models

test_that("cwexp_fit_tmb vanilla fit returns sensible output", {

  # Shared test data and vanilla DLL
  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  fit <- cwexp_fit_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll = cwexp_dll,
    control = list(iter.max = 200, eval.max = 200)
  )
  
  # Basic structure
  expect_s3_class(fit, "cwexp_tmb_fit")
  expect_true(is.numeric(fit$par$alpha))
  expect_true(is.matrix(fit$par$B))
  expect_true(is.numeric(fit$par$sigma))
  expect_gt(fit$par$sigma, 0)
  
  # Dimensions should match prepared design matrices
  prep <- cwexp_prepare(dat, cwexp_test_formula, cover_cols)
  expect_equal(length(fit$par$alpha), prep$G)
  expect_equal(dim(fit$par$B), c(prep$P, prep$G))
  
  # Predictions should be finite and positive
  mu_hat <- predict(fit, dat, type = "mu")
  expect_length(mu_hat, nrow(dat))
  expect_true(all(is.finite(mu_hat)))
  expect_true(all(mu_hat > 0))
  
  mu_true <-  predict(cwexp_dummy)
  expect_true(cor(mu_true, mu_hat) > 0.8) # recovering true values 
  
  # Optimizer should report convergence
  expect_equal(fit$tmb$opt$convergence, 0)
})

test_that("cwexp_fit_tmb elastic net returns diagnostics and shrinks B as lambda increases", {

  # Shared test data and elastic-net DLL
  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  fit_small <- cwexp_fit_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll = cwexp_dll_en,
    penalty = "elastic_net",
    en_alpha = 0.5,
    lambda = 0.001
  )
  
  fit_big <- cwexp_fit_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll = cwexp_dll_en,
    penalty = "elastic_net",
    en_alpha = 0.5,
    lambda = 0.01
  )
  
  # Both fits should succeed
  expect_s3_class(fit_small, "cwexp_tmb_fit")
  expect_s3_class(fit_big, "cwexp_tmb_fit")
  
  # REPORT() diagnostics should be available and finite
  expect_true(is.list(fit_small$report))
  expect_true(is.finite(fit_small$report$nll))
  expect_true(is.finite(fit_small$report$penalty))
  expect_true(is.finite(fit_small$report$obj))
  
  # Larger lambda should shrink the slope matrix more
  expect_lt(sum(abs(fit_big$par$B)), sum(abs(fit_small$par$B)))
})

test_that("cwexp_fit_tmb approximately recovers simulated parameters", {

  # Shared test data and vanilla DLL
  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  fit <- cwexp_fit_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll = cwexp_dll_en,
    penalty = 'elastic_net',
    # testing for no penalty (i.e. equivalent to to penalty = 'none')
    # but interested in this more used implementation
    lambda = 0
  )
  
  # Recovered alpha and B should broadly track the truth
  cor_alpha <- suppressWarnings(cor(cwexp_dummy$par$alpha, fit$par$alpha))
  cor_B <- suppressWarnings(cor(as.numeric(cwexp_dummy$par$B), as.numeric(fit$par$B)))
  
  expect_true(is.finite(cor_alpha))
  expect_true(is.finite(cor_B))
  expect_gte(cor_alpha, 0.7)
  expect_gte(cor_B, 0.7)
  
  # The fitted mean function should recover the true mean function well
  mu_true <- predict(cwexp_dummy, type = "mu")
  mu_hat <- predict(fit, dat, type = "mu")
  cor_mu <- suppressWarnings(cor(mu_true, mu_hat))
  
  expect_true(is.finite(cor_mu))
  expect_gte(cor_mu, 0.9)
})

test_that("cwexp_fit_tmb elastic net validates inputs", {

  # Shared test data and elastic-net DLL
  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  # Invalid alpha
  expect_error(
    cwexp_fit_tmb(
      data = dat,
      formula = cwexp_test_formula,
      cover_cols = cover_cols,
      dll = cwexp_dll_en,
      penalty = "elastic_net",
      en_alpha = -0.1,
      lambda = 0.01
    )
  )
  
  # Invalid lambda
  expect_error(
    cwexp_fit_tmb(
      data = dat,
      formula = cwexp_test_formula,
      cover_cols = cover_cols,
      dll = cwexp_dll_en,
      penalty = "elastic_net",
      en_alpha = 0.5,
      lambda = -1
    )
  )
  
  # Invalid smoothing constant
  expect_error(
    cwexp_fit_tmb(
      data = dat,
      formula = cwexp_test_formula,
      cover_cols = cover_cols,
      dll = cwexp_dll_en,
      penalty = "elastic_net",
      en_alpha = 0.5,
      lambda = 0.01,
      l1_eps = 0
    )
  )
})