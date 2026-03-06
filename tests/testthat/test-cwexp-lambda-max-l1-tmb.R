# test the selection of lambda for selecting range of lambda's
# to test 

test_that("cwexp_lambda_max_l1_tmb returns expected structure", {

  # Shared test data and elastic-net DLL
  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  out <- cwexp_lambda_max_l1_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll_en = cwexp_dll_en
  )
  
  prep <- cwexp_prepare(dat, cwexp_test_formula, cover_cols)
  
  # Basic structure checks
  expect_true(is.list(out))
  expect_true(is.finite(out$lambda_max_l1))
  expect_gt(out$lambda_max_l1, 0)
  
  # Gradient matrix should match B dimensions
  expect_equal(dim(out$grad_B), c(prep$P, prep$G))
  
  # lambda_max_l1 should equal the max absolute gradient
  expect_equal(
    out$lambda_max_l1,
    max(abs(as.numeric(out$grad_B))),
    tolerance = 1e-10
  )
  
  # Optimizer should converge when estimating nuisance parameters
  if (!is.null(out$opt_fixB$convergence)) {
    expect_equal(out$opt_fixB$convergence, 0)
  }
})

test_that("cwexp_lambda_max_l1_tmb gives a useful shrinkage threshold", {

  # Shared test data
  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  max_obj <- cwexp_lambda_max_l1_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll_en = cwexp_dll_en
  )
  
  # Convert L1 threshold to elastic-net lambda scale
  en_alpha <- 0.5
  lambda_en <- max_obj$lambda_max_l1 / en_alpha
  
  # Fit just below and just above the threshold
  fit_below <- cwexp_fit_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll = cwexp_dll_en,
    penalty = "elastic_net",
    en_alpha = en_alpha,
    lambda = 0.8 * lambda_en,
    l1_eps = 1e-8,
    control = list(iter.max = 200, eval.max = 200)
  )
  
  fit_above <- cwexp_fit_tmb(
    data = dat,
    formula = cwexp_test_formula,
    cover_cols = cover_cols,
    dll = cwexp_dll_en,
    penalty = "elastic_net",
    en_alpha = en_alpha,
    lambda = 1.2 * lambda_en,
    l1_eps = 1e-8,
    control = list(iter.max = 200, eval.max = 200)
  )
  
  # Because the L1 penalty is smoothed, we do not expect exact zeros.
  # But coefficients should shrink more above the threshold
  # and be close to 0 above the threshold
  expect_lt(sum(abs(fit_above$par$B)), sum(abs(fit_below$par$B)))
  expect_lt(max(abs(fit_above$par$B)), max(abs(fit_below$par$B)))
  expect_true(max(abs(fit_below$par$B)/max(abs(fit_above$par$B))) > 10)
  expect_lt(max(abs(fit_above$par$B)), 0.01)
  
})

test_that("cwexp_lambda_max_l1_tmb validates l1_eps input", {

  dat <- cwexp_dummy$data
  cover_cols <- cwexp_dummy$spec$cover_cols
  
  # l1_eps must be strictly positive
  expect_error(
    cwexp_lambda_max_l1_tmb(
      data = dat,
      formula = cwexp_test_formula,
      cover_cols = cover_cols,
      dll_en = cwexp_dll_en,
      l1_eps = 0
    )
  )
})
