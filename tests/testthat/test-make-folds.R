# testing some functions in make_folds.R


source(file.path(root, 'Functions/grouping', 'make_folds.R'))

test_that("check_cluster_folds ", {
  # functions mostly as test of the make_cluster_folds function
  env_cluster <- sample(1:10, size = 100, replace = TRUE)
  folds <- make_cluster_folds(env_cluster = env_cluster, n_folds = 3, seed = 1)
  expect_true(check_cluster_folds(folds, env_cluster))
  env_cluster2 <- sample(1:11, size = 100, replace = TRUE)
  expect_error(check_cluster_folds(folds, env_cluster2))
})



