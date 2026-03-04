# test_lm_tmb_fit.R
# for testing fitting a linear model

library(TMB)

if(FALSE) {
  # 1) confirm TMB is installed where you think
  find.package("TMB")
  
  # 2) confirm RcppEigen is installed
  install.packages("RcppEigen")  # safe to run even if it is installed
  
  # 3) confirm the header actually exists
  file.exists(file.path(find.package("RcppEigen"), "include", "Eigen", "Dense"))
}

# --- compile and load model ---
cpp <- file.path("src", "lm_tmb.cpp")
TMB::compile(cpp)

dll <- sub("\\.cpp$", "", cpp) # "src/lm_tmb"
dyn.load(TMB::dynlib(dll))


# --- simulate data ---
set.seed(1)
n <- 1e5
x <- rnorm(n)
a_true <- 2
b_true <- -0.7
sigma_true <- 1.3
y <- a_true + b_true * x + rnorm(n, sd = sigma_true)

# fit model ---

data <- list(y = y, x = x)
parameters <- list(a = 0, b = 0, log_sigma = 0)


lm_tmb <- function(data, parameters) {
  # --- build AD function object ---
  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    DLL = basename(dll)   # IMPORTANT on Windows: just "lm_tmb"
  )
  
  # --- optimize (uses obj$fn and obj$gr) ---
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
  opt

}

# comparing with vanilla optimization approach
nll_R <- function(par, y, x) {
  a <- par[1]; b <- par[2]; sigma <- exp(par[3])
  mu <- a + b*x
  -sum(dnorm(y, mu, sigma, log = TRUE))
}

lm_R_mle <- function(data) {
  optim(c(0,0,0), nll_R, y=data$y, x=data$x, method="BFGS")
}

# compare speeds
results <- microbenchmark::microbenchmark(
  lm_tmb(data, parameters),
  lm_R_mle(data),
  times = 1
)
results
