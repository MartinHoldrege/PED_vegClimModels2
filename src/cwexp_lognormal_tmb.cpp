#include <TMB.hpp>
// this version does not do regularization (just
// for optimizing log likelihood)

// Numerically-stable softplus: log(1 + exp(x))
// - avoids overflow when x is large positive
// - avoids loss of precision when x is large negative
template<class Type>
Type softplus(const Type &x) {
  return logspace_add(Type(0), x);  // log(exp(0) + exp(x)) = log(1 + exp(x))
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ---- DATA (from R) ----
  DATA_VECTOR(y);          // length N, observed biomass; must be > 0 for lognormal
  DATA_MATRIX(X);          // N x P, predictors (no intercept column expected)
  DATA_MATRIX(C);          // N x G, cover weights (>= 0)
  
  // DATA_SCALAR(cap_eta);    // cap linear predictor to avoid extreme values (optional safety)
  DATA_SCALAR(eps_mu);     // lower bound so log(mu) is defined (lognormal needs mu>0)
  
  // ---- PARAMETERS (estimated) ----
  PARAMETER_VECTOR(alpha); // length G, group-specific intercepts
  PARAMETER_MATRIX(B);     // P x G, slopes for each predictor x group
  PARAMETER(log_sigma);    // scalar; we estimate log(sigma) so sigma stays > 0
  
  const Type sigma = exp(log_sigma);
  
  // ---- Basic dimension checks (fail fast) ----
  const int N = y.size();
  const int P = X.cols();
  const int G = C.cols();
  
  if (X.rows() != N)        error("X rows must match length(y).");
  if (C.rows() != N)        error("C rows must match length(y).");
  if (alpha.size() != G)    error("alpha length must match ncol(C).");
  if (B.rows() != P)        error("B rows must match ncol(X).");
  if (B.cols() != G)        error("B cols must match ncol(C).");
  
  // ---- Linear predictor: eta = X * B + alpha ----
  // eta is N x G, with:
  //   eta_ig = alpha_g + sum_p X_ip * B_pg
  matrix<Type> eta = X * B; // N x G
  
  // Add group intercept alpha_g to each column g
  for (int g = 0; g < G; g++) {
    eta.col(g).array() += alpha(g);
  }
  
  // Optional cap to prevent absurd values (softplus is stable, but predictions can still blow up)
  // NOTE: this caps only above; if you also want to cap below, use cwiseMax(min_eta).
  // eta = eta.cwiseMin(cap_eta);
  
  // ---- Mean function: mu_i = sum_g C_ig * softplus(eta_ig) ----
  // We still do a loop over rows because we need a row-wise reduction (sum over g).
  vector<Type> mu(N);
  for (int i = 0; i < N; i++) {
    Type mui = 0;
    for (int g = 0; g < G; g++) {
      // Each group contributes a positive amount via softplus(eta_ig),
      // weighted by cover C_ig.
      mui += C(i, g) * softplus(eta(i, g));
    }
    // Keep mu positive so log(mu) is defined (lognormal model).
    if (mui < eps_mu) mui = eps_mu;
    mu(i) = mui;
  }
  
  // ---- Lognormal likelihood on log scale ----
  // log(y_i) ~ Normal(log(mu_i), sigma^2)
  vector<Type> logy  = log(y);
  vector<Type> logmu = log(mu);
  
  // dnorm(..., true) returns log density; sum() gives total log-likelihood.
  // We return negative log-likelihood because optimizers minimize.
  Type nll = -sum(dnorm(logy, logmu, sigma, true));
  
  return nll;
}
