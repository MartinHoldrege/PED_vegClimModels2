#include <TMB.hpp>
// this version does elastic net (en) regularization
// of the cover weighted lognormal model
// (combination of L1 and L2)
// this version calculates likelihood based on log(y + 1)
// to remove the extreme weight yi near zero (<<1) had on the likelihood
// doing a test version where models is alpha*alpha
// while not changing inputs, that means the estimated 'alpha' is actually
// sqrt(alpha), as a result the real alpha can't be below 0.
// which mean the minimum intercept only estimate for a group is log(1 + exp(0)) = log(2)

// Numerically-stable softplus: log(1 + exp(x))
// - avoids overflow when x is large positive
// - avoids loss of precision when x is large negative
template<class Type>
Type softplus(const Type &x) {
  return logspace_add(Type(0), x);  // log(exp(0) + exp(x)) = log(1 + exp(x))
}

// Smooth absolute value for L1 to avoid non-differentiability at 0.
// (True LASSO uses |x|; this is a common, practical approximation in AD frameworks.)
template<class Type>
Type smooth_abs(const Type &x, const Type &eps) {
  return sqrt(x * x + eps);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ---- DATA (from R) ----
  DATA_VECTOR(y);          // length N, observed biomass; must be > 0 for lognormal
  DATA_MATRIX(X);          // N x P, predictors (no intercept column expected)
  DATA_MATRIX(C);          // N x G, cover weights (>= 0)
  
  DATA_SCALAR(eps_mu);     // lower bound so log(mu) is defined (lognormal needs mu>0)
  
  // Elastic net hyperparameters (inputs; not estimated)
  DATA_SCALAR(en_alpha);   // mixing: 1 = LASSO, 0 = Ridge
  DATA_SCALAR(lambda);     // overall penalty strength
  
  DATA_SCALAR(l1_eps);     // smoothing for L1
  
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
  // Due to this change actually optimizing for alpha^2
    eta.col(g).array() += alpha(g)*alpha(g);
  }
  
  
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
  // log(y_i + 1) ~ Normal(log(mu_i + 1), sigma^2)
  // vector<Type> logy  = log(y);
  // vector<Type> logmu = log(mu);
  
  vector<Type> logyp  = log(y + Type(1));
  vector<Type> logmup = log(mu + Type(1));
  
  // dnorm(..., true) returns log density; sum() gives total log-likelihood.
  // We return negative log-likelihood because optimizers minimize.
  Type nll = -sum(dnorm(logyp, logmup, sigma, true));
  
  // Always scale by N, so the regularization params
  // are comparable across training datasets of different sizes
  nll = nll / Type(N); 
  
  // ---- Elastic net penalty on B ----
  Type l1 = 0;
  Type l2 = 0;
  
  for (int p = 0; p < P; p++) {
    for (int g = 0; g < G; g++) {
      
      const Type b = B(p, g);
      
      l1 += smooth_abs(b, l1_eps);
      l2 += b * b;
    }
  }
  
  Type penalty = lambda * (en_alpha * l1 + (Type(1) - en_alpha) * l2);
  
  // ---- Objective ----
  Type obj = nll + penalty;
  
  // ---- Report pieces ----
  REPORT(nll);
  REPORT(penalty);
  REPORT(l1);
  REPORT(l2);
  REPORT(obj);
  
  return obj;
}
