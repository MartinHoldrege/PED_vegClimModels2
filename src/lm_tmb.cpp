// src/lm_tmb.cpp
// example code for fitting a linear model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  
  // Parameters
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(log_sigma);
  
  // Negative log-likelihood
  Type nll = -sum(dnorm(y, a+b*x, exp(log_sigma), true));
  return nll;
}
