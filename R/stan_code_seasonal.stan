data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  vector[N] cYear;  // population-level design matrix
  vector[N] cday;  // population-level design matrix
  // data needed for ARMA correlations
  int<lower=0> Kar;  // AR order
  // number of lags per observation
  int<lower=0> J_lag[N];
}
transformed data {
  int max_lag = 1;
}
parameters {
  real b;  // trend
  real s;  // seasonal term
  real Intercept;  // temporary intercept for centered predictors
  vector[Kar] ar;  // autoregressive coefficients
  real<lower=0> sigma;  // residual SD
}
transformed parameters {
}
model {
  // likelihood including constants
    // matrix storing lagged residuals
    matrix[N, max_lag] Err = rep_matrix(0, N, max_lag);
    vector[N] err;  // actual residuals
    // initialize linear predictor term
    vector[N] mu = Intercept + cYear * b + cday * s;
    // include ARMA terms
    for (n in 1:N) {
      err[n] = Y[n] - mu[n];
      for (i in 1:J_lag[n]) {
        Err[n + 1, i] = err[n + 1 - i];
      }
      mu[n] += Err[n, 1:Kar] * ar;
    }
    target += normal_lpdf(Y | mu, sigma);
  // priors including constants
  target += normal_lpdf(b | 0,1);
  target += normal_lpdf(s | 0,1);
  target += normal_lpdf(Intercept | 0,10);;
  target += cauchy_lpdf(sigma | 0,1.5)
    - 1 * cauchy_lccdf(0 | 0,1.5);
}
