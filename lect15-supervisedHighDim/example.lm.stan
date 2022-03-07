/* multivariate normal regression */
data {
  int n;          // number of observations
  int p;          // number of predictors
  real y[n];      // response vector
  matrix[n, p] X; // design matrix
}

parameters {
  vector[p] beta;       // the regression parameters
  real alpha;           // intercept
  real<lower=0> sigmaSq;// the standard deviation of error
}

transformed parameters {
  real<lower=0> sigma;  // the standard deviation of error
  vector[n] linpred;
  linpred = X * beta + alpha;
  sigma = sqrt(sigmaSq);
}

model {
  for (i in 1 : p) {
    beta[i] ~ normal(0, 100);
  }
  sigmaSq ~ inv_gamma(.0001, .0001);
  y ~ normal(linpred, sigma);
}
