// Ref: https://mc-stan.org/docs/2_22/stan-users-guide/simulating-from-a-gaussian-process.html
data {
  int<lower=1> N;                   // N data points
  int<lower=1> D;                   // Dimensionality
  array[N] vector[D] X;             // N x D data
}
transformed data {
  matrix[N, N] K;                   // (realized) Kernel matrix
  K = gp_exp_quad_cov(X, 1., 1.0);  // RBF kernel
  for (i in 1:N)                    // regularizer
    K[i, i] = K[i, i] + 0.1;        //
  matrix[N, N] L;                   // Cholesky
  L = cholesky_decompose(K);        // for faster N(mu,K)
}
parameters {
  vector[N] y;
}
model {
  vector[N] mu = rep_vector(0, N);  // mean vector
  y ~ multi_normal_cholesky(mu, L); // sample
}
