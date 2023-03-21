// Ref: https://mc-stan.org/docs/2_22/stan-users-guide/fit-gp-section.html
// Ref: https://avehtari.github.io/casestudies/Motorcycle/motorcycle_gpcourse.html
data {
  int<lower=1> N;                   // N data points
  int<lower=1> D;                   // Dimensionality
  int<lower=1> Ntest;               // test data points
  array[N] vector[D] x;             // N x D X data
  array[Ntest] vector[D] xtest;         // Ntest x D X data
  vector[N] y;                      // N x 1 Y data
}
transformed data {
  real delta = 1.0;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_f;
  real<lower=0> lengthscale_f; // lengthscale of f
  vector[N] eta;
}
model {
  vector[N] f; // mean function
  lengthscale_f ~ normal(0,1);
  sigma_f ~ normal(0,1);
  sigma ~ normal(0,1);
  eta ~ normal(0,1);

  matrix[N, N] L_K;
  matrix[N, N] K = gp_exp_quad_cov(x, sigma_f, lengthscale_f);
  for (n in 1:N)
    K[n, n] = K[n, n] + delta;

  L_K = cholesky_decompose(K);

  f = L_K * eta;
  y ~ normal(f, sigma);
}
generated quantities {
  vector[Ntest] pred;
  {

    /**************************/
    /* Let's predict new data */
    /**************************/

    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(x, sigma_f, lengthscale_f);
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);

    /******************************/
    /* Compute conditional Normal */
    /******************************/

    vector[N] K_div_y;
    K_div_y = mdivide_left_tri_low(L_K, y);
    K_div_y = mdivide_right_tri_low(K_div_y', L_K)';

    matrix[N, Ntest] k_trans;
    k_trans = gp_exp_quad_cov(x, xtest, sigma_f, lengthscale_f);
    vector[Ntest] test_mu;
    test_mu = (k_trans' * K_div_y);

    matrix[N, Ntest] v_pred;
    v_pred = mdivide_left_tri_low(L_K, k_trans);

    matrix[Ntest, Ntest] cov_test;
    cov_test = gp_exp_quad_cov(xtest, sigma_f, lengthscale_f) - v_pred' * v_pred;
  
    real jitter = 1e-8;
    pred = multi_normal_rng(test_mu, add_diag(cov_test, rep_vector(jitter, Ntest)));
  }
}
