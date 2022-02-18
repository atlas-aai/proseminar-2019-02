data {
  int<lower=0> I;                       // number of items
  int<lower=0> J;                       // number of students
  int<lower=0> D;                       // number of dimensions
  int<lower=0> N;                       // number of observations
  int<lower=0,upper=I> ii[N];              // item for observation n
  int<lower=0,upper=J> jj[N];              // student for observation n
  int<lower=0,upper=D> dd[N];              // dimension for observation n
  int<lower=0,upper=1> y[N];               // score for observation n
}
transformed data {
  vector[D] theta_mean;
  vector[D] theta_scale;
  
  for (d in 1:D) {
    theta_mean[d] = 0;
    theta_scale[d] = 1;
  }
}
parameters {
  vector[D] theta[J];                   // student ability
  corr_matrix[D] theta_corr;            // theta correlation matrix
  
  vector<lower=0>[I] alpha;             // discriminations
  vector[I] beta;                       // difficulty
}
model {
  vector[N] eta;
  vector[D] theta_j;
  
  // Priors
  theta_corr ~ lkj_corr(1);
  theta ~ multi_normal(theta_mean, quad_form_diag(theta_corr, theta_scale));
  alpha ~ lognormal(0.5, 1);
  beta ~ normal(0, 10);
  
  // Likelihood
  for (n in 1:N) {
    theta_j = theta[jj[n]];
    eta[n] = alpha[ii[n]] * (theta_j[dd[n]] - beta[ii[n]]);
  }
  y ~ bernoulli_logit(eta);
}
