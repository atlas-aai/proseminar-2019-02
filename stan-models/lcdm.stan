data {
  int<lower=1> I;                                     // number of items
  int<lower=1> J;                                     // number of respondents
  int<lower=1> K;                                     // number of attributes
  int<lower=1> C;                                     // number of profiles
  int<lower=1> N;                                     // number of observations
  int<lower=1,upper=I> ii[N];                         // item for obs n
  int<lower=1,upper=J> jj[N];                         // respondent for obs n
  int<lower=0,upper=1> y[N];                          // score for obs n
  int<lower=1,upper=N> s[J];                          // row of start for j
  int<lower=1,upper=I> l[J];                          // number of items for j
  matrix[C,K] alpha;                                  // attribute profile matrix
  matrix[I,C] xi;                                     // attribute mastery indicator
}
parameters {
  simplex[C] nu;                                      // probability of profile membership
  real intercept[I];                                  // item intercepts
  real<lower=0> maineffect[I];                        // item main effects
}
transformed parameters {
  vector[C] log_nu = log(nu);
  matrix[I,C] pi;

  // Probability of correct response on each item for each class
  for (c in 1:C) {
    for (i in 1:I) {
      pi[i,c] = inv_logit(intercept[i])^(1 - xi[i,c]) * inv_logit(intercept[i] + maineffect[i])^xi[i,c];
    }
  }
}
model{
  real ps[C];

  // Priors
  intercept ~ normal(-2, 1);
  maineffect ~ lognormal(0, 1);

  // Define model
  for (j in 1:J) {
    for (c in 1:C) {
      real log_items[l[j]];
      for (m in 1:l[j]) {
        int i = ii[s[j] + m - 1];
        log_items[m] = y[s[j] + m - 1] * log(pi[i,c]) + (1 - y[s[j] + m - 1]) * log(1 - pi[i,c]);
      }
      ps[c] = log_nu[c] + sum(log_items);
    }
    target += log_sum_exp(ps);
  }
}
generated quantities {
  matrix[J,C] prob_resp_class;
  matrix[J,K] prob_resp_attr;
  row_vector[C] prob_joint;
  real prob_attr_class[C];
  vector[N] log_lik;
  int y_rep[N];

  // Generate posterior probabilities of profile membership
  for (j in 1:J) {
    for (c in 1:C) {
      real log_items[l[j]];
      for (m in 1:l[j]) {
        int i = ii[s[j] + m - 1];
        log_items[m] = y[s[j] + m - 1] * log(pi[i,c]) + (1 - y[s[j] + m - 1]) * log(1 - pi[i,c]);
        log_lik[s[j] + m - 1] = log_nu[c] + log_items[m];
      }
      prob_joint[c] = nu[c] * exp(sum(log_items));
    }
    prob_resp_class[j] = prob_joint / sum(prob_joint);
  }

  // Generate posterior probabilities of attribute mastery
  for (j in 1:J) {
    for (k in 1:K) {
      for (c in 1:C) {
        prob_attr_class[c] = prob_resp_class[j,c] * alpha[c,k];
      }
      prob_resp_attr[j,k] = sum(prob_attr_class);
    }
  }

  // Generate posterior predictive samples
  for (j in 1:J) {
    real stu_prob = uniform_rng(0, 1);
    int jrm = 1;
    for (c in 1:C) {
      if (sum(prob_resp_class[j,1:c]) < stu_prob) {
        jrm = jrm + 1;
      }
    }
    for (m in 1:l[j]) {
      int item = ii[s[j] + m - 1];
      y_rep[s[j] + m - 1] = bernoulli_rng(pi[item, jrm]);
    }
  }
}
