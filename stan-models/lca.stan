data {
  int<lower=1> I;                                     // number of items
  int<lower=1> J;                                     // number of respondents
  int<lower=1> N;                                     // number of observations
  int<lower=1,upper=I> ii[N];                         // item for obs n
  int<lower=1,upper=J> jj[N];                         // respondent for obs n
  int<lower=0,upper=1> y[N];                          // score for obs n
  int<lower=1,upper=N> s[J];                          // row of start for j
  int<lower=1,upper=I> l[J];                          // number of items for j
}
parameters {
  simplex[2] nu;
  real intercept[I];
  real<lower=0> maineffect[I];
}
transformed parameters {
  vector[2] log_nu = log(nu);
  matrix[I,2] pi;

  // Probability of correct response for each class on each item
  // class 1 = nonmaster; class 2 = master
  for (c in 1:2) {
    for (i in 1:I) {
      pi[i,c] = inv_logit(intercept[i])^(1 - (c - 1)) * inv_logit(intercept[i] + maineffect[i])^(c - 1);
    }
  }
}
model{
  real ps[2];

  // Priors
  intercept ~ normal(-2, 1);
  maineffect ~ lognormal(0, 1);

  // Define model
  for (j in 1:J) {
    for (c in 1:2) {
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
  matrix[J,2] prob_resp_class;
  row_vector[2] prob_joint;

  // Generate posterior probabilities of mastery
  for (j in 1:J) {
    for (c in 1:2) {
      real log_items[l[j]];
      for (m in 1:l[j]) {
        int i = ii[s[j] + m - 1];
        log_items[m] = y[s[j] + m - 1] * log(pi[i,c]) + (1 - y[s[j] + m - 1]) * log(1 - pi[i,c]);
      }
      prob_joint[c] = nu[c] * exp(sum(log_items));
    }
    prob_resp_class[j] = prob_joint / sum(prob_joint);
  }
}
