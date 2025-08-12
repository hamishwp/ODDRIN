data {
  int<lower=1> N;
  real<lower=0> X[N, 6];
  int<lower=0> morts[N];
 // real phi;
}

parameters {
  real<lower=1e-3, upper = 2> event_sd;
  vector[N] event_error;
  real<lower=-40, upper=-10> loc_param;
  real<lower=0> sd_param;
  //real<lower=0> error_added;
}

transformed parameters {
  real expected_morts[N];

  for (i in 1:N) {
    expected_morts[i] =  0;
    for (j in 1:6) {
      real prob =  1/(1+exp(-(loc_param + sd_param * (j+3+event_error[i]))));
      expected_morts[i] += X[i, j] * prob;
    }
  }
}

model {
  ##event_sd ~ gamma(5,0.1);
  event_error ~ normal(0,event_sd);

  target += poisson_lpmf(morts| expected_morts);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_lpmf(morts[n] | expected_morts[n]);  // replace with your likelihood
  }
}