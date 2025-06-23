functions {
  real lba_pdf(real t, real b, real A, real v, real s) {
    real z1 = (b - A - t * v) / (t * s);
    real z2 = (b -     t * v) / (t * s);
    return (1 / A) * (-v * Phi(z1) + s * exp(normal_lpdf(z1 | 0,1))
                      + v * Phi(z2) - s * exp(normal_lpdf(z2 | 0,1)));
  }
  real lba_cdf(real t, real b, real A, real v, real s) {
    real z1 = b - A - t * v;
    real z2 = b -     t * v;
    real ts = t * s;
    return 1 + (z1 / A) * Phi(z1 / ts)
             - (z2 / A) * Phi(z2 / ts)
             + (ts / A) * (exp(normal_lpdf(z1 / ts | 0,1)) - exp(normal_lpdf(z2 / ts | 0,1)));
  }

  real lba_lpdf(real rt, int ch, vector v, real b, real A, real s, real tau) {
    real t = rt - tau;
    if (t <= 0 || v[ch] <= 0) return log(1e-10);

    real pdf  = lba_pdf(t, b, A, v[ch], s);
    real surv = 1;                        // joint survival probability

    for (k in 1:num_elements(v)) if (k != ch) {
      real cdf_k = lba_cdf(t, b, A, v[k], s);
      surv *= (1 - cdf_k);
      if (surv <= 0 || is_nan(surv)) return log(1e-10);
    }
    return log( fmax(pdf * surv , 1e-10) );
  } 

  vector lba_rng(real d, real A, vector v, real s, real tau) {

    int get_pos_drift;
    int no_pos_drift;
    int get_first_pos;
    vector[num_elements(v)] drift;
    int max_iter;
    int iter;
    real start[num_elements(v)];
    real ttf[num_elements(v)];
    int resp[num_elements(v)];
    real rt;
    vector[2] pred;
    real b;

    //try to get a positive drift rate
    get_pos_drift = 1;
    no_pos_drift  = 0;
    max_iter      = 1000;
    iter          = 0;
    while(get_pos_drift) {
      for (j in 1:num_elements(v)) {
        drift[j] = normal_rng(v[j], s);
        if (drift[j] > 0) {
          get_pos_drift = 0;
        }
      }
      iter += 1;
      if (iter > max_iter) {
        get_pos_drift = 0;
        no_pos_drift  = 1;
      }
    }
    //if both drift rates are <= 0
    //return an infinite response time
    if (no_pos_drift) {
      pred[1] = -1;
      pred[2] = -1;
    } else {
      b = A + d;
      for (i in 1:num_elements(v)) {
        //start time of each accumulator
        start[i] = uniform_rng(0, A);
        //finish times
        ttf[i] = (b-start[i])/drift[i];
      }
      //rt is the fastest accumulator finish time
      //if one is negative get the positive drift
      resp          = sort_indices_asc(ttf);
      {
        real temp_ttf[num_elements(v)];
        temp_ttf    = sort_asc(ttf);
        ttf         = temp_ttf;
      }
      get_first_pos = 1;
      iter          = 1;
      while(get_first_pos) {
        if (ttf[iter] > 0) {
          pred[1]       = ttf[iter];
          pred[2]       = resp[iter];
          get_first_pos = 0;
        }
        iter += 1;
      }
    }
    return pred;
  }
}

data {
  int<lower=1> N;
  int<lower=1> Max_tr;
  int<lower=2> N_choices;
  int<lower=1> N_tr_cond[N];
  real<lower=0> RT[N, Max_tr];
  int<lower=1,upper=N_choices> choice[N, Max_tr];
  int<lower=1,upper=2> stimulus_type[N, Max_tr];
  real subtlex[N, Max_tr];
  real minRT[N]; 
}

parameters {
  vector[N] d_raw;          real<lower=0> mu_d;     real<lower=0> sigma_d;
  vector[N] A_raw;          real<lower=0> mu_A;     real<lower=0> sigma_A;
  vector[N] tau_raw;        real<lower=0> mu_tau;   real<lower=0> sigma_tau;
  vector[N] intercept_raw;  real<lower=0> mu_int;   real<lower=0> sigma_int;
  vector[N] slope_raw;      real mu_slope;          real<lower=0> sigma_slope;
}

transformed parameters {
  vector<lower=0>[N] d   = log1p_exp(mu_d + sigma_d * d_raw);
  vector<lower=0>[N] A   = log1p_exp(mu_A + sigma_A * A_raw);
  vector<lower=0>[N] tau = log1p_exp(mu_tau + sigma_tau * tau_raw);
  vector<lower=0>[N] intercept = mu_int + sigma_int * intercept_raw;
  vector[N] slope = mu_slope + sigma_slope * slope_raw;
  vector[N] b = A + d;
}

model {
  mu_d     ~ normal(0.74, 0.05);
  sigma_d  ~ normal(0.35, 0.05) T[0,];
  mu_A     ~ normal(0.44, 0.05);
  sigma_A  ~ normal(0.15, 0.05) T[0,];
  mu_tau   ~ normal(0.20, 0.02) T[0,];
  sigma_tau~ normal(0.06, 0.01) T[0,];

  d_raw     ~ normal(0, 1);
  A_raw     ~ normal(0, 1);
  tau_raw   ~ normal(0, 1);

  mu_int   ~ normal(3.1, 0.1);
  sigma_int~ normal(0.5, 0.05) T[0,];
  intercept_raw ~ normal(0, 1);

  mu_slope ~ normal(0.1, 0.05);
  sigma_slope ~ normal(0.05, 0.02) T[0,];
  slope_raw ~ normal(0, 1);

  for (i in 1:N) {
    
    int T = N_tr_cond[i];
    vector[2] v_trial;
    target += -10 * Phi(tau[i] - (minRT[i] - 0.01));

    for (t in 1:T) {
      real v_t = intercept[i] + slope[i] * subtlex[i, t];
      if (v_t <= 0)
        target += log(1e-10);

      v_trial[stimulus_type[i, t]] = v_t;
      v_trial[3 - stimulus_type[i, t]] = 1.0;

      target += lba_lpdf(RT[i, t] | choice[i, t], v_trial,
                         b[i], A[i], 1.0, tau[i]);
    }
  }
}

generated quantities {
  real log_lik[N];
  real y_pred[N, Max_tr];

  for (i in 1:N) {
    log_lik[i] = 0;
    for (t in 1:N_tr_cond[i]) {
      vector[2] v_trial;
      real v_t = intercept[i] + slope[i] * subtlex[i, t];;
      v_trial[stimulus_type[i, t]] = v_t;
      v_trial[3 - stimulus_type[i, t]] = 1.0;

      log_lik[i] += lba_lpdf(RT[i, t] | choice[i, t], v_trial, b[i], A[i], 1.0, tau[i]);

      vector[2] pred = lba_rng(d[i], A[i], v_trial, 1.0, tau[i]);
      y_pred[i, t] = pred[1] + tau[i];
    }
  }
}
