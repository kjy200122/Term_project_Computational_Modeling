functions{
  vector wiener_rng(real a, real tau, real z, real d) {
    real dt;
    real sigma;
    real p;
    real y;
    real i;
    real aa;
    real ch;
    real rt;
    vector[2] ret;
    
    dt = .0001;
    sigma = 1;
    
    y = 0.5 * a;  // starting point
    p = .5 * (1 + ((d * sqrt(dt)) / sigma));
    i = 0;
    while (y < a && y > 0){
      aa = uniform_rng(0,1);
      if (aa <= p){
        y = y + sigma * sqrt(dt);
        i = i + 1;
      } else {
        y = y - sigma * sqrt(dt);
        i = i + 1;
      }
    }
    ch = (y <= 0) * 1 + 1;  // Upper boundary choice -> 1, lower boundary choice -> 2
    rt = i * dt + tau;
    
    ret[1] = ch;
    ret[2] = rt;
    return ret;
  }
}

data {
  int<lower=1> N;                         // Number of subjects
  int<lower=1> T;                         // Maximum number of trials
  int<lower=1> Tsubj[N];                  // Number of trials for each subject
  int<lower=-1, upper=2> choice[N, T];    // Response (-1 means NA)
  real RT[N, T];                          // Reaction times
  real subtlex[N, T];                     // Word frequency (continuous)
  real greene [N, T];
  real minRT[N];                          // Minimum RT per subject
  real RTbound;                           // RT lower bound (e.g., 0.1s)
  real subtlex_mean;
  real subtlex_sd;
  real greene_mean;
  real greene_sd;
}

parameters {
  vector[5] mu_pr;               // Group means: a, tau, d1_intercept, d2_intercept, d2_slope
  vector<lower=0>[5] sigma;      // Group std devs

  vector[N] a_pr;         // Subject-level raw parameters for boundary
  vector[N] tau_pr;       // Subject-level raw parameters for non-decision time

  vector[N] d1_pr;        // Subject-level raw for d1 intercept (subtlex intercept)
  vector[N] d2_pr;        // Subject-level raw for d1 slope (subtlex slope)
  vector[N] d3_pr;        // d2 (greene) slope

}

transformed parameters {  
  vector<lower=0>[N] a;                                       
  vector<lower=RTbound, upper=max(minRT)>[N] tau;

  vector[N] intercept;
  vector[N] slope1;
  vector[N] slope2;
  
  for (i in 1:N) {
    a[i]     = exp(mu_pr[1] + sigma[1] * a_pr[i]);
    tau[i]   = Phi_approx(mu_pr[2] + sigma[2] * tau_pr[i]) * (minRT[i] - RTbound) + RTbound; 

    intercept[i] = mu_pr[3] + sigma[3] * d1_pr[i];
    slope1[i]     = mu_pr[4] + sigma[4] * d2_pr[i];
    slope2[i]     = mu_pr[5] + sigma[5] * d3_pr[i];

  }
}

model {
  // mu_pr: [1]=log(a), [2]=tau, [3]=intercept, [4]=slope1, [5]=slope2
  mu_pr[1] ~ normal(0, 1);         // log(a)
  mu_pr[2] ~ normal(0, 1);         // tau
  mu_pr[3] ~ normal(1, 1.5);     // intercept
  mu_pr[4] ~ normal(2, 1.5);     // slope1
  mu_pr[5] ~ normal(2, 1.5);     // slope2

  
  a_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);

  d1_pr ~ normal(0, 1);
  d2_pr ~ normal(0, 1);
  d3_pr ~ normal(0, 1);


  for (i in 1:N) {
    int r;

    for (t in 1:Tsubj[i]) {
      r = choice[i, t];

    
      real drift = intercept[i] + slope1[i] *((subtlex[i, t] - subtlex_mean) / subtlex_sd) + slope2[i] * ((greene[i, t] - greene_mean) / greene_sd);
    
      if (r == 1) 
        RT[i, t] ~ wiener(a[i], tau[i], 0.5, drift);
      else if (r == 2)
        RT[i, t] ~ wiener(a[i], tau[i], 0.5, -drift);
   }
  }
}

generated quantities {
  // Group-level parameter
  real<lower=0> mu_a;
  real<lower=RTbound, upper=max(minRT)> mu_tau;
  real mu_d1;
  real mu_d2;
  real mu_d3;

  // Log-likelihood for LOOIC
  real log_lik[N];

  // Posterior predictive check
  matrix[N, T] choice_os;
  matrix[N, T] RT_os;
  vector[2] tmp_os;

  // Group-level
  mu_a   = exp(mu_pr[1]);
  mu_tau = Phi_approx(mu_pr[2]) * (mean(minRT) - RTbound) + RTbound;
  mu_d1  = mu_pr[3];
  mu_d2  = mu_pr[4];
  mu_d3  = mu_pr[5];

  for (i in 1:N) {
    for (t in 1:T) {
      choice_os[i, t] = -1;
      RT_os[i, t]     = -1;
    }
  }

  // Subject loop
  for (i in 1:N) {
    log_lik[i] = 0;

    for (t in 1:Tsubj[i]) {
      real drift = intercept[i] + slope1[i] * ((subtlex[i, t] - subtlex_mean) / subtlex_sd) + slope2[i] *((greene[i, t] - greene_mean) / greene_sd);
      int r = choice[i, t];

      if (r == 1) {
        log_lik[i] += wiener_lpdf(RT[i, t] | a[i], tau[i], 0.5, drift);
      } else if (r == 2) {
        log_lik[i] += wiener_lpdf(RT[i, t] | a[i], tau[i], 0.5, -drift);
      }

      // posterior predictive
      tmp_os = wiener_rng(a[i], tau[i], 0.5, drift);
      choice_os[i, t] = tmp_os[1];
      RT_os[i, t]     = tmp_os[2];
    }
  }
}
