library(RWiener)
library(rstan)

set.seed(2025)

N <- 42        # Number of subjects
T <- 100       # Trials per subject
RTbound <- 0.1

simul_pars2 <- data.frame(
  a = rnorm(N, 1, 0.2),
  tau = rnorm(N, 0.30, 0.02),
  slope1 = rnorm(N, 0.5, 0.5),
  slope2 = rnorm(N, -1, 0.5),
  intercept = rnorm(N, -0.3, 0.5),
  subjID = 1:N
)

subtlex <- matrix(rnorm(N * T, mean = 2, sd = 0.7), N, T)
greene <- matrix(rnorm(N * T, mean = 0.7, sd = 0.7), N, T)

subtlex_all <- as.vector(subtlex)
subtlex_mean <- mean(subtlex_all)
subtlex_sd   <- sd(subtlex_all)
greene_all <- as.vector(greene)
greene_mean <- mean(greene_all)
greene_sd   <- sd(greene_all)

sim_data_list <- list()
for (i in 1:N) {
  a_i     <- simul_pars2$a[i]
  tau_i   <- simul_pars2$tau[i]
  slope1  <- simul_pars2$slope1[i]
  slope2  <- simul_pars2$slope2[i]
  intercept <- simul_pars2$intercept[i]
  
  subtlex_i <- subtlex[i, ]
  greene_i  <- greene[i, ]
  
  subtlex_z <- (subtlex_i - mean(subtlex_i)) / sd(subtlex_i)
  greene_z  <- (greene_i - mean(greene_i)) / sd(greene_i)
  
  RTs <- numeric(T)
  choices <- numeric(T)
  
  for (t in 1:T) {
    drift <- intercept + slope1 * subtlex_z[t] + slope2 * greene_z[t]
    
    sim <- rwiener(n = 1, alpha = a_i, tau = tau_i, beta = 0.5, delta = drift)
    
    if (is.na(sim$q) || is.na(sim$resp)) {
      RTs[t] <- NA
      choices[t] <- NA
    } else {
      RTs[t] <- sim$q
      choices[t] <- ifelse(sim$resp == "upper", 1, 2)
    }
  }
  
  sim_data_list[[i]] <- data.frame(
    subjID = i,
    trial = 1:T,
    subtlex = subtlex_i,
    greene = greene_i,
    RT = RTs,
    choice = choices,
    drift = drift
  )
}

sim_df <- do.call(rbind, sim_data_list)
head(sim_df)

N <- length(unique(sim_df$subjID))  
T <- max(table(sim_df$subjID))  

choice <- array(0, c(N, T))
RT <- array(0, c(N, T))
subtlex <- array(0, c(N, T))
greene <- array(0, c(N, T))


for (i in 1:N) {
  curSubj <- unique(sim_df$subjID)[i]  
  tmp <- subset(sim_df, subjID == curSubj)  
  
  nTrials <- nrow(tmp)  
  
  # choice: 0 → 2, 1 → 1 
  choice[i, 1:nTrials] <- ifelse(tmp$choice == 1, 1, 2)
  
  RT[i, 1:nTrials] <- tmp$RT
  subtlex[i, 1:nTrials] <- tmp$subtlex
  greene[i, 1:nTrials] <- tmp$greene
}

stan_data <- list(
  N = N,  
  T = T,  
  Tsubj = table(sim_df$subjID),  
  choice = choice,  
  RT = RT, 
  subtlex = subtlex,  
  greene = greene,
  RTbound = 0.1, 
  minRT = apply(RT, 1, min),  
  subtlex_mean  = subtlex_mean,  
  subtlex_sd    = subtlex_sd,    
  greene_mean  = greene_mean,  
  greene_sd    = greene_sd   
  
)

simul_fit2 <- stan(
  file = "DDM_Lexical_2.stan",
  data = stan_data,
  pars = c("a", "tau", "intercept", "slope1", "slope2", "mu_pr", "sigma"),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 2,
  init = function() {
    list(tau_pr = rep(-1.5, N))  # N은 tau_pr 길이 (예: 참가자 수)
  })

param_sim <- rstan::extract(simul_fit2)

par(mfrow = c(2, 3))

#a
a_mean = apply(param_sim$a, 2, mean)
a_sd = apply(param_sim$a, 2, sd)
plot(simul_pars2$a, a_mean, abline(0,1), xlab="true a",ylab="estimated a", main="a Parameter Recovery")
arrows(x0=simul_pars2$a, y0= a_mean - a_sd, y1= a_mean + a_sd, length=0.02, angle=90, code=3)

tau_mean = apply(param_sim$tau, 2, mean)
tau_sd = apply(param_sim$tau, 2, sd)
plot(simul_pars2$tau, tau_mean, abline(0,1), xlab="true tau",ylab="estimated tau", main="tau Parameter Recovery")
arrows(x0=simul_pars2$tau, y0= tau_mean - tau_sd, y1= tau_mean + tau_sd, length=0.02, angle=90, code=3)

intercept_mean = apply(param_sim$intercept, 2, mean)
intercept_sd = apply(param_sim$intercept, 2, sd)
plot(simul_pars2$intercept, intercept_mean, abline(0,1), , xlab="true intercept",ylab="estimated intercept", main="intercept Parameter Recovery")
arrows(x0=simul_pars2$intercept, y0= intercept_mean - intercept_sd, y1= intercept_mean + intercept_sd, length=0.02, angle=90, code=3)

slope1_mean = apply(param_sim$slope1, 2, mean)
slope1_sd = apply(param_sim$slope1, 2, sd)
plot(simul_pars2$slope1, slope1_mean, abline(0,1), , xlab="true slope1",ylab="estimated slope1", main="slope1 Parameter Recovery")
arrows(x0=simul_pars2$slope1, y0= slope1_mean - slope1_sd, y1= slope1_mean + slope1_sd, length=0.02, angle=90, code=3)

slope2_mean = apply(param_sim$slope2, 2, mean)
slope2_sd = apply(param_sim$slope2, 2, sd)
plot(simul_pars2$slope2, slope2_mean, abline(0,1), xlab="true slope2",ylab="estimated slope2", main="slope2 Parameter Recovery")
arrows(x0=simul_pars2$slope2, y0= slope2_mean - slope2_sd, y1= slope2_mean + slope2_sd, length=0.02, angle=90, code=3)

