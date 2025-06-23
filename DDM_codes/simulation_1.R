library(RWiener)

set.seed(2025)
N <- 42          # Number of subjects
T <- 100         # Trials per subject
RTbound <- 0.1   # Lower bound for RT

simul_pars <- data.frame(
  a         = rnorm(N, 1, 0.2),
  tau       = rnorm(N, 0.3, 0.02),
  slope     = rnorm(N, 0.3, 0.5),     
  intercept = rnorm(N, -0.3, 0.5),       
  subjID    = 1:N
)

subtlex <- matrix(rlnorm(N * T, meanlog = 0.7, sdlog = 0.2), N, T) # to make subtlex is always positive value
subtlex_all <- as.vector(subtlex)
subtlex_mean <- mean(subtlex_all)
subtlex_sd   <- sd(subtlex_all)


sim_data_list <- list()


for (i in 1:N) {
  a_i       <- simul_pars$a[i]
  tau_i     <- simul_pars$tau[i]
  slope     <- simul_pars$slope[i]
  intercept <- simul_pars$intercept[i]
  subtlex_i <- subtlex[i, ]
  
  RTs     <- numeric(T)
  choices <- numeric(T)
  drifts  <- numeric(T)
  
  for (t in 1:T) {
    z <- (subtlex_i[t] - subtlex_mean) / subtlex_sd  # 표준화된 예측값
    drift <- intercept + slope * z
    drifts[t] <- drift
    
    sim <- rwiener(n = 1, alpha = a_i, tau = tau_i, beta = 0.5, delta = drift)
    RTs[t] <- sim$q
    choices[t] <- ifelse(sim$resp == "upper", 1, 2)
  }
  
  sim_data_list[[i]] <- data.frame(
    subjID     = i,
    trial      = 1:T,
    subtlex    = subtlex_i,
    subtlex_z  = (subtlex_i - subtlex_mean) / subtlex_sd,
    drift      = drifts,
    RT         = RTs,
    choice     = choices
  )
}


sim_df <- do.call(rbind, sim_data_list)
head(sim_df)
summary(sim_df$drift)


sim_df <- do.call(rbind, sim_data_list)
N <- length(unique(sim_df$subjID))
T <- max(table(sim_df$subjID))

choice  <- array(0, c(N, T))
RT      <- array(0, c(N, T))
subtlex <- array(0, c(N, T))

for (i in 1:N) {
  curSubj <- unique(sim_df$subjID)[i]
  tmp     <- subset(sim_df, subjID == curSubj)
  nTrials <- nrow(tmp)
  
  choice[i, 1:nTrials]  <- tmp$choice
  RT[i, 1:nTrials]      <- tmp$RT
  subtlex[i, 1:nTrials] <- tmp$subtlex  # 비표준화된 값 전달
}

stan_data <- list(
  N             = N,
  T             = T,
  Tsubj         = table(sim_df$subjID),
  choice        = choice,
  RT            = RT,
  subtlex       = subtlex,  # 비표준화된 subtlex 전달
  RTbound       = RTbound,
  minRT         = apply(RT, 1, min),
  subtlex_mean  = subtlex_mean,  # 전체 평균
  subtlex_sd    = subtlex_sd     # 전체 표준편차
)

simul_fit <- stan(
  file = "DDM_Lexical_1.stan",
  data = stan_data,
  pars = c("a", "tau", "intercept", "slope", "mu_pr", "sigma"),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 2,
  init = function() {
    list(tau_pr = rep(-1.5, N))
  }
)

param_sim <- rstan::extract(simul_fit)

par(mfrow = c(2, 2))

#a
a_mean = apply(param_sim$a, 2, mean)
a_sd = apply(param_sim$a, 2, sd)
plot(simul_pars$a, a_mean, abline(0,1), xlab="true a",ylab="estimated a", main="a Parameter Recovery")
arrows(x0=simul_pars$a, y0= a_mean - a_sd, y1= a_mean + a_sd, length=0.02, angle=90, code=3)

tau_mean = apply(param_sim$tau, 2, mean)
tau_sd = apply(param_sim$tau, 2, sd)
plot(simul_pars$tau, tau_mean, abline(0,1), xlab="true tau",ylab="estimated tau", main="tau Parameter Recovery")
arrows(x0=simul_pars$tau, y0= tau_mean - tau_sd, y1= tau_mean + tau_sd, length=0.02, angle=90, code=3)

intercept_mean = apply(param_sim$intercept, 2, mean)
intercept_sd = apply(param_sim$intercept, 2, sd)
plot(simul_pars$intercept, intercept_mean, abline(0,1), xlab="true intercept",ylab="estimated intercept", main="intercept Parameter Recovery")
arrows(x0=simul_pars$intercept, y0= intercept_mean - intercept_sd, y1= intercept_mean + intercept_sd, length=0.02, angle=90, code=3)

slope_mean = apply(param_sim$slope, 2, mean)
slope_sd = apply(param_sim$slope, 2, sd)
plot(simul_pars$slope, slope_mean, abline(0,1), xlab="true slope",ylab="estimated slope", main="slope Parameter Recovery")
arrows(x0=simul_pars$slope, y0= slope_mean - slope_sd, y1= slope_mean + slope_sd, length=0.02, angle=90, code=3)

