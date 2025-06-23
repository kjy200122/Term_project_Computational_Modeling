library(dplyr)
library(rstan)
library(ggplot2)
library(loo)
library(bayesplot)
library(posterior)

dat <- read.table("simulated_lba1.txt", header = TRUE)

allSubjs <- unique(dat$subjID)
N <- length(allSubjs)
T <- max(table(dat$subjID))

RT            <- matrix(0, nrow = N, ncol = T)
choice        <- matrix(0, nrow = N, ncol = T)
stimulus_type <- matrix(0, nrow = N, ncol = T)
subtlex_arr   <- matrix(0, nrow = N, ncol = T)
N_tr_cond     <- integer(N)

for (i in 1:N) {
  tmp <- subset(dat, subjID == allSubjs[i])
  nTrials <- nrow(tmp)
  RT[i, 1:nTrials] <- tmp$rt
  choice[i, 1:nTrials] <- tmp$choice
  stimulus_type[i, 1:nTrials] <- tmp$stimulus_type
  subtlex_arr[i, 1:nTrials] <- tmp$subtlex
  N_tr_cond[i] <- nTrials
}

minRT <- apply(RT, 1, min)

dataList <- list(
  N = N,
  Max_tr = T,
  N_choices = 2,
  N_tr_cond = N_tr_cond,
  RT = RT,
  choice = choice,
  subtlex = subtlex_arr,
  minRT = minRT,
  stimulus_type = stimulus_type
)

model_lba <- stan_model(file = "LBA_word_hier.stan")

init_fn <- function() {
  list(
    mu_int = rnorm(1, 3.19, 0.01),
    sigma_int = runif(1, 0.55, 0.60),
    
    mu_A = rnorm(1, 0.34, 0.01),
    sigma_A = runif(1, 1.18, 1.22),
    
    mu_d = rnorm(1, 0.41, 0.01),
    sigma_d = runif(1, 0.58, 0.62),
    
    mu_tau = rnorm(1, 0.01, 0.001),
    sigma_tau = runif(1, 0.28, 0.32),
    
    mu_slope = rnorm(1, 0.11, 0.01),
    sigma_slope = runif(1, 0.28, 0.32),
    
    A_raw = rnorm(N, 0, 0.05),
    d_raw = rnorm(N, 0, 0.05),
    tau_raw = rnorm(N, 0, 0.05),
    intercept_raw = rnorm(N, 0, 0.05),
    slope_raw = rnorm(N, 0, 0.05)
  )
}

# 8. Stan 실행
output <- sampling(
  object = model_lba,
  data   = dataList,
  init   = init_fn,
  iter   = 2000,
  warmup = 1000,
  chains = 4,
  cores  = 4,                  
  #control = list(adapt_delta = 0.95, max_treedepth = 14),
  refresh = 50,
  verbose = TRUE,
  pars = c("intercept", "slope", "A", "d", "tau")
)

# 1. true parameters 불러오기
true_pars <- read.table("true_parameters_lba1.txt", header = TRUE)

# 2. posterior summary 가져오기
post <- summary(output)$summary

# 3. 함수: posterior 추정치 추출
get_par_stats <- function(name, stat = c("mean", "sd")) {
  stat <- match.arg(stat)
  return(post[grep(paste0("^", name, "\\["), rownames(post)), stat])
}

# 4. 비교 대상 파라미터
params <- c("A", "d", "tau", "intercept", "slope")

# 5. 비교 플롯
par(mfrow = c(2, 3))  # 2행 3열 그래프 레이아웃

for (param in params) {
  true_val <- true_pars[[param]]
  est_mean <- get_par_stats(param, "mean")
  est_sd   <- get_par_stats(param, "sd")
  
  # plot
  plot(true_val, est_mean,
       xlab = paste("True", param),
       ylab = "Estimated",
       main = paste("Recovery:", param),
       xlim = range(true_val), ylim = range(est_mean))
  abline(0, 1, col = "red", lty = 2)
  arrows(x0 = true_val,
         y0 = est_mean - est_sd,
         y1 = est_mean + est_sd,
         angle = 90, code = 3, length = 0.05)
}

par(mfrow = c(1, 1))  # 원래 레이아웃으로 복원