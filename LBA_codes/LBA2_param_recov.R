library(dplyr)
library(rstan)
library(ggplot2)
library(loo)
library(bayesplot)
library(posterior)

# 데이터 불러오기
dat <- read.table("simulated_lba2.txt", header = TRUE)

allSubjs <- unique(dat$subjID)
N <- length(allSubjs)
T <- max(table(dat$subjID))

RT            <- matrix(0, nrow = N, ncol = T)
choice        <- matrix(0, nrow = N, ncol = T)
stimulus_type <- matrix(0, nrow = N, ncol = T)
subtlex_arr   <- matrix(0, nrow = N, ncol = T)
greene_arr    <- matrix(0, nrow = N, ncol = T)
N_tr_cond     <- integer(N)

for (i in 1:N) {
  tmp <- subset(dat, subjID == allSubjs[i])
  nTrials <- nrow(tmp)
  RT[i, 1:nTrials] <- tmp$rt
  choice[i, 1:nTrials] <- tmp$choice
  stimulus_type[i, 1:nTrials] <- tmp$stimulus_type
  subtlex_arr[i, 1:nTrials] <- tmp$subtlex
  greene_arr[i, 1:nTrials] <- tmp$greene
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
  greene = greene_arr,
  minRT = minRT,
  stimulus_type = stimulus_type
)

# 모델 컴파일
model_lba <- stan_model(file = "LBA_word_obj_non_centered.stan")

# 초기값 함수
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
    
    mu_slope1 = rnorm(1, 0.11, 0.01),
    sigma_slope1 = runif(1, 0.28, 0.32),
    
    mu_slope2 = rnorm(1, 0.05, 0.005),
    sigma_slope2 = runif(1, 0.28, 0.32),
    
    A_raw = rnorm(N, 0, 0.05),
    d_raw = rnorm(N, 0, 0.05),
    tau_raw = rnorm(N, 0, 0.05),
    intercept_raw = rnorm(N, 0, 0.05),
    slope1_raw = rnorm(N, 0, 0.05),
    slope2_raw = rnorm(N, 0, 0.05)
  )
}

# Stan 실행
output2 <- sampling(
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
  pars = c("intercept", "slope1", "slope2", "A", "d", "tau")
)

# 추정값 비교
true_pars <- read.table("true_parameters_lba2.txt", header = TRUE)
post <- summary(output2)$summary

get_par_stats <- function(name, stat = c("mean", "sd")) {
  stat <- match.arg(stat)
  return(post[grep(paste0("^", name, "\\["), rownames(post)), stat])
}

params <- c("intercept", "slope1", "slope2", "A", "d", "tau")

par(mfrow = c(2, 3))

for (param in params) {
  true_val <- true_pars[[param]]
  est_mean <- get_par_stats(param, "mean")
  est_sd   <- get_par_stats(param, "sd")
  
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

par(mfrow = c(1, 1))  # 레이아웃 복원
