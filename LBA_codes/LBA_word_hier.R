library(dplyr)
library(rstan)
library(ggplot2)
library(loo)
library(bayesplot)
library(posterior)

# Load data
load("output.RData")

set.seed(2025)

# 1. 데이터 전처리
dat <- dat %>%
  filter(target_mod != "word", priming_type == "no_priming") %>%
  arrange(idcode, trial_index)

minRT <- sapply(allSubjs, function(s) {
  min(dat$rt[dat$idcode == s], na.rm = TRUE)
})

# 3. 비정상적인 반응시간 보정
idx <- which(dat$idcode == 782733 & dat$trial_index == 1)
if (length(idx) == 1 && dat$rt[idx] < 0.05) {
  dat$rt[idx] <- median(dat$rt[dat$idcode == 782733 & dat$trial_index != 1], na.rm = TRUE)
}

# 4. 응답 재코딩
dat$response <- ifelse(dat$nat_cond == "man_made" & dat$acc == 1, "man_made", 
                       ifelse(dat$nat_cond == "man_made" & dat$acc == 0, "natural", 
                              ifelse(dat$nat_cond == "natural" & dat$acc == 1, "natural", 
                                     ifelse(dat$nat_cond == "natural" & dat$acc == 0, "man_made", NA))))

# 5. Stan 입력용 객체 생성
allSubjs <- unique(dat$idcode)
N <- length(allSubjs)
T <- max(table(dat$idcode))  # Max_tr
N_choices <- 2

RT            <- matrix(0, nrow = N, ncol = T)
stimulus_type <- matrix(0, nrow = N, ncol = T)
choice        <- matrix(0, nrow = N, ncol = T)
subtlex_arr   <- matrix(0, nrow = N, ncol = T)
N_tr_cond     <- integer(N)

for (i in seq_len(N)) {
  tmp <- filter(dat, idcode == allSubjs[i])
  nTrials <- nrow(tmp)
  
  RT[i, 1:nTrials]            <- tmp$rt
  stimulus_type[i, 1:nTrials] <- ifelse(tmp$nat_cond == "natural", 1, 2)
  choice[i, 1:nTrials]        <- ifelse(tmp$response == "man_made", 2, 1)
  subtlex_arr[i, 1:nTrials]   <- scale(tmp$subtlex_log)
  N_tr_cond[i] <- nTrials
}

# 6. Stan 데이터 리스트
dataList <- list(
  N = N,
  Max_tr = T,
  N_choices = N_choices,
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
    sigma_A = runif(1, 0.20, 0.25),
    
    mu_d = rnorm(1, 0.41, 0.01),
    sigma_d = runif(1, 0.68, 0.73),
    
    mu_tau = rnorm(1, 0.01, 0.001),
    sigma_tau = runif(1, 0.31, 0.33),
    
    mu_slope = rnorm(1, 0.11, 0.01),
    sigma_slope = runif(1, 0.04, 0.06),
    
    A_raw = rnorm(N, 0, 0.05),
    d_raw = rnorm(N, 0, 0.05),
    tau_raw = rnorm(N, 0, 0.05),
    intercept_raw = rnorm(N, 0, 0.05),
    slope_raw = rnorm(N, 0, 0.05)
  )
}

# 8. Stan 실행
output_lba <- sampling(
  object = model_lba,
  data   = dataList,
  init   = init_fn,
  iter   = 2000,
  warmup = 1000,
  chains = 4,
  cores  = 4,                  
  control = list(adapt_delta = 0.95, max_treedepth = 14),
  refresh = 50,
  verbose = TRUE,
  
  pars = c(
    # hyper‐parameters
    "intercept",  "mu_int",  "sigma_int",
    "slope",  "mu_slope","sigma_slope",
    "A",      "mu_A",    "sigma_A",
    "d",      "mu_d",    "sigma_d",
    "tau",    "mu_tau",  "sigma_tau",
    "A_raw",  "d_raw", "tau_raw", "intercept_raw", "slope_raw",
    
    # log-likelihood vector
    "log_lik"
    # pred
    ,"y_pred"
  )
)

summary_table <- summary(output_lba)$summary
write.table(summary_table, file = "posterior_summary_lba.txt", sep = "\t", quote = FALSE)

pairs(output_lba, pars = c("mu_A",    "sigma_A",
                           "mu_d",    "sigma_d",
                           "mu_tau",  "sigma_tau"),
      np = nuts_params)

pairs(output_lba, pars = c("mu_int",  "sigma_int",
                           "slope"),
      np = nuts_params)

# Extract posterior summary
print(output_lba, pars = c(
  "mu_int",    "sigma_int",
  "mu_slope",  "sigma_slope",
  "mu_A",      "sigma_A",
  "mu_d",      "sigma_d",
  "mu_tau",    "sigma_tau"
))

print(output_lba, pars = "d")
print(output_lba, pars = "tau")
print(output_lba, pars = "A_raw")
print(output_lba, pars = "log_lik")

traceplot(output_lba, size =0.05, pars=paste0("A[", 1:21, "]"))
traceplot(output_lba, size =0.05, pars=paste0("A[", 22:42, "]"))
traceplot(output_lba, size =0.05, pars=paste0("d[", 1:21, "]"))
traceplot(output_lba, size =0.05, pars=paste0("d[", 22:42, "]"))
traceplot(output_lba, size =0.05, pars=paste0("tau[", 1:21, "]"))
traceplot(output_lba, size =0.05, pars=paste0("tau[", 22:42, "]"))
traceplot(output_lba, size =0.05, pars=paste0("intercept[", 1:21, "]"))
traceplot(output_lba, size =0.05, pars=paste0("intercept[", 22:42, "]"))
traceplot(output_lba, size =0.05, pars=paste0("slope[", 1:21, "]"))
traceplot(output_lba, size =0.05, pars=paste0("slope[", 22:42, "]"))

traceplot(output_lba, pars = c("mu_int", "mu_slope", "mu_d", "mu_A", "mu_tau",
                               "sigma_int", "sigma_slope", "sigma_d", "sigma_A", "sigma_tau"),
          inc_warmup = FALSE, size =0.05)

plot(output_lba, show_density = T, pars=paste0("A[", 1:21, "]"))
plot(output_lba, show_density = T, pars=paste0("A[", 22:42, "]"))
plot(output_lba, show_density = T, pars=paste0("d[", 1:21, "]"))
plot(output_lba, show_density = T, pars=paste0("d[", 22:42, "]"))
plot(output_lba, show_density = T, pars=paste0("tau[", 1:21, "]"))
plot(output_lba, show_density = T, pars=paste0("tau[", 22:42, "]"))
plot(output_lba, show_density = T, pars=paste0("intercept[", 1:21, "]"))
plot(output_lba, show_density = T, pars=paste0("intercept[", 22:42, "]"))
plot(output_lba, show_density = T, pars=paste0("slope[", 1:21, "]"))
plot(output_lba, show_density = T, pars=paste0("slope[", 22:42, "]"))

parameters <- rstan::extract(output_lba)

par(mfrow = c(2, 3))  # 2행 2열 배치
hist(parameters$mu_A,
     main = "mu_A", xlab = "")

hist(parameters$mu_d,
     main = "mu_d", xlab = "")

hist(parameters$mu_tau,
     main = "mu_tau", xlab = "")

hist(parameters$mu_int,
     main = "mu_intercept", xlab = "")

hist(parameters$mu_slope,
     main = "mu_slope", xlab = "")

par(mfrow = c(2, 3))  # 2행 2열 배치
hist(parameters$sigma_A,
     main = "sigma_A", xlab = "")

hist(parameters$sigma_d,
     main = "sigma_d", xlab = "")

hist(parameters$sigma_tau,
     main = "sigma_tau", xlab = "")

hist(parameters$sigma_int,
     main = "sigma_intercept", xlab = "")

hist(parameters$sigma_slope,
     main = "sigma_slope", xlab = "")

# COMPUTE LOOIC  -----------------------------------------

ll <- extract_log_lik(output_lba, parameter_name = "log_lik", merge_chains = TRUE)

looic <- loo(ll)$looic

loo(ll, moment_match = TRUE)

print(looic)
