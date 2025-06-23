set.seed(2025)

# 설정값
N <- 42               # 참가자 수
num_trials <- 100     # 각 참가자당 trial 수
s <- 1.0              # 고정 drift SD

# posterior summary 불러오기
post <- read.table("posterior_summary_lba.txt", header = TRUE)

# 함수 정의
get_scalar_par <- function(name) {
  list(
    mean  = post[name, "mean"],
    sd    = post[name, "sd"]
  )
}

# 파라미터 목록
param_names <- c("int" = "int",
                 "slope" = "slope",
                 "A" = "A",
                 "d" = "d",
                 "tau" = "tau")

# 각 파라미터에 대해 mu/sigma 추출
param_values <- lapply(names(param_names), function(base) {
  mu    <- get_scalar_par(paste0("mu_", base))$mean
  sigma <- get_scalar_par(paste0("sigma_", base))$mean
  return(c(mu = mu, sigma = sigma))
})

names(param_values) <- param_names
param_values$tau["sigma"] <- 0.3
param_values$slope["sigma"] <- 0.3
param_values$A["sigma"]<-1.2

simul_pars <- data.frame(
  d         = log1p(exp(param_values$d["mu"] + param_values$d["sigma"] * rnorm(N))),
  A         = log1p(exp(param_values$A["mu"] + param_values$A["sigma"] * rnorm(N))),
  tau       = log1p(exp(param_values$tau["mu"] + param_values$tau["sigma"] * rnorm(N))),
  intercept = param_values$int["mu"] + param_values$int["sigma"] * rnorm(N),
  slope     = param_values$slope["mu"] + param_values$slope["sigma"] * rnorm(N),
  subjID    = 1:N
)

# LBA 시뮬레이션 함수
simulate_lba <- function(A, b, v, s, tau) {
  drift <- rnorm(2, mean = v, sd = s)
  drift <- pmax(drift, 1e-3)
  start <- runif(2, 0, A)
  ttf <- (b - start) / drift
  choice <- which.min(ttf)
  rt <- ttf[choice] + tau
  return(list(rt = rt, choice = choice))
}
# trial-by-trial 데이터 생성
all_data <- NULL

for (i in 1:N) {
  pars <- simul_pars[i, ]
  d <- pars$d
  A <- pars$A
  tau <- pars$tau
  intercept <- pars$intercept
  slope <- pars$slope
  b <- A + d
  
  # 자극 정보
  subtlex <- rnorm(num_trials)
  stim_type <- sample(1:2, num_trials, replace = TRUE)
  
  for (t in 1:num_trials) {
    v_t <- intercept + slope * subtlex[t]
    v_trial <- rep(1.0, 2)
    v_trial[stim_type[t]] <- v_t
    v_trial[3 - stim_type[t]] <- 1.0
    
    sim <- simulate_lba(A, b, v_trial, s, tau)
    
    all_data <- rbind(all_data, data.frame(
      subjID = i,
      trial = t,
      rt = sim$rt,
      choice = sim$choice,
      stimulus_type = stim_type[t],
      subtlex = subtlex[t]
    ))
  }
}

# 결과 저장
write.table(all_data, file = "simulated_lba1.txt", row.names = FALSE, sep = "\t")
write.table(simul_pars, file = "true_parameters_lba1.txt", row.names = FALSE, sep = "\t")

