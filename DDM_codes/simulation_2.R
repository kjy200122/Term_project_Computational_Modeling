# --- 1. 패키지 설치 및 로드 ---
library(RWiener)
library(rstan)

# --- 2. 시뮬레이션 파라미터 설정 ---
set.seed(2025)

N <- 42        # Number of subjects
T <- 100       # Trials per subject
RTbound <- 0.1

# 그룹 평균 및 표준편차
simul_pars2 <- data.frame(
  a = rnorm(N, 1, 0.2),
  tau = rnorm(N, 0.30, 0.02),
  slope1 = rnorm(N, 0.5, 0.5),
  slope2 = rnorm(N, -1, 0.5),
  intercept = rnorm(N, -0.3, 0.5),
  subjID = 1:N
)

# subtlex: N x T 행렬
subtlex <- matrix(rnorm(N * T, mean = 2, sd = 0.7), N, T)
greene <- matrix(rnorm(N * T, mean = 0.7, sd = 0.7), N, T)

subtlex_all <- as.vector(subtlex)
subtlex_mean <- mean(subtlex_all)
subtlex_sd   <- sd(subtlex_all)
greene_all <- as.vector(greene)
greene_mean <- mean(greene_all)
greene_sd   <- sd(greene_all)
# --- 3. 시뮬레이션 실행 ---
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

# --- 4. 하나의 데이터프레임으로 결합 ---
sim_df <- do.call(rbind, sim_data_list)
head(sim_df)

# N: 전체 참가자 수, T: 최대 시험 횟수
N <- length(unique(sim_df$subjID))  # 참가자 수
T <- max(table(sim_df$subjID))  # 최대 시험 횟수

# choice, RT, subtlex 배열 초기화
choice <- array(0, c(N, T))
RT <- array(0, c(N, T))
subtlex <- array(0, c(N, T))
greene <- array(0, c(N, T))


# 각 참가자에 대해 데이터를 할당
for (i in 1:N) {
  curSubj <- unique(sim_df$subjID)[i]  # 참가자 ID 추출
  tmp <- subset(sim_df, subjID == curSubj)  # 해당 참가자의 데이터
  
  nTrials <- nrow(tmp)  # 해당 참가자의 시험 수
  
  # choice: 0 → 2, 1 → 1 변환 (원하는 방식으로 변환)
  choice[i, 1:nTrials] <- ifelse(tmp$choice == 1, 1, 2)
  
  # RT, subtlex는 trial 수만큼 할당
  RT[i, 1:nTrials] <- tmp$RT
  subtlex[i, 1:nTrials] <- tmp$subtlex
  greene[i, 1:nTrials] <- tmp$greene
}

# Stan 모델에 필요한 데이터 리스트 준비
stan_data <- list(
  N = N,  # 참가자 수
  T = T,  # 최대 시험 횟수
  Tsubj = table(sim_df$subjID),  # 각 참가자별 시험 횟수
  choice = choice,  # choice 배열
  RT = RT,  # RT 배열
  subtlex = subtlex,  # subtlex 배열
  greene = greene,
  RTbound = 0.1,  # RT lower bound
  minRT = apply(RT, 1, min),  # 각 참가자별 최소 RT
  subtlex_mean  = subtlex_mean,  # 전체 평균
  subtlex_sd    = subtlex_sd,     # 전체 표준편차
  greene_mean  = greene_mean,  # 전체 평균
  greene_sd    = greene_sd     # 전체 표준편차
  
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

## 결과
param_sim <- rstan::extract(simul_fit2)

# 그래프 4개를 한 화면에 배치
par(mfrow = c(2, 3))

## 각 변수에 대해 subplot 생성
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

