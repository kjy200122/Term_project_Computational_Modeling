library(dplyr)
library(rstan)

#pre-processing
dat <- data %>%
  filter(target_mod != "word")
dat <- dat %>%
  filter(priming_type == "no_priming")

dat <- dat %>% 
  arrange(idcode, trial_index)

set.seed(2025)

# read the data file

allSubjs = unique(dat$idcode)  # all subject IDs
N = length(allSubjs)      # number of subjects
T <- max(table(dat$idcode))
dat$rt <- dat$rt / 1000

# idcode가 782733인 첫 번째 trial 인덱스 찾기
idx <- which(dat$idcode == 782733 & dat$trial_index == 1)

# 해당 idcode의 반응시간 중 첫 번째 rt가 너무 작으면 중간값으로 대체 (첫번째 제외)
rt_median <- median(dat$rt[dat$idcode == 782733 & dat$trial_index != 1], na.rm = TRUE)

# 조건에 맞으면 값 바꾸기
if (length(idx) == 1 && dat$rt[idx] < 0.05) {
  dat$rt[idx] <- rt_median
}

dat$response <- ifelse(dat$nat_cond == "man_made" & dat$acc == 1, "man_made", 
                       ifelse(dat$nat_cond == "man_made" & dat$acc == 0, "natural", 
                              ifelse(dat$nat_cond == "natural" & dat$acc == 1, "natural", 
                                     ifelse(dat$nat_cond == "natural" & dat$acc == 0, "man_made", NA))))



minRT <- sapply(allSubjs, function(s) {
  min(dat$rt[dat$idcode == s], na.rm = TRUE)
})


choice  <- array(0, c(N, T))
RT <- array(0, c(N, T))
subtlex <- array(0, c(N, T))
greene <- array(0, c(N, T))
subtlex_mean <- mean(dat$subtlex_log, na.rm = TRUE)
subtlex_sd   <- sd(dat$subtlex_log, na.rm = TRUE)
greene_mean <- mean(dat$greene_log, na.rm = TRUE)
greene_sd   <- sd(dat$greene_log, na.rm = TRUE)


for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp = subset(dat, idcode == curSubj)
  nTrials = nrow(tmp)
  
  choice[i, 1:nTrials] <- ifelse(tmp$response == "man_made", 2, 1)
  
  # rt, subtlex_log, greene_log도 실제 trial 수만큼 할당
  RT[i, 1:nTrials] <- tmp$rt
  subtlex[i, 1:nTrials] <- tmp$subtlex_log
  greene[i, 1:nTrials] <- tmp$greene_log
}


dataList2 <- list(
  N       = N,
  T       = T,
  Tsubj   = rep(T, N),
  RT = RT,
  subtlex = subtlex,
  greene = greene,
  choice = choice,
  minRT = minRT,
  RTbound = 0.1,
  subtlex_mean = subtlex_mean,  
  subtlex_sd   = subtlex_sd,
  greene_mean = greene_mean,  
  greene_sd   = greene_sd
)


#Model 2
output2 = stan("DDM_Lexical_2.stan", data = dataList2, 
               pars = c("mu_pr","sigma","a","tau", "intercept", "slope1","slope2", "log_lik", "choice_os", "RT_os"),
               iter = 2000, warmup=1000, chains=4, cores=2,
               init = function() {
                 list(tau_pr = rep(-1.5, N))  # N은 tau_pr 길이 (예: 참가자 수)
               })



#Model1 parameters extraction
parameters_2 <- rstan::extract(output2)

# Stan Model summary
summary_2 <- summary(output2)

# Rhat extraction
summary_2$summary[grep("a", rownames(fisummary$summary)), "Rhat"]
summary_2$summary[grep("tau", rownames(fisummary$summary)), "Rhat"]
summary_2$summary[grep("slope", rownames(fisummary$summary)), "Rhat"]
summary_2$summary[grep("intercept", rownames(fisummary$summary)), "Rhat"]

# Each parameter trace plot
traceplot(output2, pars = "a")
traceplot(output2, pars = "slope1")
traceplot(output2, pars = "slope2")  
traceplot(output2, pars = "intercept")
traceplot(output2, pars = "tau")

#histogram: mu_pr
par(mfrow = c(2, 3))  # 2행 2열로 배치

hist(parameters_2$mu_pr[,1])
hist(parameters_2$mu_pr[,2])
hist(parameters_2$mu_pr[,3])
hist(parameters_2$mu_pr[,4])
hist(parameters_2$mu_pr[,5])


#histogram: sigma
hist(parameters_1$sigma[,1])
hist(parameters_1$sigma[,2])
hist(parameters_1$sigma[,3])
hist(parameters_1$sigma[,4])
hist(parameters_1$mu_pr[,5])


plot(output2, show_density = T, pars="mu_pr")
plot(output2, show_density = T, pars="sigma")
plot(output2, show_density = T, pars="a")
plot(output2, show_density = T, pars="intercept")
plot(output2, show_density = T, pars="slope")

## LOOIC
library(loo)
log_lik_2 <- extract(output2)$log_lik  # 또는 rstan::extract(output1, "log_lik")[[1]]

loo_result2 <- loo(log_lik_2)
looic2 <- loo_result2$estimates["looic", "Estimate"]
looic2_sd <- loo_result2$estimates["looic", "SE"]

print(loo_result2)

