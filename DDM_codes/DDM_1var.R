library(dplyr)
library(rstan)

## DATA pre-processing
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

# Processing shortest response of one subject [782733]
idx <- which(dat$idcode == 782733 & dat$trial_index == 1)
rt_median <- median(dat$rt[dat$idcode == 782733 & dat$trial_index != 1], na.rm = TRUE)

if (length(idx) == 1 && dat$rt[idx] < 0.05) {
  dat$rt[idx] <- rt_median
}

# Make [response] column
dat$response <- ifelse(dat$nat_cond == "man_made" & dat$acc == 1, "man_made", 
                       ifelse(dat$nat_cond == "man_made" & dat$acc == 0, "natural", 
                              ifelse(dat$nat_cond == "natural" & dat$acc == 1, "natural", 
                                     ifelse(dat$nat_cond == "natural" & dat$acc == 0, "man_made", NA))))


choice  <- array(0, c(N, T))
RT <- array(0, c(N, T))
subtlex <- array(0, c(N, T))
subtlex_mean <- mean(dat$subtlex_log, na.rm = TRUE)
subtlex_sd   <- sd(dat$subtlex_log, na.rm = TRUE)
minRT <- sapply(allSubjs, function(s) {
  min(dat$rt[dat$idcode == s], na.rm = TRUE)
})

for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp = subset(dat, idcode == curSubj)
  nTrials = nrow(tmp)
  choice[i, 1:nTrials] <- ifelse(tmp$response == "man_made", 2, 1)
  RT[i, 1:nTrials] <- tmp$rt
  subtlex[i, 1:nTrials] <- tmp$subtlex_log
}


dataList <- list(
  N       = N,
  T       = T,
  Tsubj   = rep(T, N),
  RT = RT,
  subtlex = subtlex,
  choice = choice,
  minRT = minRT,
  RTbound = 0.1,
  subtlex_mean = subtlex_mean,  
  subtlex_sd   = subtlex_sd     
)

#Model 1
output1 = stan("DDM_Lexical_1.stan",
              data = dataList,
              pars = c("mu_pr", "sigma", "a", "tau", "intercept", "slope","log_lik", "choice_os", "RT_os"),
              iter = 2000,
              warmup = 1000,
              chains = 4,
              cores = 2,
              init = function() {
                list(tau_pr = rep(-1.5, N))  # N은 tau_pr 길이 (예: 참가자 수)
              })


#Model1 parameters extraction
parameters_1 <- rstan::extract(output1)

# Stan Model summary
summary_1 <- summary(output1)

# Rhat extraction
summary_1$summary[grep("a", rownames(fisummary$summary)), "Rhat"]
summary_1$summary[grep("tau", rownames(fisummary$summary)), "Rhat"]
summary_1$summary[grep("slope", rownames(fisummary$summary)), "Rhat"]
summary_1$summary[grep("intercept", rownames(fisummary$summary)), "Rhat"]

# Each parameter trace plot
traceplot(output1, pars = "a")
traceplot(output1, pars = "slope")  
traceplot(output1, pars = "intercept")
traceplot(output1, pars = "tau")

# Histogram
par(mfrow = c(2, 2)) 

hist(parameters_1$mu_pr[,1])
hist(parameters_1$mu_pr[,2])
hist(parameters_1$mu_pr[,3])
hist(parameters_1$mu_pr[,4])

# histogram: sigma
hist(parameters_1$sigma[,1])
hist(parameters_1$sigma[,2])
hist(parameters_1$sigma[,3])
hist(parameters_1$sigma[,4])

plot(output1, show_density = T, pars="mu_pr")
plot(output1, show_density = T, pars="sigma")
plot(output1, show_density = T, pars="a")
plot(output1, show_density = T, pars="intercept")
plot(output1, show_density = T, pars="slope")

## LOOIC
library(loo)
log_lik_1 <- extract(output1)$log_lik  # 또는 rstan::extract(output1, "log_lik")[[1]]

loo_result <- loo(log_lik_1)
looic1 <- loo_result$estimates["looic", "Estimate"]
looic_sd <- loo_result$estimates["looic", "SE"]

print(loo_result)

