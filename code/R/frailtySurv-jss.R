################################################################################
# Standalone script for JSS submission:
# 
# General Semiparametric Shared Frailty Model: 
# Estimation and Simulation with frailtySurv
# 
# Authors:
# John V. Monaco, Malka Gorfine, Li Hsu
# 
# **Note** 
# This script is very computationally demanding. It takes several days to
# complete using: 8 core Intel(R) Xeon(R) CPU E5-2603 0 @ 1.80GHz
# Some results are saved as RData files and loaded in knitr chunks
################################################################################
library(frailtySurv)
SEED <- 2015 # Seed used for data gen and bootstrap functions
################################################################################
# Section 2: Data generation

# Specify the baseline hazard
set.seed(SEED)
dat <- genfrail(N=300, K=2, beta=c(log(2),log(3)), 
                frailty="gamma", theta=2, 
                lambda_0=function(t, tau=4.6, C=0.01) (tau*(C*t)^tau)/t)
head(dat)

# Specify the cumulative baseline hazard
set.seed(SEED)
dat.cbh <- genfrail(N=300, K=2, beta=c(log(2),log(3)), 
                    frailty="gamma", theta=2, 
                    Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau)
head(dat.cbh, 3)

# Specify the inverse cumulative baseline hazard
set.seed(SEED)
dat.inv <- genfrail(N=300, K=2, beta=c(log(2),log(3)), 
                    frailty="gamma", theta=2, 
                    Lambda_0_inv=function(t, tau=4.6, C=0.01) (t^(1/tau))/C)
head(dat.inv, 3)

# PVF frailty
set.seed(SEED)
dat.pvf <- genfrail(N=300, K=2, beta=c(log(2),log(3)), 
                    frailty="pvf", theta=0.3, censor.rate=0.40,
                    Lambda_0_inv=function(t, tau=4.6, C=0.01) (t^(1/tau))/C)
summary(dat.pvf)

################################################################################
# Section 3: Model estimation

# Fit a model by maximizing the log-likelihood
fit <- fitfrail(Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                dat, frailty="gamma", fitmethod="loglik")
fit

# Fit a model by solving the score equations
fit.score <- fitfrail(Surv(time, status) ~ Z1 + Z2 + cluster(family), 
                      dat, frailty="gamma", fitmethod="score")
fit.score

# Summarize the survival curve
head(summary(fit), 3)
tail(summary(fit), 3)

# Summarize the baseline hazard at specific time points
summary(fit, type="hazard", Lambda.times=c(20,50,80,110))

# Estimated variance
COV.est <- vcov(fit)
sqrt(diag(COV.est))

# Weighted-bootstrap variance
set.seed(SEED)
COV.boot <- vcov(fit, boot=TRUE, B=500)
sqrt(diag(COV.boot))[1:8]

save(fit, COV.est, COV.boot, fit.score, file="fit-example.RData")

################################################################################
# Section 4: Simulation

# Note that 1000 reps are used for the simulation in the main text
set.seed(SEED)
sim <- simfrail(1000, 
                genfrail.args=alist(beta=c(log(2),log(3)), frailty="gamma", 
                                    censor.rate=0.30, N=300, K=2, theta=2,
                                    covar.distr="uniform", covar.param=c(0, 1),
                                    Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau), 
                fitfrail.args=alist(formula=Surv(time, status) ~ Z1 + Z2 
                                    + cluster(family), 
                                    frailty="gamma"),
                Lambda.times=1:120)
summary(sim)

# Using coxph
set.seed(SEED)
sim.coxph <- frailtySurv:::simcoxph(1000, 
                                    genfrail.args=alist(beta=c(log(2),log(3)), frailty="gamma",
                                                        censor.rate=0.30, N=300, K=2, theta=2,
                                                        covar.distr="uniform", covar.param=c(0, 1),
                                                        Lambda_0=function(t, tau=4.6, C=0.01) (C*t)^tau),
                                    coxph.args=alist(formula=Surv(time, status) ~ Z1 + Z2
                                                     + frailty.gamma(family)),
                                    Lambda.times=1:120)
summary(sim.coxph)

# Correlation between coefficients and frailty parameter
sapply(names(sim)[grepl("^hat.beta|^hat.theta", names(sim))],
       function(name) cor(sim[[name]], sim.coxph[[name]]))

# Mean correlation between cumulative baseline hazard
mean(sapply(names(sim)[grepl("^hat.Lambda", names(sim))],
            function(name) cor(sim[[name]], sim.coxph[[name]])), na.rm=TRUE)

save(sim, sim.coxph, file="sim1.example.RData")

################################################################################
# Section 5: Case study

# Diabetic Retinopathy Study
data(drs)
fit.drs <- fitfrail(Surv(time, status) ~ treated + cluster(subject_id), 
                    drs, frailty="gamma")
COV.drs <- vcov(fit.drs)

fit.drs
sqrt(diag(COV.drs))

## Warning: computationally expensive
set.seed(SEED)
COV.boot.drs <- vcov(fit.drs, boot=TRUE)

# Parameter estimation trace
plot(fit.drs, "trace") 

set.seed(SEED) # Seed the weights generated by the bootstrap procedure
plot(fit.drs, "hazard", CI=0.95)

# Compare to coxph estimates
library(survival)
coxph(Surv(time, status) ~ treated + frailty.gamma(subject_id), drs)


# Hard drive failure
data(hdfail)
hdfail.sub <- subset(hdfail, grepl("WDC", model))
fit.hd <- fitfrail(Surv(time, status) ~ temp + rer + rsc 
                   + psc + cluster(model), 
                   hdfail.sub, frailty="gamma", fitmethod="score")
fit.hd

## Warning: computationally expensive
set.seed(SEED)
COV <- vcov(fit.hd, boot=TRUE)

# Coefficients standard errors and pvalues
se <- sqrt(diag(COV)[c("temp","rsc","rer","psc","theta.1")])
se

pvalues <- pnorm(abs(c(fit.hd$beta, fit.hd$theta))/se, lower.tail=FALSE)*2
pvalues

# This should run quickly since the bootstraped variance was already calculated
# above and cached in the fit.hd object
set.seed(SEED)
plot(fit.hd, "hazard", CI=0.95, end=365*6)


save(fit.drs, COV.drs, COV.boot.drs, file="drs.RData")
save(fit.hd, file="hdfail.RData")

################################################################################
# Appendix B: Simulation results

# Simulation defaults
REPS <- 200 # Number of repetitions for each simulation
LAMBDA.TIMES <- 0:120 # where to evaluate the baseline hazard
CORES <- 0 # run in parallel, use all cores
N <- c(300) # Number of clusters
K <- 2 # default cluster sizes

BETA <- c(log(2),log(3))
THETA <- 2
THETA.PVF <- 0.3
FRAILTY <- "gamma" 
CENSOR.PARAM <- c(130,15)
CENSOR.RATE <- 0.30
COVAR.DISTR <- "uniform"

# Simple baseline hazard
lambda_0 = function(t, tau=4.6, C=0.01) 
  (tau*(C*t)^tau)/t

Lambda_0 = function(t, tau=4.6, C=0.01) 
  (C*t)^tau

# Increasing oscillating baseline hazard
lambda_0.osc = function(t, tau=4.6, C=0.01, A=2, f=0.1) 
  A^sin(f*pi*t) * (tau*(C*t)^tau)/t

Lambda_0.osc = Vectorize(function(t, tau=4.6, C=0.01, A=2, f=0.1) 
  ifelse(t > 0, integrate(lambda_0.osc, 0, t, subdivisions = 1000L)$value, 0))

# Function to save the simulation results
simfrail.output <- function(name, seed=SEED, ...) {
  set.seed(seed)
  assign(name, simfrail(...))
  save(list=name, file=paste(name, ".RData", sep=""))
  invisible()
}

sim1.ext <- frailtySurv:::simfrail.enum(
  reps=REPS,
  seed=SEED,
  genfrail.args=alist(
    beta=BETA,
    frailty="gamma", 
    censor.rate=CENSOR.RATE,
    K=K, theta=THETA, 
    covar.distr=COVAR.DISTR,
    Lambda_0=Lambda_0), 
  fitfrail.args=alist(
    formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
    frailty="gamma"),
  Lambda.times=LAMBDA.TIMES, 
  param.name="N", param.values=c(25,50,250,500))
save(sim1.ext, file="sim1.ext.RData")

simfrail.output("sim1",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="gamma", 
                  censor.rate=CENSOR.RATE,
                  N=N, K=K, theta=THETA, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="gamma"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim2",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="gamma", 
                  censor.rate=CENSOR.RATE,
                  N=100, K=6, theta=THETA, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="gamma"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim3",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="gamma", 
                  censor.rate=CENSOR.RATE,
                  N=N, K=K, theta=THETA, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0,
                  round.base=10), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="gamma"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim4",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="gamma", 
                  censor.rate=CENSOR.RATE,
                  N=N, K=K, theta=THETA, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0.osc), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="gamma"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim5",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="pvf", 
                  censor.rate=CENSOR.RATE,
                  N=N, K=K, theta=THETA.PVF, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="pvf"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim6",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="pvf", 
                  censor.rate=CENSOR.RATE,
                  N=N, K="poisson", K.param=c(2,0), theta=THETA.PVF, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="pvf"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim7",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="lognormal", 
                  censor.rate=CENSOR.RATE,
                  N=N, K=K, theta=THETA, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="lognormal"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)

simfrail.output("sim8",
                reps=REPS, 
                genfrail.args=alist(
                  beta=BETA,
                  frailty="invgauss", 
                  censor.rate=CENSOR.RATE,
                  N=N, K=K, theta=THETA, 
                  covar.distr=COVAR.DISTR,
                  Lambda_0=Lambda_0), 
                fitfrail.args=alist(
                  formula=Surv(time, status) ~ Z1 + Z2  + cluster(family), 
                  frailty="invgauss"),
                Lambda.times=LAMBDA.TIMES, cores = CORES)
