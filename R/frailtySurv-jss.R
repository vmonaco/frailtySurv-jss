################################################################################
# Standalone script for JSS submission 2387:
# 
# General Semiparametric Shared Frailty Model: 
# Estimation and Simulation with frailtySurv
# 
# Authors:
# John V. Monaco, Malka Gorfine, Li Hsu
# 
# **Note** 
# This script is very computationally demanding: each section takes about a day.
#
# Some results are saved as RData files and then loaded during manuscript 
# compilation. This includes the case studies and simulation results.
# 
# To run as a background script and log the results:
# 
# $ nohup Rscript frailtySurv-jss.R > log 2> log.err < /dev/null
#
################################################################################

library("frailtySurv")
SEED <- 2015 # Seed used for data gen and bootstrap functions


################################################################################
# Section 2: Data generation
################################################################################

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
################################################################################

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
summary(fit, type="cumhaz", Lambda.times=c(20,50,80,110))

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
################################################################################

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
################################################################################

# Diabetic Retinopathy Study
data("drs")
fit.drs <- fitfrail(Surv(time, status) ~ treated + cluster(subject_id),
                    drs, frailty="gamma")
COV.drs <- vcov(fit.drs)

fit.drs
sqrt(diag(COV.drs))

## Warning: computationally expensive
set.seed(SEED)
COV.boot.drs <- vcov(fit.drs, boot=TRUE)

# Parameter estimation trace
plot(fit.drs, type="trace")

set.seed(SEED) # Seed the weights generated by the bootstrap procedure
plot(fit.drs, type="cumhaz", CI=0.95)

# Compare to coxph estimates
library("survival")
coxph(Surv(time, status) ~ treated + frailty.gamma(subject_id), drs)


# Hard drive failure
data("hdfail")
hdfail.sub <- subset(hdfail, grepl("WDC", model))
fit.hd <- fitfrail(Surv(time, status) ~ temp + rer + rsc
                   + psc + cluster(model),
                   hdfail.sub, frailty="gamma", fitmethod="score")
fit.hd

## Warning: computationally expensive weighted bootstrap procedure
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
plot(fit.hd, type="cumhaz", CI=0.95, end=365*6)


save(fit.drs, COV.drs, COV.boot.drs, file="drs.RData")
save(fit.hd, file="hdfail.RData")


################################################################################
# Appendix B: Simulation results
################################################################################

library("frailtySurv")
SEED <- 2015

# Default simulation parameters unless otherwise noted
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

# Includes simulations for various values of N to show the decreasing residuals
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

# See Appendix B for simulation details

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


################################################################################
# Appendix C: Performance analysis
################################################################################

library("frailtySurv")
library("survival")
library("frailtypack")
library("plyr")

SEED <- 2015

FITMETHOD <- "score"

REPS <- 100
N <- 300
K <- 2
BETA <- c(log(2), log(3))
KENDALLS.TAU <- 0.3

FRAILTY.ALL <- c("gamma", "pvf", "lognormal", "invgauss")
FRAILTY.INT <- c("lognormal", "invgauss")

lambda_0 <- function(t, c = 0.01, d = 4.6)
  (d * (c * t) ^ d) / t
Lambda_0 <-  function(t, c = 0.01, d = 4.6)
  (c * t) ^ d
Lambda_0_inv <- function(t, c = 0.01, d = 4.6)
  (t ^ (1 / d)) / c

TOLS <- c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)

# Function to compare the timings of core frailtySurv functions
compare.frailtySurv.timings <- function(reps, N) {
  params <- expand.grid(
    seed = sample(1:1e7, reps, replace = TRUE),
    frailty = FRAILTY.ALL,
    N = N,
    stringsAsFactors = FALSE
  )
  
  fn <- function(p) {
    set.seed(p[["seed"]])
    genfrail.bh.time <- system.time(dat.bh <- do.call(
      "genfrail",
      list(
        N = as.integer(p[["N"]]),
        K = K,
        beta = BETA,
        frailty = p[["frailty"]],
        theta = frailtySurv:::theta.given.tau(KENDALLS.TAU, p[["frailty"]]),
        lambda_0 = lambda_0
      )
    ))[["elapsed"]]
    set.seed(p[["seed"]])
    genfrail.cbh.time <- system.time(dat.cbh <- do.call(
      "genfrail",
      list(
        N = as.integer(p[["N"]]),
        K = K,
        beta = BETA,
        frailty = p[["frailty"]],
        theta = frailtySurv:::theta.given.tau(KENDALLS.TAU, p[["frailty"]]),
        Lambda_0 =
          Lambda_0
      )
    ))[["elapsed"]]
    set.seed(p[["seed"]])
    genfrail.icbh.time <-
      system.time(dat.icbh <- do.call(
        "genfrail",
        list(
          N = as.integer(p[["N"]]),
          K = K,
          beta = BETA,
          frailty =
            p[["frailty"]],
          theta = frailtySurv:::theta.given.tau(KENDALLS.TAU, p[["frailty"]]),
          Lambda_0_inv =
            Lambda_0_inv
        )
      ))[["elapsed"]]
    
    set.seed(p[["seed"]])
    fitfrail.loglik.time <-
      system.time(fit.loglik <- do.call(
        "fitfrail",
        list(
          formula = Surv(time, status) ~ Z1 + Z2 + cluster(family),
          dat =
            dat.icbh,
          frailty =
            p[["frailty"]],
          fitmethod =
            "loglik"
        )
      ))[["elapsed"]]
    set.seed(p[["seed"]])
    fitfrail.score.time <-
      system.time(fit.score <- do.call(
        "fitfrail",
        list(
          formula = Surv(time, status) ~ Z1 + Z2 + cluster(family),
          dat = dat.icbh,
          frailty =
            p[["frailty"]],
          fitmethod =
            "score"
        )
      ))[["elapsed"]]
    set.seed(p[["seed"]])
    vcov.time <-
      system.time(v <- do.call("vcov", list(fit.score)))[["elapsed"]]
    
    c(
      seed = p[["seed"]],
      N = as.integer(p[["N"]]),
      frailty = p[["frailty"]],
      genfrail.bh.time = genfrail.bh.time,
      genfrail.cbh.time = genfrail.cbh.time,
      genfrail.icbh.time = genfrail.icbh.time,
      fitfrail.loglik.time = fitfrail.loglik.time,
      fitfrail.score.time = fitfrail.score.time,
      vcov.time = vcov.time
    )
  }
  
  timings <-
    parallel::mclapply(
      split(params, 1:nrow(params)),
      fn,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores(),
      mc.silent = FALSE
    )
  data.frame(t(simplify2array(timings)))
}

# function to compare the timings of three different shared frailty packages
compare.package.timings <- function(reps, N) {
  params <- expand.grid(
    seed = sample(1:1e7, reps, replace = TRUE),
    frailty = c("gamma", "lognormal"),
    N = N,
    stringsAsFactors = FALSE
  )
  
  fn <- function(p) {
    if (p[["frailty"]] == "gamma") {
      coxph.frailty <- "gamma"
      frailtyPenal.frailty <- "Gamma"
    } else if (p[["frailty"]] == "lognormal") {
      coxph.frailty <- "gaussian"
      frailtyPenal.frailty <- "LogN"
    }
    
    set.seed(p[["seed"]])
    dat <- genfrail(
      N = as.integer(p[["N"]]),
      K = K,
      beta = BETA,
      frailty = p[["frailty"]],
      theta = frailtySurv:::theta.given.tau(KENDALLS.TAU, p[["frailty"]]),
      Lambda_0_inv = Lambda_0_inv
    )
    
    set.seed(p[["seed"]])
    fitfrail.time <- system.time(fit.loglik <- do.call(
      "fitfrail",
      list(
        formula = Surv(time, status) ~ Z1 + Z2 + cluster(family),
        dat = dat,
        frailty = p[["frailty"]],
        fitmethod =
          "score"
      )
    ))[["elapsed"]]
    set.seed(p[["seed"]])
    coxph.time <- system.time(fit.score <- do.call("coxph",
      list(
       formula = Surv(time, status) ~ Z1 + Z2  + frailty(family, coxph.frailty),
       data = dat
      )))[["elapsed"]]
    
    set.seed(p[["seed"]])
    frailtyPenal.time <-
      system.time(fit.score <- do.call(
        "frailtyPenal",
        list(
          formula = Surv(time, status) ~ Z1 + Z2 + cluster(family),
          data = dat,
          n.knots = 10,
          kappa = 2,
          RandDist =
            frailtyPenal.frailty
        )
      ))[["elapsed"]]
    
    c(
      seed = p[["seed"]],
      N = as.integer(p[["N"]]),
      frailty = p[["frailty"]],
      fitfrail.time = fitfrail.time,
      coxph.time = coxph.time,
      frailtyPenal.time = frailtyPenal.time
    )
  }
  
  timings <-
    parallel::mclapply(
      split(params, 1:nrow(params)),
      fn,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores(),
      mc.silent = FALSE
    )
  data.frame(t(simplify2array(timings)))
}


# function to compare the timings of three different shared frailty packages
compare.fitfrail.accuracy <-
  function(reps,
           frailty,
           fitfrail.param,
           fitfrail.param.values) {
    params <- expand.grid(
      seed = sample(1:1e7, reps, replace = TRUE),
      frailty = frailty,
      tol = fitfrail.param.values,
      stringsAsFactors = FALSE
    )
    
    fn <- function(p) {
      set.seed(p[["seed"]])
      dat <- genfrail(
        N = N,
        K = K,
        beta = BETA[1],
        frailty = p[["frailty"]],
        theta = frailtySurv:::theta.given.tau(KENDALLS.TAU, p[["frailty"]]),
        Lambda_0_inv = Lambda_0_inv
      )
      
      set.seed(p[["seed"]])
      
      fitfrail.args <-
        list(
          formula = Surv(time, status) ~ Z1 + cluster(family),
          dat = dat,
          frailty = p[["frailty"]],
          fitmethod = FITMETHOD
        )
      fitfrail.args[[fitfrail.param]] <- as.numeric(p[["tol"]])
      
      if (fitfrail.param == "abstol") {
        fitfrail.args[["reltol"]] <- 0
      } else if (fitfrail.param == "reltol") {
        fitfrail.args[["abstol"]] <- 0
      } else if (fitfrail.param == "int.abstol") {
        fitfrail.args[["int.reltol"]] <- 0
      } else if (fitfrail.param == "int.reltol") {
        fitfrail.args[["int.abstol"]] <- 0
      }
      
      fitfrail.time <-
        system.time(fit <- do.call("fitfrail", fitfrail.args))[["elapsed"]]
      
      fitfrail.beta.1.res <- (fit$beta[1] - attr(dat, "beta")[1])
      fitfrail.theta.res <- (fit$theta - attr(dat, "theta"))
      
      c(
        seed = p[["seed"]],
        frailty = p[["frailty"]],
        tol = p[["tol"]],
        fitfrail.time = fitfrail.time,
        fitfrail.beta.1.res = fitfrail.beta.1.res,
        fitfrail.theta.res = fitfrail.theta.res
      )
    }
    
    timings <-
      parallel::mclapply(
        split(params, 1:nrow(params)),
        fn,
        mc.preschedule = FALSE,
        mc.set.seed = TRUE,
        mc.cores = parallel::detectCores(),
        mc.silent = FALSE
      )
    data.frame(t(simplify2array(timings)))
  }

set.seed(SEED)
frailtySurv.timings <-
  compare.frailtySurv.timings(reps = REPS, N = seq(50, 200, 10))
save(frailtySurv.timings, file = "frailtySurv.timings.RData")

set.seed(SEED)
package.timings <-
  compare.package.timings(reps = REPS, N = seq(50, 200, 10))
save(package.timings, file = "package.timings.RData")

set.seed(SEED)
fitfrail.accuracy.abstol <- compare.fitfrail.accuracy(reps = REPS,
                                                      frailty = FRAILTY.ALL,
                                                      "abstol", TOLS)
save(fitfrail.accuracy.abstol,
     file = sprintf("fitfrail.accuracy.abstol.%s.RData", FITMETHOD))

set.seed(SEED)
fitfrail.accuracy.reltol <- compare.fitfrail.accuracy(reps = REPS,
                                                      frailty = FRAILTY.ALL,
                                                      "reltol", TOLS)
save(fitfrail.accuracy.reltol,
     file = sprintf("fitfrail.accuracy.reltol.%s.RData", FITMETHOD))

set.seed(SEED)
fitfrail.accuracy.int.abstol <- compare.fitfrail.accuracy(reps = REPS,
                                                          frailty = FRAILTY.INT,
                                                          "int.abstol", TOLS)
save(fitfrail.accuracy.int.abstol,
     file = sprintf("fitfrail.accuracy.int.abstol.%s.RData", FITMETHOD))

set.seed(SEED)
fitfrail.accuracy.int.reltol <- compare.fitfrail.accuracy(reps = REPS,
                                                          frailty = FRAILTY.INT,
                                                          "int.reltol", TOLS)
save(fitfrail.accuracy.int.reltol,
     file = sprintf("fitfrail.accuracy.int.reltol.%s.RData", FITMETHOD))

