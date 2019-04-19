#' @title Simulation of a series of three-arm trials with an interim analysis with a stopping rule for futility based on O'Brien-Fleming spending function.
#' @author Mohamed Amine BAYAR
#' \code{SimulateThreeArmOBFFutility} simulates a series of three-arm trials with an interim analysis with a stopping rule for futility based on O'Brien-Fleming spending function.
#' 
#' @param accrual  An integer specifying the accrual rate per year.
#' @param alpha    A real specifying the Type I error.
#' @param beta     A real specifying the Type II error.
#' @param hr       A real specifying the Hazard ratio Exp/ctl under the alternative hypothesis (target hazard ration \delta).
#' @param distrib  An integer specifying the index of the distribution of futur treatments effects, See DistributionFutureTreatmentEffects.R for more details
#' @param baseline A real specifying the A positive real number specifying the hazard rate associated with the control arm of the first trial of the series..
#' @param FU       A real specifying the follow-up period in years.
#' @param nbrepeat An integer specifying the number of repetitions.
#' @param Tmax     A real specifying the series duration ('on average').
#' 
#' @return a table of the results, trial by trial, series by series, and repetition by repetition.

rm(list = ls())
setwd(dir = 'G:/MohamedAmineBayar/SmallTrialsStrategy/SeriesIAMAMS/GitHub')
# ----- SOURCES - PACKAGES - SEED -----
library(gsDesign)
library(survival)
library(asd)
library(sqldf)
source("DistributionFutureTreatmentEffects.R")
source("MultipleTestingCombination.R")

# ----- DECLARE CONSTANCES -----
ns         <- 2  # Number of stages
nexp       <- 2  # Import number of experimental treatment tested
timing     <- 1  # Stages are equally spaced (of equal length) 
test.type  <- 4  # 4=two-sided, asymmetric, beta-spending with non-binding lower bound (not used), alpha-spending upper bound (not used)
hr0        <- 1  # Hazard ratio experimental/control under the null hypothesis
ratio      <- 1  # Randomization ratio experimental/control
sided      <- 1  # 1-sided testing
ia.type    <- "OBFFutility"

SimulateThreeArmOBFFutility  <- function(accrual, alpha, beta, hr, distrib, baseline, FU, nbrepeat, Tmax){  # ----- START OF SIMUL FUNCTION -----
  
  # ----- RESULTS STORAGE -----
  ## Simulation parameters ID  
  accrual.rate.trial   <- numeric()  # Vector of the accrual rate per year within the series
  lambda.trial         <- numeric()  # Vector of the baseline exponential survival distribution hazard rate in the control group within the series
  FU.trial             <- numeric()  # Vector of the fu within the series
  distribution.trial   <- numeric()  # Vector of the hypothetical distributions of futur treatments effect within the series
  beta.trial           <- numeric()  # Vector of the type II error within the series
  hr.alternative.trial <- numeric()  # Vector of the hazard ratio of alternative hypothsesis within the series
  ## Design parameters ID
  alpha.trial     <- numeric()  # Vector of the type I error within the series
  weight.trial    <- numeric()  # Vector of the weight of the first stage, fixed a priori to 0.5 vs. 0.5
  nexp.trial      <- numeric()  # Vector of the number of experimental treatment tested
  ia.type.trial   <- numeric()  # Vector of the type of interim analysis within the series
  ## Series (repetition) ID 
  series.trial <- numeric()  # Vector of the series ID
  ## Trial ID
  ntrial.trial     <-  numeric()  # Vector of the trial ID
  time.start.trial <-  numeric()  # Vector of the time the trial starts
  ## Previuous trial legacy
  lambdaC.inherited.trial     <- numeric()  # Vector of the control group hazard rate (inherited from the previous trial) of each trial within the series state of nature
  lambdaC.inherited.obs.trial <- numeric()  # Vector of the control group hazard rate (inherited from the previous trial) of each trial within the series estimated accorindg to simulation (in the previous trial) 
  ## Design characteristics
  n.trial    <- numeric()  # Vector of the required sample size of each trial within the series
  Texp.trial <- numeric()  # Vector of the expected trial duration of each trial within the series   
  Ts1.trial  <- numeric()  # Vector of the expected first stage ending time (50% of expected events) of each trial within the series  
  Tacc.trial <- numeric()  # Vector of the expected accrual duration of each trial within the series
  ## Interim analysis bounds
  pvalue.futility.ia.trial    <- numeric()  # Vector of the critical values in terms of hazard ratios for futility at the interim analysis of each trial within the series
  pvalue.efficacy.ia.trial    <- numeric()  # Vector of the critical values in terms of hazard ratios for efficacy at the interim analysis of each trial within the series
  pvalue.final.decision.trial <- numeric()  # Vector of the critical values in terms of hazard ratios for efficacy at the final analysis of each trial within the series
  ## Survival analysis
  nb.pts.C.s1.trial     <- numeric()  # Vector of the number of patients accrued in control group at the end of the first stage of each trial within the series
  nb.pts.E1.s1.trial    <- numeric()  # Vector of the number of patients accrued in experimental group 1 at the end of the first stage of each trial within the series
  nb.pts.E2.s1.trial    <- numeric()  # Vector of the number of patients accrued in experimental group 2 at the end of the first stage of each trial within the series
  nb.pts.C.s2.trial     <- numeric()  # Vector of the number of patients accrued in control group at the final analysis of each trial within the series
  nb.pts.E.s2.trial     <- numeric()  # Vector of the number of patients accrued in (carried over) experimental group at the final analysis of each trial within the series
  nb.death.C.s1.trial   <- numeric()  # Vector of the number of events observed in control group at the end of the first stage of each trial within the series
  nb.death.E1.s1.trial  <- numeric()  # Vector of the number of events observed in experimental group 1 at the end of the first stage of each trial within the series
  nb.death.E2.s1.trial  <- numeric()  # Vector of the number of events observed in experimental group 2 at the end of the first stage of each trial within the series
  nb.death.C.s2.trial   <- numeric()  # Vector of the number of events observed in control group at the final analysis of each trial within the series
  nb.death.E.s2.trial   <- numeric()  # Vector of the number of events observed in (carried over) experimental group at the final analysis of each trial within the series
  pts.year.C.s1.trial   <- numeric()  # Vector of the number of the sum of patient-years in control group at the end of the first stage of each trial within the series
  pts.year.E1.s1.trial  <- numeric()  # Vector of the number of the sum of patient-years in experimental group 1 at the end of the first stage of each trial within the series
  pts.year.E2.s1.trial  <- numeric()  # Vector of the number of the sum of patient-years in experimental group 2 at the end of the first stage of each trial within the series  
  pts.year.C.s2.trial   <- numeric()  # Vector of the number of the sum of patient-years in control group at the final analysis of each trial within the series
  pts.year.E.s2.trial   <- numeric()  # Vector of the number of the sum of patient-years in (carried over) experimental group at the final analysis of each trial within the series
  lambdaE1.trial        <- numeric()  # Vector of the experimental group 1 hazard rate of each trial within the series state of nature
  lambdaE2.trial        <- numeric()  # Vector of the experimental group 2 hazard rate of each trial within the series state of nature  
  lambdaC.obs.s1.trial  <- numeric()  # Vector of the control group hazard rate of each trial within the series estimated accorindg to simulation at the end of the first stage
  lambdaE1.obs.s1.trial <- numeric()  # Vector of the experimental group 1 hazard rate of each trial within the series estimated accorindg to simulation at the end of the first stage
  lambdaE2.obs.s1.trial <- numeric()  # Vector of the experimental group 2 hazard rate of each trial within the series estimated accorindg to simulation at the end of first stage
  lambdaC.obs.s2.trial  <- numeric()  # Vector of the control group hazard rate of each trial within the series estimated accorindg to simulation at the final analysis
  lambdaE.obs.s2.trial  <- numeric()  # Vector of the (carried over) experimental group hazard rate of each trial within the series estimated accorindg to simulation at the final analysis
  chisqE1.s1.trial      <- numeric()  # Vector of the log-rank test statistic at the end of the first stage experimental group 1 vs control   
  chisqE2.s1.trial      <- numeric()  # Vector of the log-rank test statistic at the end of the first stage experimental group 2 vs control
  chisqE.s2.trial       <- numeric()  # Vector of the log-rank test statistic at the final analysis experimental group 1 vs control 
  pvalue1.s1.trial      <- numeric()  # Vector of the p-value at the end of the first stage experimental group 1 vs control   
  pvalue2.s1.trial      <- numeric()  # Vector of the p-value at the end of the first stage experimental group 2 vs control
  pvalue12.s1.trial     <- numeric()  # Vector of the p-value at the end of the first stage of the intersection hypotheses      
  pvalue1.s2.trial      <- numeric()  # Vector of the p-value at the final analysis experimental group 1 vs control 
  pvalue2.s2.trial      <- numeric()  # Vector of the p-value at the final analysis experimental group 2 vs control
  pvalue12.s2.trial     <- numeric()  # Vector of the p-value at the final analysis of the intersection hypotheses  
  pvalue1.trial         <- numeric()  # Vector of the p-value at the final analysis experimental group 1 combined between stage1 and stage2
  pvalue2.trial         <- numeric()  # Vector of the p-value at the final analysis experimental group 2 combined between stage1 and stage2
  pvalue12.trial        <- numeric()  # Vector of the p-value at the final analysis of the intersection hypotheses combined between stage1 and stage2  
  hr.E1.trial           <- numeric()  # Vector of the hazard ratio of each trial within the series state of nature experimental group 1 vs control
  hr.E2.trial           <- numeric()  # Vector of the hazard ratio of each trial within the series state of nature experimental group 2 vs control 
  hr.E1.obs.s1.trial    <- numeric()  # Vector of the hazard ratio observed at the end of the first stage experimental group 1 vs control
  hr.E2.obs.s1.trial    <- numeric()  # Vector of the hazard ratio observed at the end of the first stage experimental group 2 vs control
  hr.E.obs.s2.trial     <- numeric()  # Vector of the hazard ratio observed at the final analysis (carried over) experimental group vs control
  ## First stage outcome
  end.s1.trial       <- numeric()  # Vector of whether the trial stopped at interim 
  time.end.trial     <- numeric()  # Vector of the time the trial ends 
  exp.retained.trial <- numeric()  # Vector of the experimental treatment TO BE carried on to the second stage
  ## Trial outcome
  weight.trial          <- numeric()  # Vector of the weight of the first stage, fixed a priori to 0.5 vs. 0.5
  end.s2.trial          <- numeric()  # Vector of the outcome at the end of the second stage (0=null hypothesis) 
  hr.retained.obs.trial <- numeric()  # Vector of the hazard ratio selected after each trial within the series estimated accorindg to simulation
  hr.retained.trial     <- numeric()  # Vector of the hazard ratio selected after each trial within the series state of nature
  
  # Define parameters for the mean of the experimental group hazard rate according to time since series initiation
  a      <- as.numeric(DistributionFutureTreatmentEffects(distrib, baseline)[1])
  b      <- as.numeric(DistributionFutureTreatmentEffects(distrib, baseline)[2])
  sigma  <- as.numeric(DistributionFutureTreatmentEffects(distrib, baseline)[3])
  
  set.seed(1234)
  
  for (i in 1:nbrepeat){  # ----- START LOOP 10 000 REPETITIONS -----
    time    <- 0  # Initialize time within the series 
    ntrial  <- 0  # Initialize the number of trials perfomred within the series
    
    lambdaC     <- baseline  # The drawn hazard rate of the baseline control arm 
    lambda.obs  <- baseline  # The observed hazard rate of the baseline control arm (considered the same for the first trial)
    
    # Whenever there is still time left for a trial before the end of the research horizon
    while (ifelse(is.na(time < Tmax), FALSE,  time < Tmax)){  # ----- START SERIES -----
      ## For sample size and IA time computation (Bonferroni corrected)
      try(design <- gsSurv(k = ns,
                           test.type = test.type, 
                           alpha = alpha/nexp,
                           beta = beta,
                           timing = timing,
                           sfl=sfExponential,
                           sflpar=0.7849295,
                           sfu=sfPoints,
                           sfupar = 0,  
                           lambdaC = lambda.obs,
                           hr = hr,
                           hr0 = hr0,
                           eta = 0,  # Dropout hazard rates for the control group
                           etaE = NULL,
                           gamma = accrual,
                           R = 1,  # Durations of time periods for recruitment rates 
                           S = NULL,
                           T = NULL,
                           minfup = FU,
                           ratio = ratio,
                           sided = sided))  # Try to design the trial which is not always possible
      
      n    <- ceiling((design$eNC + design$eNE)[2])  # Required sample size
      Texp <- design$T[2]  # Expected trial duration   
      Ts1  <- design$T[1]  # Expected first stage ending time (50% of the expected events) 
      Tacc <- Texp - FU  # Expected accrual period
      
      ## For bounds computation (Not Bonferroni corrected) but the trial duration is fixed as computed previously
      try(design <- gsSurv(k = ns,
                           test.type = test.type, 
                           alpha = alpha,
                           beta = beta,
                           timing = timing,
                           sfl=sfExponential,
                           sflpar=0.7849295,
                           sfu=sfPoints,
                           sfupar = 0,  
                           lambdaC = lambda.obs,
                           hr = hr,
                           hr0 = hr0,
                           eta = 0,  # Dropout hazard rates for the control group
                           etaE = NULL,
                           gamma = accrual,
                           R = 1,  # Durations of time periods for recruitment rates 
                           S = NULL,
                           T = Texp,
                           minfup = FU,
                           ratio = ratio,
                           sided = sided))  # Try to design the trial which is not always possible
      
      # p-value threshold for futility interim analysis
      pvalue.futility.ia    <- 1-pnorm(design$lower$bound)[1]
      # p-value threshold for efficacy interim analysis
      pvalue.efficacy.ia    <- 1-pnorm(design$upper$bound)[1]
      # p-value threshold for the final analysis
      pvalue.final.decision <- 1-pnorm(design$upper$bound)[2]
      
      if (abs(time + Texp - Tmax) >= abs(Tmax - time)) break  # Stop the series if the remaining time before the end of the horizon is shorter than the exceeding time if the trial would be initiated (considering that the IA will not result in ending the trial).
      
      # ----- START TRIAL -----
      ntrial <- ntrial + 1  # Increment the number of trials performed within the series if the algorithm doesn't break    
      
      ## Simulation parameters ID  
      accrual.rate.trial   <- c(accrual.rate.trial, accrual)
      lambda.trial         <- c(lambda.trial, baseline)
      FU.trial             <- c(FU.trial, FU)
      distribution.trial   <- c(distribution.trial, distrib)
      beta.trial           <- c(beta.trial, beta)
      hr.alternative.trial <- c(hr.alternative.trial, hr)
      ## Design parameters ID
      alpha.trial     <- c(alpha.trial, alpha)
      ia.type.trial   <- c(ia.type.trial, ia.type)
      nexp.trial      <- c(nexp.trial, nexp)
      ## Series (repetition) ID 
      series.trial <- c(series.trial, i)
      ## Trial ID
      ntrial.trial     <- c(ntrial.trial, ntrial)
      time.start.trial <- c(time.start.trial, time)
      ## Previuous trial legacy
      lambdaC.inherited.trial     <- c(lambdaC.inherited.trial, lambdaC)
      lambdaC.inherited.obs.trial <- c(lambdaC.inherited.obs.trial, lambda.obs) 
      ## Design characteristics
      n.trial    <- c(n.trial, n)
      Texp.trial <- c(Texp.trial, Texp)   
      Ts1.trial  <- c(Ts1.trial, Ts1)  
      Tacc.trial <- c(Tacc.trial, Tacc)
      ## Interim analysis bounds
      pvalue.futility.ia.trial     <- c(pvalue.futility.ia.trial, pvalue.futility.ia)
      pvalue.efficacy.ia.trial     <- c(pvalue.efficacy.ia.trial, pvalue.efficacy.ia)
      pvalue.final.decision.trial  <- c(pvalue.final.decision.trial, pvalue.final.decision)
      
      # ----- SIMULATION OF PATIENTS' SURVIVAL DATA -----        
      lambdaE      <- rlnorm(nexp, mu(time,a,b), sigma)  # Random sample of k hazard rates for k experimental treatments from the log-normal distribution, fixed variance and mean varying over time (from the series initiation) 
      lambda       <- c(lambdaC, lambdaE)
      time.accrual <- 1:n*Tacc/n  # Patients' accrual time assuming uniform accrual during the accrual duration
      ## State of nature
      lambdaE1.trial <- c(lambdaE1.trial, lambda[2])
      lambdaE2.trial <- c(lambdaE2.trial, lambda[3])         
      hr.E1.trial    <- c(hr.E1.trial, lambda[2]/lambda[1])
      hr.E2.trial    <- c(hr.E2.trial, lambda[3]/lambda[1]) 
      
      # ----- SIMULATION OF PATIENTS' SURVIVAL DATA FOR THE FIRST STAGE -----    
      time.accrual.s1    <- time.accrual[time.accrual<=Ts1]  # Accrual time of patients included in stage 1
      group.s1           <- gl(nexp+1,1,length(time.accrual.s1))  # Randomization group 1:1, 1: ctl, 2-3: exp         
      survival.time.s1   <- numeric(length(time.accrual.s1))  # Vector of patients' survival time (patients accrued during stage 1) 
      del.s1             <- numeric(length(time.accrual.s1))  # Vector of patients' survival time at the end of stage 1 (patients accrued during stage 1) 
      status.s1          <- numeric(length(time.accrual.s1))  # Vector of patients' status at the end of stage 1 (patients accrued during stage 1), censored if alive at the date of stage 1 analysis
      survival.data.s1   <- cbind.data.frame(time.accrual.s1, group.s1, survival.time.s1, del.s1, status.s1)  # Survival data of stage 1 (patients accrued during stage 1)
      
      for (k in 1:(nexp+1)){
        survival.data.s1$survival.time.s1[which(survival.data.s1$group.s1==k)]  <-
          rexp(NROW(survival.data.s1$survival.time.s1[which(survival.data.s1$group.s1==k)]),lambda[k])
      }
      
      survival.data.s1$del.s1    <- pmin(Ts1-survival.data.s1$time.accrual.s1, survival.data.s1$survival.time.s1)
      survival.data.s1$status.s1 <- ifelse(survival.data.s1$del.s1<(Ts1-survival.data.s1$time.accrual.s1),1,0)  # Event or censor
      
      nb.pts.C.s1  <- length(survival.data.s1[which(survival.data.s1$group.s1 == 1), "group.s1"])  # Number of patients accrued in the control group at the end of the first stage
      nb.pts.E1.s1 <- length(survival.data.s1[which(survival.data.s1$group.s1 == 2), "group.s1"])  # Number of patients accrued in experimental group 1 at the end of the first stage
      nb.pts.E2.s1 <- length(survival.data.s1[which(survival.data.s1$group.s1 == 3), "group.s1"])  # Number of patients accrued in experimental group 2 at the end of the first stage
      
      nb.death.C.s1  <- sum(survival.data.s1[which(survival.data.s1$group.s1 == 1), "status.s1"])  # Number of events observed in control group at the end of the first stage
      nb.death.E1.s1 <- sum(survival.data.s1[which(survival.data.s1$group.s1 == 2), "status.s1"])  # Number of events observed in experimental group 1 at the end of the first stage
      nb.death.E2.s1 <- sum(survival.data.s1[which(survival.data.s1$group.s1 == 3), "status.s1"])  # Number of events observed in experimental group 2 at the end of the first stage
      
      pts.year.C.s1  <- sum(survival.data.s1[which(survival.data.s1$group.s1 == 1), "del.s1"])  # Sum of patient-years in control group at the end of the first stage
      pts.year.E1.s1 <- sum(survival.data.s1[which(survival.data.s1$group.s1 == 2), "del.s1"])  # Sum of patient-years in experimental group 1 at the end of the first stage
      pts.year.E2.s1 <- sum(survival.data.s1[which(survival.data.s1$group.s1 == 3), "del.s1"])  # Sum of patient-years in experimental group 2 at the end of the first stage
      
      # ----- INFERENCE ON PATIENTS' SURVIVAL DATA FOR THE FIRST STAGE -----
      p1     <- numeric()
      chisq1 <- numeric()
      pp1     <- as.list(rep(NA,nexp))
      
      for (k in 1:nexp)
        length(pp1[[k]])     <- choose(nexp,k)
      
      for (k in 2:(nexp+1)){
        survival.data.E.s1          <- survival.data.s1[which(survival.data.s1$group.s1 %in% c(1,k)),] 
        survival.data.E.s1$group.s1 <- factor(survival.data.E.s1$group.s1)
        logrank.E                   <- survdiff(Surv(del.s1, status.s1) ~ group.s1, data=survival.data.E.s1)
        chisq1                      <- c(chisq1, logrank.E$chisq)
        p.logrank.E                 <- 1-pchisq(logrank.E$chisq, df=1)
        hazard                      <- exp(coxph(Surv(del.s1, status.s1) ~ group.s1, data=survival.data.E.s1)$coefficients)
        if (hazard >= 1)
          p1 <- c(p1, 1-p.logrank.E/2) else
            p1 <- c(p1, p.logrank.E/2)
      }
      
      pp1[[1]]     <- c(combn(p1,1))
      
      for (k in 2:nexp)
        pp1[[k]] <- SimesFirstStage(combn(p1,k))
      
      lambdaC.obs.s1  <- nb.death.C.s1/pts.year.C.s1 # the observed hazard rate in the control group, estimated from the simulated data at the end of the first stage
      lambdaE1.obs.s1 <- nb.death.E1.s1/pts.year.E1.s1 # the observed hazard rate in the experimental group 1, estimated from the simulated data at the end of the first stage   
      lambdaE2.obs.s1 <- nb.death.E2.s1/pts.year.E2.s1 # the observed hazard rate in the experimental group 2 , estimated from the simulated data at the end of the first stage
      
      chisqE1.s1 <- chisq1[1] 
      chisqE2.s1 <- chisq1[2]
      
      pvalue1.s1  <- pp1[[1]][1]
      pvalue2.s1  <- pp1[[1]][2]
      pvalue12.s1 <- pp1[[2]][1]
      
      hr.E1.obs.s1 <- lambdaE1.obs.s1/lambdaC.obs.s1
      hr.E2.obs.s1 <- lambdaE2.obs.s1/lambdaC.obs.s1
      
      ## Survival analysis
      nb.pts.C.s1.trial     <- c(nb.pts.C.s1.trial, nb.pts.C.s1)
      nb.pts.E1.s1.trial    <- c(nb.pts.E1.s1.trial, nb.pts.E1.s1)
      nb.pts.E2.s1.trial    <- c(nb.pts.E2.s1.trial, nb.pts.E2.s1)
      
      nb.death.C.s1.trial   <- c(nb.death.C.s1.trial, nb.death.C.s1)
      nb.death.E1.s1.trial  <- c(nb.death.E1.s1.trial, nb.death.E1.s1)
      nb.death.E2.s1.trial  <- c(nb.death.E2.s1.trial, nb.death.E2.s1)
      
      pts.year.C.s1.trial   <- c(pts.year.C.s1.trial, pts.year.C.s1)
      pts.year.E1.s1.trial  <- c(pts.year.E1.s1.trial, pts.year.E1.s1)
      pts.year.E2.s1.trial  <- c(pts.year.E2.s1.trial, pts.year.E2.s1)  
      
      lambdaC.obs.s1.trial  <- c(lambdaC.obs.s1.trial, lambdaC.obs.s1)
      lambdaE1.obs.s1.trial <- c(lambdaE1.obs.s1.trial, lambdaE1.obs.s1)
      lambdaE2.obs.s1.trial <- c(lambdaE2.obs.s1.trial, lambdaE2.obs.s1)
      
      chisqE1.s1.trial <- c(chisqE1.s1.trial, chisqE1.s1) 
      chisqE2.s1.trial <- c(chisqE2.s1.trial, chisqE2.s1)
      
      pvalue1.s1.trial   <- c(pvalue1.s1.trial, pvalue1.s1)   
      pvalue2.s1.trial   <- c(pvalue2.s1.trial, pvalue2.s1)
      pvalue12.s1.trial  <- c(pvalue12.s1.trial, pvalue12.s1)      
      
      hr.E1.obs.s1.trial <- c(hr.E1.obs.s1.trial, hr.E1.obs.s1)
      hr.E2.obs.s1.trial <- c(hr.E2.obs.s1.trial, hr.E2.obs.s1)
      
      ## First stage outcome
      p.s1           <- c(combn(p1,1))
      select.s1.ind  <- select.rule(-p.s1, type=1)$select # Selection of the smallest p-value
      
      end.s1         <- ifelse((((pvalue1.s1 >= pvalue.futility.ia) & (pvalue2.s1 >= pvalue.futility.ia)) |
                                  ((pvalue1.s1 >= pvalue.futility.ia) & (pvalue12.s1 >= pvalue.futility.ia)) |
                                  ((pvalue2.s1 >= pvalue.futility.ia) & (pvalue12.s1 >= pvalue.futility.ia)) |
                                  ((pvalue12.s1 >= pvalue.futility.ia))), 1, 0)  
      
      # If the trial stops at the end of the first stage, the time from the series initiation is incremented with the first stage ending time,
      # Unless the trial continue until the expected trial duration and the time from the series inition is incremented with the second stage ending time
      time <- ifelse(end.s1 == 1, time + Ts1, time + Texp)
      
      ## Interim analysis outcome
      end.s1.trial       <- c(end.s1.trial, end.s1)
      time.end.trial     <- c(time.end.trial, time)
      exp.retained.trial <- c(exp.retained.trial, which(select.s1.ind!=0))
      
      # ----- IF THE TRIAL STOPPED AT INTERIM (FOR FUTILITY)  -----        
      lambda.s1     <- lambdaC
      lambda.obs.s1 <- lambdaC.obs.s1
      
      # ----- SIMULATION OF PATIENTS' SURVIVAL DATA FOR THE SECOND STAGE -----    
      time.accrual.s2  <- time.accrual[time.accrual>Ts1]  # Accrual time of patients included in stage 2
      group.s2         <- gl(2,1,length(time.accrual.s2), c(1,which(select.s1.ind != 0)+1))  # Randomization group 1:1, 1: ctl, 2: exp 
      survival.time.s2 <- numeric(length(time.accrual.s2))  # Vector of patients' survival time (patients accrued during stage 2) 
      del.s2           <- numeric(length(time.accrual.s2))  # Vector of patients' survival time at the end of stage 2 (patients accrued during stage 2) 
      status.s2        <- numeric(length(time.accrual.s2))  # Vector of patients' status at the end of stage 2 (patients accrued during stage 2), censored if alive at the date of stage 2 analysis
      survival.data.s2 <- cbind.data.frame(time.accrual.s2, group.s2, survival.time.s2, del.s2, status.s2) # Survival data of stage 2 (patients accrued during stage 2)
      
      survival.data.s2$survival.time.s2[which(survival.data.s2$group.s2==1)]                            <-
        rexp(NROW(survival.data.s2$survival.time.s2[which(survival.data.s2$group.s2==1)]),lambdaC)
      survival.data.s2$survival.time.s2[which(survival.data.s2$group.s2==which(select.s1.ind != 0)+1)]  <-
        rexp(NROW(survival.data.s2$survival.time.s2[which(survival.data.s2$group.s2==which(select.s1.ind != 0)+1)]),lambdaE[which(select.s1.ind != 0)])
      
      survival.data.s2$del.s2    <- pmin(Texp-survival.data.s2$time.accrual.s2,survival.data.s2$survival.time.s2)
      survival.data.s2$status.s2 <- ifelse(survival.data.s2$del.s2<(Texp-survival.data.s2$time.accrual.s2),1,0)  # Event or censor       
      
      # ----- FIRST STAGE SURVIVAL DATA USED IN THE SECOND STAGE -----    
      survival.data.s1.s2 <- survival.data.s1[which(survival.data.s1$status.s1 != 1), ] #  Patients who did not experience an event at the end of the first stage        
      survival.data.s1.s2 <- survival.data.s1.s2[which(survival.data.s1.s2$group.s1 %in% c(1, which(select.s1.ind != 0)+1)), ] #   Patients allocated to the treatment carried over to the second stage or the control group
      survival.data.s1.s2$survival.time.s1 <- survival.data.s1.s2$survival.time.s1 - (Ts1 - survival.data.s1.s2$time.accrual.s1) #  Left truncation at the second stage    
      survival.data.s1.s2$del.s1           <- pmin(Texp-Ts1,survival.data.s1.s2$survival.time.s1)
      survival.data.s1.s2$status.s1        <- ifelse(survival.data.s1.s2$del.s1<Texp-Ts1,1,0)  # Event or censor 
      survival.data.s1.s2$time.accrual.s1  <- Ts1          
      
      colnames(survival.data.s1.s2) <- colnames(survival.data.s2)
      survival.data.s2              <- rbind(survival.data.s2, survival.data.s1.s2)          
      
      ## Survival analysis
      nb.pts.C.s2 <- length(survival.data.s2[which(survival.data.s2$group.s2 == 1), "group.s2"])  # Number of patients accrued in the control group at the end of the second stage
      nb.pts.E.s2 <- length(survival.data.s2[which(survival.data.s2$group.s2 != 1), "group.s2"])  # Number of patients accrued in the (carried over) experimental group at the end of the second stage
      
      nb.death.C.s2 <- sum(survival.data.s2[which(survival.data.s2$group.s2 == 1), "status.s2"])  # Number of events observed in the control group at the end of the second stage
      nb.death.E.s2 <- sum(survival.data.s2[which(survival.data.s2$group.s2 != 1), "status.s2"])  # Number of events observed in the (carried over) experimental group at the end of the second stage
      
      pts.year.C.s2 <- sum(survival.data.s2[which(survival.data.s2$group.s2 == 1), "del.s2"])  # Sum of patient-years in the control group at the end of the second stage 
      pts.year.E.s2 <- sum(survival.data.s2[which(survival.data.s2$group.s2 != 1), "del.s2"])  # Sum of patient-years in the (carried over) experimental group at the end of the second stage
      
      # ----- INFERENCE ON PATIENTS' SURVIVAL DATA FOR THE SECOND STAGE -----        
      p2 <- numeric()
      chisq2 <- numeric()
      pp2 <- as.list(rep(NA,nexp))
      
      for (k in 1:nexp)
        length(pp2[[k]]) <- choose(nexp,k)
      
      for (k in 2:(nexp+1)){
        survival.data.E.s2          <- survival.data.s2
        survival.data.E.s2$group.s2 <- factor(survival.data.E.s2$group.s2)
        logrank.E                   <- survdiff(Surv(del.s2, status.s2) ~ group.s2, data=survival.data.E.s2)
        chisq2                      <- c(chisq2, logrank.E$chisq)
        p.logrank.E                 <- 1-pchisq(logrank.E$chisq, df=1)
        hazard                      <- exp(coxph(Surv(del.s2, status.s2) ~ group.s2, data=survival.data.E.s2)$coefficients)
        if (hazard >= 1)
          p2 <- c(p2, 1-p.logrank.E/2) else
            p2 <- c(p2, p.logrank.E/2)
      }
      
      p2[which(select.s1.ind==0)] <- 1
      pp2[[1]] <- c(combn(p2,1))
      for (k in 2:nexp)
        pp2[[k]] <- SimesSecondStage(combn(p2,k))
      
      lambdaC.obs.s2 <- nb.death.C.s2/pts.year.C.s2 # the observed hazard rate in the control group, estimated from the simulated data at the end of the second stage
      lambdaE.obs.s2 <- nb.death.E.s2/pts.year.E.s2 # the observed hazard rate in the (carried over) experimental group, estimated from the simulated data at the end of the second stage
      
      chisqE.s2 <- chisq2[which(select.s1.ind!=0)]
      
      pvalue1.s2  <- pp2[[1]][1]
      pvalue2.s2  <- pp2[[1]][2]
      pvalue12.s2 <- pp2[[2]][1]
      
      hr.E.obs.s2 <- lambdaE.obs.s2/lambdaC.obs.s2
      
      ## Survival analysis
      nb.pts.C.s2.trial     <- c(nb.pts.C.s2.trial, nb.pts.C.s2)
      nb.pts.E.s2.trial     <- c(nb.pts.E.s2.trial, nb.pts.E.s2)
      
      nb.death.C.s2.trial   <- c(nb.death.C.s2.trial, nb.death.C.s2)
      nb.death.E.s2.trial   <- c(nb.death.E.s2.trial, nb.death.E.s2)
      
      pts.year.C.s2.trial   <- c(pts.year.C.s2.trial, pts.year.C.s2)
      pts.year.E.s2.trial   <- c(pts.year.E.s2.trial, pts.year.E.s2)
      
      lambdaC.obs.s2.trial  <- c(lambdaC.obs.s2.trial, lambdaC.obs.s2)
      lambdaE.obs.s2.trial  <- c(lambdaE.obs.s2.trial, lambdaE.obs.s2)
      
      chisqE.s2.trial <- c(chisqE.s2.trial, chisqE.s2)
      
      pvalue1.s2.trial      <- c(pvalue1.s2.trial, pvalue1.s2)   
      pvalue2.s2.trial      <- c(pvalue2.s2.trial, pvalue2.s2)
      pvalue12.s2.trial     <- c(pvalue12.s2.trial, pvalue12.s2)      
      
      hr.E.obs.s2.trial    <- c(hr.E.obs.s2.trial, hr.E.obs.s2)
      
      ## Trial outcome
      comb <- as.list(rep(NA,nexp))
      for (k in 1:nexp)
        length(comb[[k]]) <- choose(nexp,k)
      
      weight       <- 0.5 # The weights have to be predefined so as to control the FWER 
      weight.trial <- c(weight.trial, weight)
      
      for (ii in 1:nexp)
        for(jj in 1:length(comb[[ii]]))
          comb[[ii]][jj] <- InverseNormalCombination(pp1[[ii]][jj], pp2[[ii]][jj], weight)
      
      comb <- unlist(comb)
      
      pvalue1  <- comb[1]
      pvalue2  <- comb[2]
      pvalue12 <- comb[3]
      
      pvalue1.trial  <- c(pvalue1.trial, pvalue1)
      pvalue2.trial  <- c(pvalue2.trial, pvalue2)
      pvalue12.trial <- c(pvalue12.trial, pvalue12)  
      
      if (max(comb[which(comb != 1)]) < pvalue.final.decision){
        end.s2        <- 1
        lambda.s2     <- lambdaE[which(select.s1.ind != 0)]
        lambda.obs.s2 <- lambdaE.obs.s2         
      } else{
        end.s2        <- 0
        lambda.s2     <- lambdaC
        lambda.obs.s2 <- lambdaC.obs.s2
      }
      # ----- END OF THE CURRENT TRIAL -----
      # -----  ----- ----- Summary of the current trial -----  ----- -----
      # Selected Hazard rate according to the interim or final analysis (lambda.s1 and lambda.s2 could be lambdaC or lambda2)  
      end.s2.trial          <- c(end.s2.trial, end.s2)
      hr.retained.obs.trial <- c(hr.retained.obs.trial, ifelse(end.s1==1, lambda.obs.s1/lambdaC.obs.s1, lambda.obs.s2/lambdaC.obs.s2))
      hr.retained.trial     <- c(hr.retained.trial, ifelse(end.s1==1, lambda.s1/lambdaC, lambda.s2/lambdaC))
      # -----  ----- ----- Preparation of the next trial -----  ----- -----           
      lambdaC    <- ifelse(end.s1 == 1, lambda.s1, lambda.s2)  # The selected arm for the next trial as 'drawn'
      lambda.obs <- ifelse(end.s1 == 1, lambda.obs.s1, lambda.obs.s2)   # The selected arm for the next trial as observed 
      
    } # ----- END OF THE CURRENT SERIES -----
  } # ----- END OF SIMULATION FOR A SPECIFIC DESIGN ----- # ----- END LOOP 10 000 REPETITIONS -----
  results   <- cbind.data.frame(accrual.rate.trial,
                                lambda.trial,
                                FU.trial,
                                distribution.trial,
                                beta.trial,
                                hr.alternative.trial,
                                alpha.trial,
                                weight.trial,
                                nexp.trial,
                                ia.type.trial,
                                series.trial,
                                ntrial.trial,
                                time.start.trial,
                                lambdaC.inherited.trial,
                                lambdaC.inherited.obs.trial,
                                n.trial,
                                Texp.trial,
                                Ts1.trial,
                                Tacc.trial,
                                pvalue.futility.ia.trial,
                                pvalue.efficacy.ia.trial,
                                pvalue.final.decision.trial,
                                nb.pts.C.s1.trial,
                                nb.pts.E1.s1.trial,
                                nb.pts.E2.s1.trial,
                                nb.pts.C.s2.trial,
                                nb.pts.E.s2.trial,
                                nb.death.C.s1.trial,
                                nb.death.E1.s1.trial,
                                nb.death.E2.s1.trial,
                                nb.death.C.s2.trial,
                                nb.death.E.s2.trial,
                                pts.year.C.s1.trial,
                                pts.year.E1.s1.trial,
                                pts.year.E2.s1.trial,
                                pts.year.C.s2.trial,
                                pts.year.E.s2.trial,
                                lambdaE1.trial,
                                lambdaE2.trial,
                                lambdaC.obs.s1.trial,
                                lambdaE1.obs.s1.trial,
                                lambdaE2.obs.s1.trial,
                                lambdaC.obs.s2.trial,
                                lambdaE.obs.s2.trial,
                                chisqE1.s1.trial,
                                chisqE2.s1.trial,
                                chisqE.s2.trial,
                                pvalue1.s1.trial,
                                pvalue2.s1.trial,
                                pvalue12.s1.trial,
                                pvalue1.s2.trial,
                                pvalue2.s2.trial,
                                pvalue12.s2.trial,
                                pvalue1.trial,
                                pvalue2.trial,
                                pvalue12.trial,
                                hr.E1.trial,
                                hr.E2.trial,
                                hr.E1.obs.s1.trial,
                                hr.E2.obs.s1.trial,
                                hr.E.obs.s2.trial,
                                end.s1.trial,
                                time.end.trial,
                                end.s2.trial,
                                exp.retained.trial,
                                hr.retained.obs.trial,
                                hr.retained.trial)
  
  return(results)
} # ----- END OF SIMUL FUNCTION -----

results <- SimulateThreeArmOBFFutility(
  accrual  = 100,
  alpha    = 0.5,
  beta     = 0.1,
  hr       = 0.75,
  distrib  = 1,
  baseline = log(2),
  FU       = 0,
  nbrepeat = 10,
  Tmax     = 15
)

#' computing metrics 
final.results <- NULL
metric        <- data.frame()

##### Measure : Overall Hazard Ratio -----
overall.hazard.ratio <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          exp(sum(log([hr.retained.trial]))) as [overall.hazard.ratio]

          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- overall.hazard.ratio

##### Metric : Expected Overall Hazard Ratio -----
expected.overall.hazard.ratio <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
          
          avg([overall.hazard.ratio]) as [expected.overall.hazard.ratio]
  
          FROM measure
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial]"
)

metric <- expected.overall.hazard.ratio

##### Metric : Median Overall Hazard Ratio -----
metric$median.overall.hazard.ratio <- median(measure$overall.hazard.ratio, na.rm = TRUE) 

##### Metric : Expected Total Survival Benefit -----
metric$expected.total.survival.benefit <- mean(100*((1/measure$overall.hazard.ratio)-1), na.rm = TRUE)

##### Metric : Median Total Survival Benefit -----
metric$median.total.survival.benefit <- median(100*((1/measure$overall.hazard.ratio)-1), na.rm = TRUE)

##### Metric : Number of Repetition -----
metric$number.repeat <- nrow(measure)

##### Metric : Number of Detrimental Effect -----
metric$number.detrimental.effect <- length(which(measure$overall.hazard.ratio > 1))

##### Metric : Probability of a Detrimental Effect -----
metric$probability.detrimental.effect <- 100*metric$number.detrimental.effect/metric$number.repeat

##### Metric : Number of Achieving a Major Improvement -----
metric$number.achieving.major.improvement <- length(which(measure$overall.hazard.ratio < 0.7))

##### Metric : Probability of Achieving a Major Improvement -----
metric$probability.achieving.major.improvement <- 100*metric$number.achieving.major.improvement/metric$number.repeat

##### Metric : Number of clinically meaningful detriment to survival -----
metric$number.meaningful.detriment <- length(which(measure$overall.hazard.ratio > 1.1))

##### Metric : Probability of clinically meaningful detriment to survival -----
metric$probability.meaningful.detriment <- 100*metric$number.meaningful.detriment/metric$number.repeat

##### Results : Observed Number of Deaths -----
results$observed.number.deaths.trial <- ifelse(results$end.s1.trial==1,
                                               results$nb.death.C.s1.trial + results$nb.death.E1.s1.trial + results$nb.death.E2.s1.trial,
                                               results$nb.death.C.s1.trial + results$nb.death.E1.s1.trial + results$nb.death.E2.s1.trial + 
                                                 results$nb.death.C.s2.trial + results$nb.death.E.s2.trial
)

##### Measure : Observed Number of Deaths -----
observed.number.deaths <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],

          sum([observed.number.deaths.trial]) as [observed.number.deaths]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=observed.number.deaths,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Result : Expected Number of deaths if the patients had been treated with the control treatment of the first trial -----
for (i in 1:nrow(results)){
  n       <- results[i, "n.trial"]
  Tacc    <- results[i, "Tacc.trial"]
  Texp    <- results[i, "Texp.trial"]
  Ts1     <- results[i, "Ts1.trial"]
  lambdaC <- results[i, "lambda.trial"]
  
  # ----- SIMULATION OF PATIENTS' SURVIVAL DATA AS IF ALL PATIENTS HAD BEEN TREATED WITH THE INITIAL CONTROL TREATMENT -----
  time.accrual <- 1:n*Tacc/n  # Patients' accrual time assuming uniform accrual during the trial duration 
  # ----- SIMULATION OF PATIENTS' SURVIVAL DATA FOR THE FIRST STAGE -----    
  time.accrual.s1    <- time.accrual[time.accrual<=Ts1]  # Accrual time of patients included in stage 1
  group.s1           <- gl(3,1,length(time.accrual.s1))  # Randomization group 1:1, 1: ctl, 2-3: exp       
  survival.time.s1   <- numeric(length(time.accrual.s1))  # Vector of patients' survival time (patients accrued during stage 1) 
  del.s1             <- numeric(length(time.accrual.s1))  # Vector of patients' survival time at the end of stage 1 (patients accrued during stage 1) 
  status.s1          <- numeric(length(time.accrual.s1))  # Vector of patients' status at the end of stage 1 (patients accrued during stage 1), censored if alive at the date of stage 1 analysis
  survival.data.s1   <- cbind.data.frame(time.accrual.s1, group.s1, survival.time.s1, del.s1, status.s1)  # Survival data of stage 1 (patients accrued during stage 1)
  
  survival.data.s1$survival.time.s1 <- rexp(NROW(survival.data.s1),lambdaC)
  
  survival.data.s1$del.s1    <- pmin(Ts1-survival.data.s1$time.accrual.s1, survival.data.s1$survival.time.s1)
  survival.data.s1$status.s1 <- ifelse(survival.data.s1$del.s1<(Ts1-survival.data.s1$time.accrual.s1),1,0)  # Event or censor
  
  nb.death.s1  <- sum(survival.data.s1[, "status.s1"])  
  
  # ----- SIMULATION OF PATIENTS' SURVIVAL DATA FOR THE SECOND STAGE -----    
  time.accrual.s2  <- time.accrual[time.accrual>Ts1]  # Accrual time of patients included in stage 2
  group.s2         <- gl(2,1,length(time.accrual.s2))  # Randomization group 1:1, 1: ctl, 2: exp 
  survival.time.s2 <- numeric(length(time.accrual.s2))  # Vector of patients' survival time (patients accrued during stage 2) 
  del.s2           <- numeric(length(time.accrual.s2))  # Vector of patients' survival time at the end of stage 2 (patients accrued during stage 2) 
  status.s2        <- numeric(length(time.accrual.s2))  # Vector of patients' status at the end of stage 2 (patients accrued during stage 2), censored if alive at the date of stage 2 analysis
  survival.data.s2 <- cbind.data.frame(time.accrual.s2, group.s2, survival.time.s2, del.s2, status.s2) # Survival data of stage 2 (patients accrued during stage 2)
  
  survival.data.s2$survival.time.s2 <- rexp(NROW(survival.data.s2),lambdaC)
  
  survival.data.s2$del.s2    <- pmin(Texp-survival.data.s2$time.accrual.s2,survival.data.s2$survival.time.s2)
  survival.data.s2$status.s2 <- ifelse(survival.data.s2$del.s2<(Texp-survival.data.s2$time.accrual.s2),1,0)  # Event or censor       
  
  # ----- FIRST STAGE SURVIVAL DATA USED IN THE SECOND STAGE -----    
  survival.data.s1.s2 <- survival.data.s1[which(survival.data.s1$status.s1 != 1), ] #  Patients who did not experience an event at the end of the first stage        
  survival.data.s1.s2 <- survival.data.s1.s2[which(survival.data.s1.s2$group.s1 %in% c(1,2)), ] #   Patients allocated to the treatment carried over to the second stage or the control group (defined as 1 and 2 but it is alll the same)
  survival.data.s1.s2$survival.time.s1 <- survival.data.s1.s2$survival.time.s1 - (Ts1 - survival.data.s1.s2$time.accrual.s1) #  Left truncation at the second stage    
  survival.data.s1.s2$del.s1           <- pmin(Texp-Ts1,survival.data.s1.s2$survival.time.s1)
  survival.data.s1.s2$status.s1        <- ifelse(survival.data.s1.s2$del.s1<Texp-Ts1,1,0)  # Event or censor 
  survival.data.s1.s2$time.accrual.s1  <- Ts1          
  
  colnames(survival.data.s1.s2) <- colnames(survival.data.s2)
  survival.data.s2              <- rbind(survival.data.s2, survival.data.s1.s2)          
  
  nb.death.s2 <- sum(survival.data.s2[, "status.s2"])  # Number of events observed in the control group at the end of the second stage
  
  # ----- FINALLY -----    
  results$expected.number.deaths.trial[i] <- ifelse(results$end.s1.trial[i]==1, nb.death.s1, nb.death.s1 + nb.death.s2)               
}

##### Measure : Expected Number of deaths if the patients had been treated with the control treatment of the first trial -----
expected.number.deaths <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([expected.number.deaths.trial]) as [expected.number.deaths]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=expected.number.deaths,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : Proportion of Lives Saved -----
measure$proportion.lives.saved <- (measure$expected.number.deaths - measure$observed.number.deaths)/measure$expected.number.deaths  

##### Metric : Expected Proportion of Lives Saved -----
metric$expected.proportion.lives.saved <- mean(100*measure$proportion.lives.saved, na.rm = TRUE)

##### Metric : Number of Any Mortality Benefit -----
metric$number.mortality.benefit <- length(which(measure$proportion.lives.saved > 0))

##### Metric : Probability of Any Mortality Benefit -----
metric$probability.mortality.benefit <- 100*metric$number.mortality.benefit/metric$number.repeat

##### Measure : Number of Trials Performed -----
number.trials <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          max([ntrial.trial]) as [number.trials]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.trials,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Metric : Expected Number of Trials Performed -----
metric$expected.number.trials <- mean(measure$number.trials, na.rm = TRUE)

##### Measure : Research Horizon -----
research.horizon <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          max([time.end.trial]) as [research.horizon]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=research.horizon,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Metric : Expected Research Horizon -----
metric$expected.research.horizon <- mean(measure$research.horizon, na.rm = TRUE)

##### Measure : Number of Trials Stopped at Interim -----
number.trials.stopped.interim <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([end.s1.trial]) as [number.trials.stopped.interim]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.trials.stopped.interim,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : Probability of Stopping the Trial at Interim  -----
measure$probability.stopping.interim <- 100*(measure$number.trials.stopped.interim/measure$number.trials)  

##### Metric : Expected  Probability of Stopping the Trial at Interim -----
metric$expected.probability.stopping.interim <-  mean(measure$probability.stopping.interim, na.rm = TRUE)

##### Result : Boolean False Positive -----
results$best.E1.E2     <- ifelse(results$lambdaE1.trial<=results$lambdaE2.trial, 1, 2)
results$hr.best.E1.E2  <- ifelse(results$lambdaE1.trial<=results$lambdaE2.trial, results$hr.E1.trial, results$hr.E2.trial)
results$is.false.positive <- ifelse(results$hr.retained.trial != 1 & results$hr.best.E1.E2 >= 1, 1, 0)          
results$is.real.negative  <- ifelse(results$hr.best.E1.E2 >= 1, 1, 0)

##### Measure : Number of False Positive -----
number.false.positive <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([is.false.positive]) as [number.false.positive]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.false.positive,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : Number of Real Negative -----
number.real.negative <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([is.real.negative]) as [number.real.negative]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.real.negative,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : False Positive Rate  -----
measure$false.positive.rate <- 100*(measure$number.false.positive/measure$number.real.negative)  

##### Metric : Expected False Positive Rate -----
metric$expected.false.positive.rate <- mean(measure$false.positive.rate, na.rm=TRUE)

##### Result : Boolean False Negative -----
cutoff.positive            <- 0.8 
results$is.false.negative1 <- ifelse(((results$hr.retained.trial == 1) 
                                      | ((results$hr.retained.trial != 1) & (results$best.E1.E2 =! results$exp.retained.trial)))              
                                     & (results$hr.best.E1.E2 <= cutoff.positive), 1, 0) 
results$is.real.positive1  <- ifelse(results$hr.best.E1.E2 <= cutoff.positive, 1, 0)

percentage.positive        <- 0.8
results$is.false.negative2 <- ifelse(((results$hr.retained.trial == 1) 
                                      | ((results$hr.retained.trial != 1) & (results$best.E1.E2 =! results$exp.retained.trial)))              
                                     & (results$hr.best.E1.E2 <= results$hr.alternative.trial/percentage.positive), 1, 0) 
results$is.real.positive2  <- ifelse(results$hr.best.E1.E2 <= results$hr.alternative.trial/percentage.positive, 1, 0)

##### Measure : Number of False Negative -----
number.false.negative1 <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([is.false.negative1]) as [number.false.negative1]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.false.negative1,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : Number of Real Positive -----
number.real.positive1 <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([is.real.positive1]) as [number.real.positive1]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.real.positive1,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : False negative Rate  -----
measure$false.negative.rate1 <- 100*(measure$number.false.negative1/measure$number.real.positive1)  

##### Metric : Expected False negative Rate -----
metric$expected.false.negative.rate1 <- mean(measure$false.negative.rate1, na.rm = TRUE)

##### Measure : Number of False Negative -----
number.false.negative2 <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([is.false.negative2]) as [number.false.negative2]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.false.negative2,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : Number of Real Positive -----
number.real.positive2 <- sqldf("
          SELECT [accrual.rate.trial],
                 [lambda.trial],
                 [FU.trial],
                 [distribution.trial],
                 [beta.trial],
                 [hr.alternative.trial],
                 [alpha.trial],
                 [ia.type.trial],
                 [series.trial],
           
          sum([is.real.positive2]) as [number.real.positive2]
  
          FROM results
          GROUP BY [accrual.rate.trial],
                   [lambda.trial],
                   [FU.trial],
                   [distribution.trial],
                   [beta.trial],
                   [hr.alternative.trial],
                   [alpha.trial],
                   [ia.type.trial],
                   [series.trial]"
)

measure <- merge(x=measure, y=number.real.positive2,
                 by=c("accrual.rate.trial",
                      "lambda.trial",
                      "FU.trial",
                      "distribution.trial",
                      "beta.trial",
                      "hr.alternative.trial",
                      "alpha.trial",
                      "ia.type.trial",
                      "series.trial")
)

##### Measure : False negative Rate  -----
measure$false.negative.rate2 <- 100*(measure$number.false.negative2/measure$number.real.positive2)  

##### Metric : Expected False negative Rate -----
metric$expected.false.negative.rate2 <- mean(measure$false.negative.rate2, na.rm = TRUE)

##### Final step    
final.results <- rbind.data.frame(final.results, metric)
