#' @title Simulation of a series of two-arm trials with no interim analysis (fixed.sample).
#' @author Mohamed Amine BAYAR
#' \code{SimulateTwoArmNoInterimAnalysis} simulates a series of two-arm trials with no interim analysis (fixed.sample).
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
library(sqldf)
source("DistributionFutureTreatmentEffects.R")

# ----- DECLARE CONSTANCES -----
hr0        <- 1  # Hazard ratio experimental/control under the null hypothesis
ratio      <- 1  # Randomization ratio experimental/control
sided      <- 1  # 1-sided testing
ia.type    <- "None"

SimulateTwoArmNoInterimAnalysis  <- function(accrual, alpha, beta, hr, distrib, baseline, FU, nbrepeat, Tmax){  # ----- START OF SIMUL FUNCTION -----
  
  # ----- RESULTS STORAGE -----
  ## Simulation parameters ID  
  accrual.rate.trial   <- numeric()  # Vector of the accrual rate per year within the series
  lambda.trial         <- numeric()  # Vector of the baseline exponential survival distribution hazard rate in the control group within the series
  FU.trial             <- numeric()  # Vector of the fu within the series
  distribution.trial   <- numeric()  # Vector of the hypothetical distributions of futur treatments effect within the series
  beta.trial           <- numeric()  # Vector of the type II error within the series
  hr.alternative.trial <- numeric()  # Vector of the hazard ratio of alternative hypothsesis within the series
  ## Design parameters ID
  alpha.trial          <- numeric()  # Vector of the type I error within the series
  ia.type.trial        <- numeric()  # Vector of the type of interim analysis within the series
  ## Series (repetition) ID 
  series.trial <- numeric()  # Vector of the series ID
  ## Trial ID
  ntrial.trial     <-  numeric()  # Vector of the trial ID
  time.start.trial <-  numeric()  # Vector of the time the trial starts
  ## Previuous trial legacy
  lambdaC.inherited.trial     <- numeric()  # Vector of the control group hazard rate (inherited from the previous trial) of each trial within the series state of nature
  lambdaC.inherited.obs.trial <- numeric()  # Vector of the control group hazard rate (inherited from the previous trial) of each trial within the series estimated accorindg to simulation (in the previous trial) 
  ## Design characteristics
  n.trial                 <- numeric()  # Vector of the required sample size of each trial within the series
  Texp.trial              <- numeric()  # Vector of the expected trial duration of each trial within the series   
  Tacc.trial              <- numeric()  # Vector of the expected accrual duration of each trial within the series
  ## Survival analysis
  nb.pts.C.final.trial    <- numeric()  # Vector of the number of patients accrued in control group at the final analysis of each trial within the series
  nb.pts.E.final.trial    <- numeric()  # Vector of the number of patients accrued in experimental group at the final analysis of each trial within the series
  nb.death.C.final.trial  <- numeric()  # Vector of the number of events observed in control group at the final analysis of each trial within the series
  nb.death.E.final.trial  <- numeric()  # Vector of the number of events observed in experimental group at the final analysis of each trial within the series
  pts.year.C.final.trial  <- numeric()  # Vector of the number of the sum of patient-years in control group at the final analysis of each trial within the series
  pts.year.E.final.trial  <- numeric()  # Vector of the number of the sum of patient-years in experimental group at the final analysis of each trial within the series
  lambdaE.trial           <- numeric()  # Vector of the experimental group hazard rate of each trial within the series state of nature
  lambdaC.obs.final.trial <- numeric()  # Vector of the control group hazard rate of each trial within the series estimated accorindg to simulation at the final analysis
  lambdaE.obs.final.trial <- numeric()  # Vector of the experimental group hazard rate of each trial within the series estimated accorindg to simulation at the final analysis
  chisq.final.trial       <- numeric()  # Vector of the log-rank test statistic at the final analysis
  hr.trial                <- numeric()  # Vector of the hazard ratio of each trial within the series state of nature
  hr.obs.final.trial      <- numeric()  # Vector of the hazard ratio observed at the final analysis
  pvalue.final.trial      <- numeric()  # Vector of the log-rank test p-value at the final analysis
  ## Trial outcome
  time.end.trial        <- numeric()  # Vector of the time the trial ends 
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
    
    lambdaC     <- baseline  # The drawn hazard rate of the baseline control group 
    lambda.obs  <- baseline  # The observed hazard rate of the baseline control group (considered the same for the first trial)
    
    # Whenever there is still time left for a trial before the end of the research horizon
    while (ifelse(is.na(time < Tmax), FALSE,  time < Tmax)){  # ----- START SERIES -----
      try(design <- nSurv(alpha = alpha,
                          beta = beta,
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
      
      n     <- ceiling(design$n)  # Required sample size
      Texp  <- design$T           # Expected trial duration   
      Tacc  <- Texp - FU        # Expected accrual duration 
      
      if (abs(time + Texp - Tmax) >= abs(Tmax - time)) break  # Stop the series if the remaining time before the end of the horizon is shorter than the exceeding time if the trial would be initiated (considering that the the interim analysiswill not result in ending the trial).
      
      # ----- START TRIAL -----
      ntrial   <- ntrial + 1  # Increment the number of trials performed within the series if the algorithm doesn't break
      
      ## Simulation parameters ID  
      accrual.rate.trial   <- c(accrual.rate.trial, accrual)
      lambda.trial         <- c(lambda.trial, baseline)
      FU.trial             <- c(FU.trial, FU)
      distribution.trial   <- c(distribution.trial, distrib)
      beta.trial           <- c(beta.trial, beta)
      hr.alternative.trial <- c(hr.alternative.trial, hr)
      ## Design parameters ID
      alpha.trial   <- c(alpha.trial, alpha)
      ia.type.trial <- c(ia.type.trial, ia.type)
      ## Series (repetition) ID 
      series.trial <- c(series.trial, i)
      ## Trial ID
      ntrial.trial     <- c(ntrial.trial, ntrial)
      time.start.trial <- c(time.start.trial, time) 
      ## Previuous trial legacy
      lambdaC.inherited.trial     <- c(lambdaC.inherited.trial, lambdaC)
      lambdaC.inherited.obs.trial <- c(lambdaC.inherited.obs.trial, lambda.obs) 
      ## Design characteristics
      n.trial                 <- c(n.trial, n)
      Texp.trial              <- c(Texp.trial, Texp)   
      Tacc.trial              <- c(Tacc.trial, Tacc)
      
      # ----- SIMULATION OF PATIENTS' SURVIVAL DATA -----
      del.final    <- numeric(n)  # Vector of patients' survival data in the trial (??? participation time)  
      status.final <- numeric(n)  # Vector of patients' status at last fu, censored if alive at the date of last fu
      lambdaE      <- rlnorm(1, mu(time,a,b), sigma)  # Random sample from the log-normal distribution of the experimental group hazard rate, fixed varaiance and mean varing over time (from the series initiation)
      time.accrual <- 1:n*Tacc/n  # Patients' accrual time assuming uniform accrual during the trial duration 
      group        <- gl(2,1,n)   # Randomization group 1:1, 1: control, 2: experimental 
      
      # -----  ----- FOR THE FINAL ANALYSIS -----  -----  
      del.final[which(group==1)] <- pmin(Texp-time.accrual,rexp(n,lambdaC))[which(group==1)]  # Time to censor or time to event in the control group 
      del.final[which(group==2)] <- pmin(Texp-time.accrual,rexp(n,lambdaE))[which(group==2)]  # Time to censor or time to event in the experimental group
      status.final               <- ifelse(del.final<(Texp-time.accrual),1,0)  # Event or censor 
      survival.data.final        <- data.frame(cbind(group, del.final, status.final, time.accrual))  # Final analysis survival data
      
      nb.pts.C.final <- length(survival.data.final[which(survival.data.final$group == 1), "group"])  # Number of patients accrued at the final analysis in the control group
      nb.pts.E.final <- length(survival.data.final[which(survival.data.final$group == 2), "group"])  # Number of patients accrued at the final analysis in the experimental group
      
      nb.death.C.final <- sum(survival.data.final[which(survival.data.final$group == 1), "status.final"])  # Number of events at the final analysis in the control group
      nb.death.E.final <- sum(survival.data.final[which(survival.data.final$group == 2), "status.final"])  # Number of events at the final analysis in the experimental group
      
      pts.year.C.final <- sum(survival.data.final[which(survival.data.final$group == 1), "del.final"])  # Sum of patient-years at the final analysis in the control group
      pts.year.E.final <- sum(survival.data.final[which(survival.data.final$group == 2), "del.final"])  # Sum of patient-years at the final analysis in the experimental group
      
      # ----- INFERENCE ON PATIENTS' SURVIVAL DATA -----
      # -----  ----- FOR THE FINAL ANALYSIS -----  -----
      chisq.final        <- survdiff(Surv(del.final,status.final) ~ group, data=survival.data.final)$chisq  # Compute the log-rank test statistic at the final analysis
      hr.obs.final       <- exp(coxph(Surv(del.final,status.final) ~ group, data=survival.data.final)$coefficients)  # Compute the hazard ratio observed at the final analysis
      pvalue.logrank.final    <- 1-pchisq(chisq.final, df=1)
      if (hr.obs.final >= 1)
        pvalue.final <- 1-pvalue.logrank.final/2 else
          pvalue.final <- pvalue.logrank.final/2
      lambdaC.obs.final  <- nb.death.C.final/pts.year.C.final  # the observed hazard rate in the control group, estimated from the sumilated data at the final analysis
      lambdaE.obs.final  <- nb.death.E.final/pts.year.E.final  # the observed hazard rate in the experimental group, estimated from the sumilated data at the final analysis
      
      ## Survival analysis
      nb.pts.C.final.trial    <- c(nb.pts.C.final.trial, nb.pts.C.final)
      nb.pts.E.final.trial    <- c(nb.pts.E.final.trial, nb.pts.E.final)
      nb.death.C.final.trial  <- c(nb.death.C.final.trial, nb.death.C.final)
      nb.death.E.final.trial  <- c(nb.death.E.final.trial, nb.death.E.final)
      pts.year.C.final.trial  <- c(pts.year.C.final.trial, pts.year.C.final)
      pts.year.E.final.trial  <- c(pts.year.E.final.trial, pts.year.E.final)
      lambdaE.trial           <- c(lambdaE.trial, lambdaE)
      lambdaC.obs.final.trial <- c(lambdaC.obs.final.trial, lambdaC.obs.final)
      lambdaE.obs.final.trial <- c(lambdaE.obs.final.trial, lambdaE.obs.final)
      chisq.final.trial       <- c(chisq.final.trial, chisq.final)
      hr.trial                <- c(hr.trial, lambdaE/lambdaC) 
      hr.obs.final.trial      <- c(hr.obs.final.trial, hr.obs.final)
      pvalue.final.trial      <- c(pvalue.final.trial, pvalue.final)
      
      # ----- HYPOTHESES TESTING  -----
      # -----  ----- FOR THE INTERIM ANALYSIS -----  -----
      # The trial never stops at the interim analysis, the trial continue until the expected trial duation.
      # The time from the series inition is incremented with the final analyis time
      time <- time + Texp
      time.end.trial <- c(time.end.trial, time) 
      
      # -----  ----- FOR THE FINAL ANALYSIS -----  -----
      # According to the log-rank test p-value retain control or experimental for the next trial and test of HR
      lambda.obs.final <- ifelse(pvalue.final <= alpha, lambdaE.obs.final, lambdaC.obs.final)
      lambda.final     <- ifelse(pvalue.final <= alpha, lambdaE, lambdaC)
      # ----- END OF THE CURRENT TRIAL -----
      
      # -----  ----- ----- Summary of the current trial -----  ----- -----
      hr.retained.trial     <- c(hr.retained.trial, lambda.final/lambdaC)  # State of nature 
      hr.retained.obs.trial <- c(hr.retained.obs.trial, lambda.obs.final/lambdaC.obs.final)  # according to simulation 
      # -----  ----- ----- Preparation of the next trial -----  ----- -----
      lambdaC             <- lambda.final  # The selected group for the next trial as 'drawn'
      lambda.obs          <- lambda.obs.final   # The selected group for the next trial as observed 
      
    } # ----- END OF THE CURRENT SERIES -----
  } # ----- END OF SIMULATION FOR A SPECIFIC DESIGN ----- # ----- END LOOP 10 000 REPETITIONS -----
  results   <- cbind.data.frame(accrual.rate.trial,
                                lambda.trial,
                                FU.trial,
                                distribution.trial,
                                beta.trial,
                                hr.alternative.trial,
                                alpha.trial,
                                ia.type.trial,
                                series.trial,
                                ntrial.trial,
                                time.start.trial,
                                lambdaC.inherited.trial,
                                lambdaC.inherited.obs.trial,
                                n.trial,
                                Texp.trial,
                                Tacc.trial,
                                nb.pts.C.final.trial,
                                nb.pts.E.final.trial,
                                nb.death.C.final.trial,
                                nb.death.E.final.trial,
                                pts.year.C.final.trial,
                                pts.year.E.final.trial,
                                lambdaE.trial,
                                lambdaC.obs.final.trial,
                                lambdaE.obs.final.trial,
                                chisq.final.trial,
                                hr.trial,
                                hr.obs.final.trial,
                                pvalue.final.trial,
                                time.end.trial,
                                hr.retained.obs.trial,
                                hr.retained.trial)
  
  return(results)
} # ----- END OF SIMUL FUNCTION -----

results <- SimulateTwoArmNoInterimAnalysis(
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
results$observed.number.deaths.trial <- results$nb.death.C.final.trial + results$nb.death.E.final.trial

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
  Tia     <- results[i, "Tia.trial"]
  lambdaC <- results[i, "lambda.trial"]
  
  # ----- SIMULATION OF PATIENTS' SURVIVAL DATA -----
  del.final    <- numeric(n)  # Vector of patients' survival data in the trial (??? participation time)  
  status.final <- numeric(n)  # Vector of patients' status at last fu, censored if alive at the date of last fu
  time.accrual <- 1:n*Tacc/n  # Patients' accrual time assuming uniform accrual during the trial duration 
  # -----  ----- FOR THE FINAL ANALYSIS -----  -----  
  del.final           <- pmin(Texp-time.accrual,rexp(n,lambdaC))  # Time to censor or time to event in the control group 
  status.final        <- ifelse(del.final<(Texp-time.accrual),1,0)  # Event or censor 
  survival.data.final <- data.frame(cbind(del.final, status.final, time.accrual))  # Final analysis survival data
  nb.death.final      <- sum(survival.data.final[, "status.final"])  # Number of events at the final analysis in the control group
  results$expected.number.deaths.trial[i] <- nb.death.final               
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
                                       
                                       sum(0) as [number.trials.stopped.interim]
                                       
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
#results$is.false.positive <- ifelse(results$hr.retained.trial > 1, 1, 0)
results$is.false.positive <- ifelse(results$hr.retained.trial != 1 & results$hr.trial >= 1, 1, 0)          
results$is.real.negative  <- ifelse(results$hr.trial >= 1, 1, 0)

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
results$is.false.negative1 <- ifelse(results$hr.retained.trial == 1 & results$hr.trial <= cutoff.positive, 1, 0) 
results$is.real.positive1  <- ifelse(results$hr.trial <= cutoff.positive, 1, 0)

percentage.positive        <- 0.8
results$is.false.negative2 <- ifelse(results$hr.retained.trial == 1 & results$hr.trial <= results$hr.alternative.trial/percentage.positive, 1, 0)
results$is.real.positive2  <- ifelse(results$hr.trial <= results$hr.alternative.trial/percentage.positive, 1, 0)

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
