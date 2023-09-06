#Load required R.matlab package

library(R.matlab)
library(forecast)
library(MLmetrics)

RunTBATS <- function(name) {
#Load data saved from Matlab
  #source('temp/MissingSamples.R')

  source(name)

#Run one sample t test and save results

  si_ts = ts(si,start = c(1), frequency=T)
  tbatsfit <- tbats(si_ts, use.parallel = FALSE)

  tbatsfcast<- forecast(tbatsfit,h=Li) 
  #Save results to Matlab data file

  imputed_values = tbatsfcast$mean

  #writeMat('temp/Imputevalues.mat',imp=imputed_values)
  writeMat(paste(name,'.mat',sep=''),imp=imputed_values)
}
args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
RunTBATS(name)