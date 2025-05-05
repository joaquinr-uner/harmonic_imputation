#Load required R.matlab package

library(R.matlab)
library(mvLSWimpute)
library(MLmetrics)

RunLSWimpute <- function(name) {
  #Load data saved from Matlab
  #source('temp/MissingSamples.R')
  
  source(name)
  
  #Run one sample t test and save results
  
  s_ts = ts(s,start = c(1))
  
  s_ts[miss_ints, ] <- NA
  newdata<- mv_impute(s_ts, index = miss_ints, p=steps, type = "forward-backward") 
  #Save results to Matlab data file
  
  imputed_values = newdata
  
  #writeMat('temp/Imputevalues.mat',imp=imputed_values)
  writeMat(paste(name,'.mat',sep=''),imp=imputed_values$ImputedData)
}
args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
RunLSWimpute(name)