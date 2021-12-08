# set seed
rm(list=ls())
source("Stochastic control model/storage model.R")

set.seed(10)

# define model
StorageModel <- list(
  # storage characteristics
  capa.min = 0, # minimum capacity (MWh)
  capa.max = 100, # maximum capacity (MWh)
  capa.steps = 100, # capacity discretisation (steps)
  charging.rate = 10, # rate of injection MWh/time step
  discharging.rate = 10, # rate of withdrawal MWh/time step
  charging.efficiency = 0.9, # efficiency when charging (%)
  discharging.efficiency = 0.9, #efficiency when discharging (%)
  
  # model framework
  T = 1, # total periods in years or months
  dt = 1/30, # length of time step as fraction of period
  r = 0.05, # interest rate
  spread = function(x){0.02*x}, # function describing bid-ask spread
  
  
  # simulation parameters
  pilot.nsims = 1000, # pilot simulations see documentation of mlOSP
  sim.func = sim.gbm, # dynamics of price
  sigma = 0.3, # volatility parameter to GBM
  div = 0, #  dividends
  dim = 1, # number of exogenous factors
  x0 = 40, # starting price
  N = 800, # number of unique simulations
  batch.nrep =25, # number of replications
  nk = 16, # number of splines
  
  placeholder = NULL 
)

test <- storage.fixed.design(StorageModel,input.domain = 0.04,method = "spline")


