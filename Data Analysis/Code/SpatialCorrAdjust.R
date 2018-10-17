setwd("/vol/data/zhuz/haozhe/Zillow/")

library(parallel)

source("earth_dist.R")
source("EstMean_Fun.R")
source("Est3DCov_Fun.R")
source("Est2DCov_Fun.R")
source("EstNugget_Fun.R")
source("fpca_Fun.R")
source("Est2DCovNugget_Fun.R")
source("EstSpatialCorr_Fun.R")
source("Est2DNugget_Fun.R")
source("EstFPC_Fun.R")

numCores <- detectCores()

Est3DCov_output <- readRDS("output/Est3DCov_output_zillow_0927.rds")
Est2DCov_output <- Est2DCov_Fun(Est3DCov_output, 0, 1)
fpca_output <- fpca_Fun(Est2DCov_output)
EstSpatialCorr_output <- EstSpatialCorr_Fun(Est3DCov_output, fpca_output)

######## Positive semidefinite adjustment for the second spatial correlation function
SpatialCorr_Hankel_Fun <- function(theta){
  tmp <- function(u){
    EstSpatialCorr_output$SpatialCorr_Fun(u,2)*besselJ(theta*u, 0)*u
  }
  tmp <- Vectorize(tmp)
  integrate(tmp, 0, 3)$value
}
SpatialCorr_positivedefinite_Fun <- function(u){
  tmp1 <- function(theta){
    max(SpatialCorr_Hankel_Fun(theta),0)*besselJ(theta*u, 0)*theta
  }
  tmp1 <- Vectorize(tmp1)
  #mean(sapply(seq(0,10, 0.05), FUN = tmp))*10
  integrate(tmp1, 0, 10)$value
}
SpatialCorr_positivedefinite_Fun <- Vectorize(SpatialCorr_positivedefinite_Fun)

distance <- seq(0,3.5,0.1)
adjust <- mclapply(distance, SpatialCorr_positivedefinite_Fun, mc.cores = numCores - 2)
print(adjust)

saveRDS(list(distance=distance, adjust=adjust), file="output/Adjusted_SpatialCorr_2.rds")
