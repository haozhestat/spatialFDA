#library(RandomFields)
library(data.table)
library(reshape2)
library(splines)
library(splines2)
library(fda)
library(ggplot2)
library(fdapace)
library(geoR)
library(fdagstat)
#library(latex2exp)
#library(ggplot2)
library (plyr)

#setwd("~/Dropbox/projects/SPCA/Zillow/")
setwd("/vol/data/zhuz/haozhe/Zillow/")
source("earth_dist.R")
source("EstMean_Fun.R")
source("Est2DCov_Fun.R")
source("EstFPC_Fun.R")
source("EstSpatialCov_Fun.R")
source("Est2DCovYao_Fun.R")
source("Est2DNugget_Fun.R")
source("EstMeasurementError_Fun.R")

zillow_data <- na.omit(readRDS("zillow_data.rds"))
location <- readRDS("location.rds")
location <- location[, c("lon", "lat")]

zillow_data <- subset(zillow_data, t_coord>0)

n = max(zillow_data$location_id)
m = max(zillow_data$t_coord) - min(zillow_data$t_coord) + 1
K_m = 10
K_t = 4
K_s = 3
deg_t = 3
deg_s = 3
t_limit = c(min(zillow_data$t_coord), max(zillow_data$t_coord))
Delta = 3
delta = 2.5
delta_lower = 0.2
delta_upper = 1

county_list <- unique(zillow_data$CountyName)
mean_est <- NULL
for(i in 1:length(county_list)){
  zillow_data$value_extmean[zillow_data$CountyName==county_list[i]] <-
    zillow_data$value[zillow_data$CountyName==county_list[i]] -
    EstMean_Fun(data = zillow_data[zillow_data$CountyName==county_list[i],],
                K_m = K_m, t_limit = t_limit)$EstMean(zillow_data$t_coord[zillow_data$CountyName==county_list[i]])
  mean_est <- rbind(mean_est,
                    data.frame(t = seq(t_limit[1], t_limit[2],0.1),
                               value = EstMean_Fun(data = zillow_data[zillow_data$CountyName==county_list[i],],
                                                   K_m = K_m,
                                                   t_limit = t_limit)$EstMean(seq(t_limit[1],t_limit[2],0.1)),
                               CountyName=county_list[i]))
}

Flies <- MakeFPCAInputs(zillow_data$location_id, 
                        zillow_data$t_coord, 
                        zillow_data$value_extmean)

fpcaObjFlies <- FPCA(Flies$Ly, Flies$Lt,
                     list(methodMuCovEst = 'smooth',
                          methodSelectK = 3))

pace_pred_error <- rep(NA, nrow(location))
for(j in 1:nrow(location)){
  g <- fstat(NULL, 'FPC1', location[-j,], as.data.frame(fpcaObjFlies$xiEst[-j,1]), scalar=TRUE)
  g <- fstat(g, 'FPC2', location[-j,], as.data.frame(fpcaObjFlies$xiEst[-j,2]), scalar=TRUE)
  g <- estimateDrift("~.", g, Intercept = TRUE)
  g <- fvariogram("~.", g, Nlags = 100, LagMax = 1, useResidual = TRUE, comments=FALSE)
  g <- fitVariograms(g, model = vgm('Mat'), forceNugget = TRUE, fitRanges = FALSE)
  g <- addCovariance(g, 'omni')

  forecast.pc1 <- predictFstat(g, .newCoordinates = location[j,], .what = "FPC1", .type = "UcoK", algIndependent = TRUE)
  forecast.pc2 <- predictFstat(g, .newCoordinates = location[j,], .what = "FPC2", .type = "UcoK", algIndependent = TRUE)

  true_curve <- zillow_data$value_extmean[zillow_data$location_id == j]
  
  pace_pred_curve <- c(forecast.pc1$Forecast)*fpcaObjFlies$phi[,1]+
    c(forecast.pc2$Forecast)*fpcaObjFlies$phi[,2]
  
  pace_pred_error[j] <- mean((true_curve - pace_pred_curve)^2)
  if(length(true_curve)!=length(pace_pred_curve)){
    pace_pred_error[j] <- mean((true_curve - pace_pred_curve[zillow_data$t_coord[zillow_data$location_id == j]])^2)
  }
  print(pace_pred_error[j])
}

print(pace_pred_error[unique(zillow_data$location_id[zillow_data$CountyName == "San Francisco"]]))









