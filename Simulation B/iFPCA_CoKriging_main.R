#install.packages("fdapace")
library(fdapace)
library(geoR)
library(fdagstat)

# Turn the original data into a list of paired amplitude and timing lists
setwd("/vol/data/zhuz/haozhe/simulation2")
#setwd("~/Dropbox/projects/SPCA/simulation2")

source("Simulation_Setting.R")

library(doParallel)
library(foreach)

cl <- makeCluster(10)
registerDoParallel(cl)

foreach(i=1:200, .combine=list,.packages=c("fdapace","geoR","fdagstat")) %dopar% {
  sim_data <- read.csv(paste0("sim_data/sim_data_", i, ".csv"))
  sim_data_location <- read.csv(paste0("sim_data_location/sim_data_location_", i, ".csv"))
  sim_data_fpcscore <- read.csv(paste0("sim_data_fpcscore/sim_data_fpcscore_", i, ".csv"))
  
  test <- sample(1:nrow(sim_data_location), 100)
  train <- setdiff(1:nrow(sim_data_location), test)
  
  #sim_data <- sim_data[sim_data$location_id<300,]
  sim_data <- sim_data[order(sim_data$t_coord),]
  Flies <- MakeFPCAInputs(sim_data$location_id, sim_data$t_coord, sim_data$value)
  
  fpcaObjFlies <- FPCA(Flies$Ly, Flies$Lt,
                       list(methodMuCovEst = 'smooth',
                            methodSelectK = 3))
  
  g <- fstat(NULL, 'FPC1', sim_data_location[train,], as.data.frame(fpcaObjFlies$xiEst[train,1]), scalar=TRUE)
  g <- fstat(g, 'FPC2', sim_data_location[train,], as.data.frame(fpcaObjFlies$xiEst[train,2]), scalar=TRUE)
  g <- fstat(g, 'FPC3', sim_data_location[train,], as.data.frame(fpcaObjFlies$xiEst[train,3]), scalar=TRUE)
  #g <- fstat(g, 'FPC4', sim_data_location[train,], as.data.frame(fpcaObjFlies$xiEst[train,4]), scalar=TRUE)
  #g <- fstat(g, 'FPC5', sim_data_location[train,], as.data.frame(fpcaObjFlies$xiEst[train,5]), scalar=TRUE)
  g <- estimateDrift("~.", g, Intercept = TRUE)
  g <- fvariogram("~.", g, Nlags = 100, LagMax = 10, useResidual = TRUE, comments=FALSE)
  g <- fitVariograms(g, model = vgm('Mat'), forceNugget = TRUE, fitRanges = FALSE)
  g <- addCovariance(g, 'omni')
    
  forecast.pc1 <- predictFstat(g, .newCoordinates = sim_data_location[test,], .what = "FPC1", .type = "UcoK", algIndependent = TRUE)
  forecast.pc2 <- predictFstat(g, .newCoordinates = sim_data_location[test,], .what = "FPC2", .type = "UcoK", algIndependent = TRUE)
  forecast.pc3 <- predictFstat(g, .newCoordinates = sim_data_location[test,], .what = "FPC3", .type = "UcoK", algIndependent = TRUE)
  #forecast.pc4 <- predictFstat(g, .newCoordinates = sim_data_location[test,], .what = "FPC4", .type = "UcoK", algIndependent = TRUE)
  #forecast.pc5 <- predictFstat(g, .newCoordinates = sim_data_location[test,], .what = "FPC5", .type = "UcoK", algIndependent = TRUE)
  
  Grid <- fpcaObjFlies$workGrid
  pred_error <- rep(NA, length(test))
  for(k in 1:length(test)){
    j = test[k]
    
    true_curve <- mean_fun(Grid) +
      phi1_fun(Grid)*sim_data_fpcscore$score1[j]+
      phi2_fun(Grid)*sim_data_fpcscore$score2[j]+
      phi3_fun(Grid)*sim_data_fpcscore$score3[j]
    
    pred_curve <- fpcaObjFlies$mu +
      forecast.pc1$Forecast[k]*fpcaObjFlies$phi[,1]+
      forecast.pc2$Forecast[k]*fpcaObjFlies$phi[,2]+
      forecast.pc3$Forecast[k]*fpcaObjFlies$phi[,3]#+
      #forecast.pc4$Forecast[k]*fpcaObjFlies$phi[,4]+
      #forecast.pc5$Forecast[k]*fpcaObjFlies$phi[,5]
    
    pred_error[k] <- mean((true_curve - pred_curve)^2, na.rm = TRUE)
  }
  
  write.csv(data.frame(pred_error = pred_error), file=paste0("CoKriging_output/CoKriging_pred_error_",i,".csv"), row.names = FALSE)
}
stopCluster(cl)
