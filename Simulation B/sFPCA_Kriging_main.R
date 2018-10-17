#setwd("~/Dropbox/projects/SPCA/simulation/")
setwd("/vol/data/zhuz/haozhe/simulation/")

source("Simulation_Setting.R")
source("EstMean_Fun.R")
source("Est2DCov_Fun.R")
source("EstFPC_Fun.R")
source("EstSpatialCov_Fun.R")
source("Est2DCovYao_Fun.R")
source("Est2DNugget_Fun.R")
source("EstMeasurementError_Fun.R")

library(splines)

num_samp <- 500
Grid <- seq(t_limit[1], t_limit[2], length.out = num_samp)

#for(i in 71:80){
library(doParallel)
library(foreach)

cl <- makeCluster(5)
registerDoParallel(cl)

foreach(i=20:50, .combine=list,.packages=c("splines")) %dopar% {
tryCatch({
  sim_data <- read.csv(paste0("sim_data/sim_data_", i, ".csv"))
  sim_data_location <- read.csv(paste0("sim_data_location/sim_data_location_", i, ".csv"))
  sim_data_fpcscore <- read.csv(paste0("sim_data_fpcscore/sim_data_fpcscore_", i, ".csv"))
  
  sim_data$value_extmean <-
    sim_data$value -
    EstMean_Fun(sim_data, K_m = K_m, t_limit = t_limit)$EstMean(sim_data$t_coord)
  
  Est3DCov_output <- readRDS(paste0("output/Est3DCov_output_", i,".rds"))
  Est2DCov_output <- Est2DCov_Fun(Est3DCov_output, delta_lower=0.1, delta_upper=delta_upper)
  EstFPC_output <- EstFPC_Fun(Est2DCov_output)
  EstSpatialCov_output <- EstSpatialCov_Fun(Est3DCov_output, EstFPC_output)
  Est2DCovYao_output <- Est2DCovYao_Fun(sim_data, K_t, deg_t, t_limit)
  Est2DNugget_output <- Est2DNugget_Fun(Est2DCovYao_output, Est3DCov_output)
  Est2DNuggetPhi_output <- EstFPC_Fun(Est2DNugget_output)
  EstMeasurementError_output <- EstMeasurementError_Fun(Est2DCovYao_output, sim_data,0.1,0.9)
  measurement_error <- mean(EstMeasurementError_output)
  print(paste0("measurement_error=",measurement_error))
  
  spcoef3D <- Est3DCov_output$spcoef3D
  s_knots <- Est3DCov_output$s_knots
  t_knots <- Est3DCov_output$t_knots
  cov3D <- function(index){
    sbs <- bs(x = index[1], knots = s_knots, intercept=TRUE, degree = deg_s, Boundary.knots = c(0, Delta))
    t1bs <- bs(x = index[2], knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    t2bs <- bs(x = index[3], knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    tbs <- rep(t1bs, each = (K_t+deg_t+1))*rep(t2bs, (K_t+deg_t+1))
    basis3D <- rep(sbs, each = (K_t+deg_t+1)^2)*rep(tbs, (K_s+deg_s+1))
    return(sum(basis3D*spcoef3D))
  }
  
  nugget_cov2D <- function(t){
    t1bs <- bs(x = t[1], knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    t2bs <- bs(x = t[2], knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    c(t1bs%*%Est2DNuggetPhi_output$pd_spcoef2D_mat%*%t(t2bs))
  }
  
  cov_fun <- function(index){
 	 if(index[1]>Delta) 
    	    return(0)
  	 else
    	    return(cov3D(index) + (index[1]==0)*nugget_cov2D(index[2:3]))
  }
    
  pred_error <- rep(NA, nrow(sim_data_location))
  baseline_error <- rep(NA, nrow(sim_data_location))
  pred_xi <- matrix(NA, nrow(sim_data_location), 3) 
  colnames(pred_xi) <- c("xi1", "xi2", "xi3")
  
  sample_location <- sample(1:nrow(sim_data_location), 100)
  for(j in sample_location){
    distance <- as.matrix(dist(sim_data_location[,c("x_coord", "y_coord")]))
    
    index <- which(distance[j,] > 0 & distance[j,] <0.5)
    sub_data <- subset(sim_data, location_id %in% index)
    
    sub_data$dist <- sqrt((sub_data$x_coord - sim_data_location$x_coord[j])^2+
                            (sub_data$y_coord - sim_data_location$y_coord[j])^2)
    sub_data_cov <- matrix(
      apply(
        cbind(
          c(as.matrix(dist(sub_data[,c("x_coord", "y_coord")]))),
          sub_data$t_coord[rep(1:nrow(sub_data), each = nrow(sub_data))],
          sub_data$t_coord[rep(1:nrow(sub_data), nrow(sub_data))]),
        1,
        cov_fun), nrow(sub_data), nrow(sub_data)
    )
    diag(sub_data_cov) <- diag(sub_data_cov) + measurement_error
    while(min(abs(eigen(sub_data_cov)$values)) < 0.1){
      diag(sub_data_cov) <- diag(sub_data_cov) + 0.01
    }
    sub_data_cov_inv <- solve(sub_data_cov)
    
    # shrinkage_ratio <- c(EstSpatialCov_output$SpatialCov_Fun(0,1)/
    #                        (EstSpatialCov_output$SpatialCov_Fun(0,1)+measurement_error),
    #                      EstSpatialCov_output$SpatialCov_Fun(0,2)/
    #                        (EstSpatialCov_output$SpatialCov_Fun(0,2)+measurement_error),
    #                      EstSpatialCov_output$SpatialCov_Fun(0,3)/
    #                        (EstSpatialCov_output$SpatialCov_Fun(0,3)+measurement_error))
    
    pred_xi[j,1] <- t(EstFPC_output$EstPhi_Fun(sub_data$t_coord,1)*
        sapply(sub_data$dist,function(x) EstSpatialCov_output$SpatialCov_Fun(x,1)))%*%
      sub_data_cov_inv%*%sub_data$value_extmean
    
    pred_xi[j,2] <- t(EstFPC_output$EstPhi_Fun(sub_data$t_coord,2)*
                     sapply(sub_data$dist,function(x) EstSpatialCov_output$SpatialCov_Fun(x,2)))%*%
      sub_data_cov_inv%*%sub_data$value_extmean
    
    pred_xi[j,3] <- t(EstFPC_output$EstPhi_Fun(sub_data$t_coord,3)*
                     sapply(sub_data$dist,function(x) EstSpatialCov_output$SpatialCov_Fun(x,3)))%*%
      sub_data_cov_inv%*%sub_data$value_extmean
    
    true_curve <- mean_fun(Grid) +
      phi1_fun(Grid)*sim_data_fpcscore$score1[j]+
      phi2_fun(Grid)*sim_data_fpcscore$score2[j]+
      phi3_fun(Grid)*sim_data_fpcscore$score3[j]+
      nugget_phi1_fun(Grid)*rnorm(1,0,sqrt(nugget_omega1))+
      nugget_phi2_fun(Grid)*rnorm(1,0,sqrt(nugget_omega2))+
      rnorm(length(Grid),0, sqrt(res_sigma))
    
    pred_curve <- 
      EstMean_Fun(sim_data, K_m = K_m, t_limit = t_limit)$EstMean(Grid)+
      pred_xi[j,1]*EstFPC_output$EstPhi_Fun(Grid,1)+
      pred_xi[j,2]*EstFPC_output$EstPhi_Fun(Grid,2)+
      pred_xi[j,3]*EstFPC_output$EstPhi_Fun(Grid,3)
    
    pred_error[j] <- mean((true_curve - pred_curve)^2)
    baseline_error[j] <- mean((true_curve - EstMean_Fun(sim_data, K_m = K_m, t_limit = t_limit)$EstMean(Grid))^2)
    print(c(j,baseline_error[j],pred_error[j]))
    
    write.csv(pred_xi, 
              file=paste0("spatialpace_output_new/spatialpace_output_",i,".csv"),
              row.names = FALSE)
    write.csv(pred_error, 
              file=paste0("spatialpace_output_new/spatialpace_pred_error_",i,".csv"),
              row.names = FALSE)
    write.csv(baseline_error, 
              file=paste0("spatialpace_output_new/spatialpace_baseline_error_",i,".csv"),
              row.names = FALSE)
  }
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
  stopCluster(cl)
