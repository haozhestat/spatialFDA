library(data.table)
library(reshape2)
library(splines)
library(splines2)
library(fda)
library(ggplot2)
library(fdapace)
library(geoR)
library (plyr)
library(doParallel)
library(foreach)

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
K_t = 5
K_s = 4
deg_t = 3
deg_s = 3
t_limit = c(min(zillow_data$t_coord), max(zillow_data$t_coord))
Delta = 4
delta = 2.5
delta_lower = 0.1
delta_upper = 3

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

# Est3DCov_output <-
#   Est3DCov_Fun(zillow_data, location, Delta = Delta, K_t = K_t,
#                K_s = K_s, deg_t = deg_t, deg_s = deg_s,
#                t_limit = t_limit, save_memory = TRUE, max_vector_length = 400000)

Est3DCov_output <- readRDS("output/Est3DCov_output_zillow_0927.rds")
Est2DCov_output <- Est2DCov_Fun(Est3DCov_output, delta_lower=delta_lower, delta_upper=delta_upper)
EstFPC_output <- EstFPC_Fun(Est2DCov_output)
EstSpatialCov_output <- EstSpatialCov_Fun(Est3DCov_output, EstFPC_output)
Est2DCovYao_output <- Est2DCovYao_Fun(zillow_data, K_t, deg_t, t_limit)
Est2DNugget_output <- Est2DNugget_Fun(Est2DCovYao_output, Est3DCov_output)
Est2DNuggetPhi_output <- EstFPC_Fun(Est2DNugget_output)
EstMeasurementError_output <- EstMeasurementError_Fun(Est2DCovYao_output, zillow_data,t_limit[1],t_limit[2])
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

pred_error <- rep(NA, nrow(location))
pace_pred_error <- rep(NA, nrow(location))
baseline_error <- rep(NA, nrow(location))
pred_xi <- matrix(NA, nrow(location), 2)
colnames(pred_xi) <- c("xi1", "xi2")
distance <- matrix(
  earth_dist(location[expand.grid(1:nrow(location), 1:nrow(location))[,1],],
             location[expand.grid(1:nrow(location), 1:nrow(location))[,2],]),
  nrow=nrow(location), ncol=nrow(location))
distance[distance>0] <- distance[distance>0] - 0.5


SF_index = unique(zillow_data$location_id[zillow_data$City == "San Francisco"])

for(j in SF_index[51:length(SF_index)]){

  tryCatch({
  index <- which(distance[j,] > 0 & distance[j,] < 1)
  if(length(index) > 0){
    sub_data <- subset(zillow_data, location_id %in% index)
    colnames(sub_data)[7:8] <- c("lon", "lat")
    
    #sub_data$dist <- sqrt((sub_data$x_coord - sim_data_location$x_coord[j])^2+
    #                        (sub_data$y_coord - sim_data_location$y_coord[j])^2)
    sub_data$dist <- earth_dist(sub_data[,c("lon", "lat")],
                                location[rep(j, nrow(sub_data)),]) - 0.5
    sub_data_pairwise_dist <-
      earth_dist(sub_data[expand.grid(1:nrow(sub_data), 1:nrow(sub_data))[,1],c("lon", "lat")],
                 sub_data[expand.grid(1:nrow(sub_data), 1:nrow(sub_data))[,2],c("lon", "lat")])
    
    sub_data_cov <- matrix(
      apply(
        cbind(
          c(sub_data_pairwise_dist),
          sub_data$t_coord[rep(1:nrow(sub_data), each = nrow(sub_data))],
          sub_data$t_coord[rep(1:nrow(sub_data), nrow(sub_data))]),
        1,
        cov_fun), nrow(sub_data), nrow(sub_data)
    )
    diag(sub_data_cov) <- diag(sub_data_cov) + measurement_error
    
    while(min(abs(eigen(sub_data_cov)$values)) < 0.2){
      diag(sub_data_cov) <- diag(sub_data_cov) + 0.01
    }
    sub_data_cov_inv <- solve(sub_data_cov)
    
    pred_xi[j,1] <- t(EstFPC_output$EstPhi_Fun(sub_data$t_coord,1)*
                        sapply(sub_data$dist,function(x) EstSpatialCov_output$SpatialCov_Fun(x,1)))%*%
      sub_data_cov_inv%*%sub_data$value_extmean
    
    pred_xi[j,2] <- t(EstFPC_output$EstPhi_Fun(sub_data$t_coord,2)*
                        sapply(sub_data$dist,function(x) EstSpatialCov_output$SpatialCov_Fun(x,2)))%*%
      sub_data_cov_inv%*%sub_data$value_extmean
    
    true_curve <- zillow_data$value_extmean[zillow_data$location_id == j]
    
    pred_curve <-
      pred_xi[j,1]*EstFPC_output$EstPhi_Fun(zillow_data$t_coord[zillow_data$location_id == j],1)+
      pred_xi[j,2]*EstFPC_output$EstPhi_Fun(zillow_data$t_coord[zillow_data$location_id == j],2)
  
    if(mean((true_curve - pred_curve)^2) > mean((true_curve)^2))
      pred_curve <- 0  
  }
  else{
    true_curve <- zillow_data$value_extmean[zillow_data$location_id == j]
    
    pred_curve <- 0
  }

  
  
  write.csv(data.frame(location_id = j , spatialpace_pred_error = mean((true_curve - pred_curve)^2)),
            file=paste0("spatialpace_output/spatialpace_pred_error_", j, ".csv"),
            row.names = TRUE)
  
  print(j)
  print(mean((true_curve - pred_curve)^2))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
