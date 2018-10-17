
#library(RandomFields)
library(data.table)
library(reshape2)
library(splines)
library(splines2)
library(fda)
#library(latex2exp)
#library(ggplot2)
library (plyr)

#setwd("~/Dropbox/projects/SPCA/Zillow/")
setwd("/vol/data/zhuz/haozhe/Zillow/")
source("earth_dist.R")
source("EstMean_Fun.R")
source("Est3DCov_Fun.R")
source("Est2DCov_Fun.R")
source("EstNugget_Fun.R")
source("fpca_Fun.R")
source("Est2DCovNugget_Fun.R")
source("EstSpatialCorr_Fun.R")

zillow_data <- readRDS("zillow_data.rds")
location <- readRDS("location.rds")
location <- location[,c("lon", "lat")]
zillow_data <- na.omit(zillow_data)

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

#zillow_data$value_extmean <- zillow_data$value -
#  EstMean_Fun(data = zillow_data, K_m = K_m, t_limit = t_limit)$EstMean(zillow_data$t_coord)

Est3DCov_output <- Est3DCov_Fun(zillow_data, location, Delta = Delta, K_t = K_t, K_s = K_s, deg_t = deg_t, deg_s = deg_s, 
                                t_limit = t_limit, save_memory = TRUE, 
                                max_vector_length = 400000)

saveRDS(Est3DCov_output, "output/Est3DCov_output_zillow_0927.rds")

Est2DCovNugget_output <- Est2DCovNugget_Fun(zillow_data, K_t = K_t, deg_t = deg_t, t_limit = t_limit, 
                                            max_vector_length = 300000)

saveRDS(Est2DCovNugget_output, "output/Est2DCovNugget_output_0927.rds")


