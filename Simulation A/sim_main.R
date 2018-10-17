#setwd("~/Dropbox/projects/spatialFDA/Simulation A")
setwd("/vol/data/zhuz/haozhe/Simulation A/")

source("Simulation_Setting.R")
source("EstMean_Fun.R")
source("Est3DCov_Fun.R")

library(doParallel)
library(foreach)

cl <- makeCluster(10)
registerDoParallel(cl)

foreach(i=1:num_rep) %dopar% {
  sim_data <- read.csv(paste0("sim_data/sim_data_", i, ".csv"), header = TRUE)
  sim_data$value_extmean <-
    sim_data$value -
    EstMean_Fun(sim_data, K_m = K_m, t_limit = t_limit)$EstMean(sim_data$t_coord)
  
  print(paste0("Start Loop: ", i))
  
  Est3DCov_output <- Est3DCov_Fun(data = sim_data,
                                  Delta = Delta,
                                  K_t = K_t,
                                  K_s = K_s,
                                  t_limit = t_limit,
                                  deg_t = deg_t,
                                  deg_s = deg_s,
                                  save_memory = save_memory,
                                  max_vector_length = max_vector_length)
  
  saveRDS(Est3DCov_output, paste0("output/Est3DCov_output_", i,".rds"))
}
stopCluster(cl)
