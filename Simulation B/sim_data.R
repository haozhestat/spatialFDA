setwd("~/Dropbox/projects/SPCA/Simulation B/")

source("SimData_Fun.R")
source("Simulation_Setting.R")

for(i in 1:num_rep){
  
  csv_output <- list(paste0("sim_data/sim_data_", i, ".csv"),
                     paste0("sim_data_location/sim_data_location_", i, ".csv"),
                     paste0("sim_data_fpcscore/sim_data_fpcscore_", i, ".csv"))
  
  sim_data <- SimData_Fun(lambda_spatial = lambda_spatial,
                          lambda_temporal = lambda_temporal,
                          x_limit = x_limit,
                          y_limit = y_limit,
                          t_limit = t_limit,
                          omega1 = omega1,
                          omega2 = omega2,
                          omega3 = omega3,
                          nu1 = nu1,
                          nu2 = nu2,
                          nu3 = nu3,
                          res_sigma = res_sigma,
                          mean_fun = mean_fun,
                          phi1_fun = phi1_fun,
                          phi2_fun = phi2_fun,
                          phi3_fun = phi3_fun,
                          nugget = TRUE,
                          csv_output = csv_output, 
                          plot = FALSE)
  print(i)
}

