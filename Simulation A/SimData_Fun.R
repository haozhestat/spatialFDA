SimData_Fun <- function(lambda_spatial = 10,
                        lambda_temporal = 10,
                        x_limit = c(0, 10),
                        y_limit = c(0, 10),
                        t_limit = c(0, 1),
                        omega1 = 4,
                        omega2 = 2,
                        omega3 = 1,
                        nugget_omega1 = 2,
                        nugget_omega2 = 1,
                        nu1 = 5.5,
                        nu2 = 3.5,
                        nu3 = 1.5,
                        scale1 = 1,
                        scale2 = 0.5,
                        scale3 = 0.5,
                        res_sigma = 0.25,
                        mean_fun = function(t) {2*t*sin(2*pi*t)},
                        phi1_fun = function(t) {rep(1,length(t))},
                        phi2_fun = function(t) {sqrt(2)*cos(2*pi*t)},
                        phi3_fun = function(t) {sqrt(2)*sin(2*pi*t)},
                        nugget_phi1_fun = function(t) {sqrt(3)*t},
                        nugget_phi2_fun = function(t) {-2*sqrt(15/7)*t^2+sqrt(15/7)},
                        nugget = TRUE,
                        csv_output = NULL, 
                        plot = FALSE){
  
  library(RandomFields)
  library(spatstat)
  
  RFoptions(seed=NA)
  
  rand_loc <- rpoispp(lambda = lambda_spatial,  win = owin(xrange=x_limit, yrange=y_limit))
  
  x_coord <- rand_loc$x
  y_coord <- rand_loc$y
  n <- rand_loc$n
  distance <- dist(cbind(x_coord,y_coord))
  
  z1_data <- RFsimulate(model = RMtrend(mean=0) + RMmatern(nu = nu1, var=omega1, scale = scale1),
                        x=x_coord, y=y_coord, grid=FALSE)$variable1
  z2_data <- RFsimulate(model = RMtrend(mean=0) + RMmatern(nu = nu2, var=omega2, scale = scale2),
                        x=x_coord, y=y_coord, grid=FALSE)$variable1
  z3_data <- RFsimulate(model = RMtrend(mean=0) + RMmatern(nu = nu3, var=omega3, scale = scale3),
                        x=x_coord, y=y_coord, grid=FALSE)$variable1
  
  if(nugget){
    nugget1_z_data <- rnorm(n,0, sd = sqrt(nugget_omega1))
    nugget2_z_data <- rnorm(n,0, sd = sqrt(nugget_omega2))
  }
  else{
    nugget1_z_data <- rep(0,n)
    nugget2_z_data <- rep(0,n)
  }
  
  sim_data <- list()
  
  for(i in 1:n){
    t_coord <- rpoisppOnLines(lambda = lambda_temporal, 
                              L =psp(0, t_limit[1], 0, t_limit[2], 
                                     window = owin(xrange=c(0,1), yrange=t_limit)))$y
    m <- length(t_coord)
    while(m == 0){
      t_coord <- rpoisppOnLines(lambda = lambda_temporal, 
                                L =psp(0, t_limit[1], 0, t_limit[2], 
                                       window = owin(xrange=c(0,1), yrange=t_limit)))$y
      m <- length(t_coord)
    }
      y_data <- mean_fun(t_coord) +
        phi1_fun(t_coord)*z1_data[i] +
        phi2_fun(t_coord)*z2_data[i] +
        phi3_fun(t_coord)*z3_data[i] +
        nugget_phi1_fun(t_coord)*nugget1_z_data[i] +
        nugget_phi2_fun(t_coord)*nugget2_z_data[i] +
        rnorm(m,0, sd=sqrt(res_sigma))
      
      sim_data[[i]] <- data.frame(location_id = rep(i, m),
                                  x_coord = rep(x_coord[i], m),
                                  y_coord = rep(y_coord[i], m),
                                  t_coord = t_coord,
                                  value = y_data)
  }
  
  sim_data <- do.call("rbind", sim_data)
  
  if(plot){
    library(ggplot2)
    ggplot(data=sim_data) + geom_point(aes(x=x_coord,y=y_coord))+
      theme_bw()+
      theme(panel.grid=element_blank())
    ggplot(data=sim_data) + geom_line(aes(x=t_coord,y=value, color=as.factor(location_id)))+
      theme(legend.position = "none")
    plot(value~t_coord,data=sim_data, pch=".", col="blue", xlab="x", ylab="l")
    lines(seq(0,1,0.01),predict(loess(valeue~t_coord,data=sim_data),seq(0,1,0.01)), col="red", lwd=1)
    lines(seq(0,1,0.01), mean_fun(seq(0,1,0.01)), col="black", lwd=2)
  }
  
  sim_data_fpcscore <- data.frame(location_id = 1:n,
                                  score1 = z1_data,
                                  score2 = z2_data,
                                  score3 = z3_data,
                                  nug_score1 = nugget1_z_data,
                                  nug_score2 = nugget2_z_data)
  
  sim_data_location <- data.frame(location_id = 1:n,
                                  x_coord = x_coord,
                                  y_coord = y_coord)
  
  if(!is.null(csv_output)){
    write.csv(sim_data, file=csv_output[[1]], row.names = FALSE)
    write.csv(sim_data_location, file=csv_output[[2]], row.names = FALSE)
    write.csv(sim_data_fpcscore, file=csv_output[[3]], row.names = FALSE)
  }
  
  return(list(sim_data, sim_data_location, sim_data_fpcscore))
}

