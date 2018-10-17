EstMeasurementError_Fun <- function(Est3DCov_output, 
                                    data,
                                    lower_bound = 0.2,
                                    upper_bound = 0.8){
  
  library(splines)
  library(splines2)
  
  Delta <- Est3DCov_output$Delta
  s_knots <- Est3DCov_output$s_knots
  t_knots <- Est3DCov_output$t_knots
  deg_s <- Est3DCov_output$deg_s
  deg_t <- Est3DCov_output$deg_t
  K_s <- Est3DCov_output$K_s
  K_t <- Est3DCov_output$K_t
  spcoef3D <- Est3DCov_output$spcoef3D
  t_limit <- Est3DCov_output$t_limit
  
  index_grid <- expand.grid(unique(data$location_id),unique(data$location_id))
  index_grid <- subset(index_grid, Var2==Var1)
  
  index <- list()
  for(i in 1:nrow(index_grid)){
    index[[i]] <- expand.grid(which(data$location_id == index_grid[i,1]),
                              which(data$location_id == index_grid[i,2]))
  }
  
  index_grid <- do.call(rbind, index)
  rm(index)
  
  pair_data <- data.frame(location_id_1 = data$location_id[index_grid$Var1],
                          location_id_2 = data$location_id[index_grid$Var2],
                          t1 = data$t_coord[index_grid$Var1],
                          t2 = data$t_coord[index_grid$Var2],
                          pair_value = data$value_extmean[index_grid$Var1]*
                            data$value_extmean[index_grid$Var2])
  
  t_knots = quantile(sim_data$t_coord, seq(0,1,length.out = K_t+2)[-c(1,K_t+2)])
  
  tdesMat1 <- bs(x = pair_data$t1, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  tdesMat2 <- bs(x = pair_data$t2, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  
  index2d_func <- function(index){
    index_t1 <- ceiling(index/((K_t+deg_t+1)^1))
    index_t2 <- index - (index_t1-1)*(K_t+deg_t+1)
    return(c(index_t1, index_t2))
  }
  
  X <- matrix(NA, nrow(tdesMat1), ncol(tdesMat1)*ncol(tdesMat2))
  for(i in 1:(ncol(tdesMat1)*ncol(tdesMat2))){
    index2d <- index2d_func(i)
    X[,i] <- tdesMat1[,index2d[1]]*tdesMat2[,index2d[2]]
  }
  XTX <- t(X)%*%X
  
  XTX_inv <- solve(XTX)
  
  XTY <- matrix(0, (K_t+deg_t+1)^2, 1)
  for(i in 1:((K_t+deg_t+1)^2)){
    index2d <- index2d_func(i)
    Xcol <- tdesMat1[,index2d[1]]*tdesMat2[,index2d[2]]
    XTY[i,1] <- sum(Xcol*pair_data$pair_value)
  }
  
  spcoef2D <- XTX_inv%*%XTY
  
  cov2D <- function(t1, t2){
    t1bs <- bs(x = t1, knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    t2bs <- bs(x = t2, knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    basis2D <- rep(t1bs, each = (K_t+deg_t+1))*rep(t2bs, (K_t+deg_t+1))
    return(sum(basis2D*spcoef2D))
  }
  
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
  
  return(sapply(seq(lower_bound,upper_bound,0.01), 
                     function(t) max(cov2D(t,t) - cov3D(c(0,t,t)),0)))
}

