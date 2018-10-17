Est3DCov_Fun <- function(data = sim_data,
                         Delta = 1,
                         K_t = 5,
                         K_s = 3,
                         t_limit = c(0, 1),
                         deg_t = 3,
                         deg_s = 3,
                         save_memory = TRUE,
                         max_vector_length = 300000){
  
  library(splines)
  
  index_grid <- expand.grid(1:nrow(data), 1:nrow(data))
  index_grid <- subset(index_grid, Var2!=Var1)
  index_grid$dist <- sqrt((data$x_coord[index_grid$Var1]-
                             data$x_coord[index_grid$Var2])^2+
                            (data$y_coord[index_grid$Var1]-
                               data$y_coord[index_grid$Var2])^2)
  index_grid <- subset(index_grid, dist < Delta)
  
  pair_data <- data.frame(location_id_1 = data$location_id[index_grid$Var1],
                          location_id_2 = data$location_id[index_grid$Var2],
                          dist = index_grid$dist,
                          t1 = data$t_coord[index_grid$Var1],
                          t2 = data$t_coord[index_grid$Var2],
                          pair_value = data$value_extmean[index_grid$Var1]*
                            data$value_extmean[index_grid$Var2])
  rm(index_grid)
  
  s_knots = quantile(pair_data$dist, seq(0,1,length.out = K_s+2)[-c(1,K_s+2)])
  t_knots = quantile(data$t_coord, seq(0,1,length.out = K_t+2)[-c(1,K_t+2)])
  
  sdesMat <- bs(x = pair_data$dist, knots = s_knots, intercept=TRUE, degree = deg_s,
                Boundary.knots = c(0, Delta))
  tdesMat1 <- bs(x = pair_data$t1, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  tdesMat2 <- bs(x = pair_data$t2, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  
  index3d_func <- function(index){
    index_s <- ceiling(index/((K_t+deg_t+1)^2))
    index_t1 <- ceiling((index-(index_s-1)*(K_t+deg_t+1)^2)/(K_t+deg_t+1))
    index_t2 <- index-(index_s-1)*(K_t+deg_t+1)^2-(index_t1-1)*(K_t+deg_t+1)
    return(c(index_s, index_t1, index_t2))
  }
  
  if(save_memory){
    sep_space <- floor(nrow(sdesMat)/max_vector_length)
    ind_list <- matrix(NA, sep_space +1 , 2)
    for(i in 1:(sep_space+1)){
      ind_list[i,] <- c(1+(i-1)*max_vector_length, i*max_vector_length)
    }
    if(ind_list[sep_space+1,2] == nrow(sdesMat))
      ind_list <- ind_list[1:sep_space,]
    if(ind_list[sep_space+1,2] > nrow(sdesMat))
      ind_list[sep_space+1,2] = nrow(sdesMat)
    
    XTX <- matrix(0, (K_t+deg_t+1)^2*(K_s+deg_s+1), (K_t+deg_t+1)^2*(K_s+deg_s+1))
    X <- matrix(NA, max_vector_length, ncol(sdesMat)*ncol(tdesMat1)*ncol(tdesMat2))
    for(loop in 1:nrow(ind_list)){
      print(paste0("Calculating the ", loop, " / ", nrow(ind_list), " step..."))
      if(loop == nrow(ind_list))
        X <- matrix(NA, ind_list[loop,2] - ind_list[loop,1] + 1,
                    ncol(sdesMat)*ncol(tdesMat1)*ncol(tdesMat2))
      for(i in 1:ncol(X)){
        index3d <- index3d_func(i)
        X[,i] <- sdesMat[ind_list[loop,1]:ind_list[loop,2],index3d[1]]*
          tdesMat1[ind_list[loop,1]:ind_list[loop,2],index3d[2]]*
          tdesMat2[ind_list[loop,1]:ind_list[loop,2],index3d[3]]
      }
      XTX <- XTX + t(X)%*%X
    }
  }
  
  if(!save_memory){
    X <- matrix(NA, nrow(sdesMat), ncol(sdesMat)*ncol(tdesMat1)*ncol(tdesMat2))
    for(i in 1:(ncol(sdesMat)*ncol(tdesMat1)*ncol(tdesMat2))){
      index3d <- index3d_func(i)
      X[,i] <- sdesMat[,index3d[1]]*tdesMat1[,index3d[2]]*tdesMat2[,index3d[3]]
    }
    XTX <- t(X)%*%X
  }
  
  XTX_inv <- solve(XTX)
  
  XTY <- matrix(0, (K_t+deg_t+1)^2*(K_s+deg_s+1), 1)
  for(i in 1:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
    index3d <- index3d_func(i)
    Xcol <- sdesMat[,index3d[1]]*tdesMat1[,index3d[2]]*tdesMat2[,index3d[3]]
    XTY[i,1] <- sum(Xcol*pair_data$pair_value)
  }
  
  spcoef3D <- XTX_inv%*%XTY
  
  return(list(#cov3D = cov3D,
              spcoef3D = spcoef3D,
              t_knots = t_knots,
              s_knots = s_knots,
              index3d_func = index3d_func,
              Delta = Delta,
              K_t = K_t,
              K_s = K_s,
              deg_t = deg_t,
              deg_s = deg_s,
              t_limit = t_limit,
              save_memory = save_memory,
              max_vector_length = max_vector_length))
}
