Est2DCovYao_Fun <- function(data,
                            K_t = 5,
                            deg_t = 3,
                            t_limit = c(0, 1),
                            knot_selection = FALSE){
  
  library(splines)
  
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
  
  pair_data <- subset(pair_data, t1>t2|t1<t2)
  
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
  
  # cov2D <- function(t1, t2){
  #   t1bs <- bs(x = t1, knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
  #   t2bs <- bs(x = t2, knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
  #   basis2D <- rep(t1bs, each = (K_t+deg_t+1))*rep(t2bs, (K_t+deg_t+1))
  #   return(sum(basis2D*spcoef2D))
  # }
  
  spcoef2D_mat <- matrix(0, K_t+deg_t+1, K_t+deg_t+1)
  for(i in 1:nrow(spcoef2D)){
    spcoef2D_mat[index2d_func(i)[1], index2d_func(i)[2]] <- spcoef2D[i]
  }
  
  cov2D <- function(t){
    t1bs <- bs(x = t[1], knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    t2bs <- bs(x = t[2], knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)
    return(t2bs%*%spcoef2D_mat%*%t(t1bs))
  }
  
  if(knot_selection){
    pred <- apply(pair_data[,c("t1","t2")], 1, cov2D)
    
    BIC <- log(sum((pred-pair_data$pair_value)^2))*length(pred)+
      (K_t+deg_t+1)^2*log(length(pred))
  }
  else{
    pred <- NULL
    BIC <- NULL
  }
  
  
  return(list(#cov2D = cov2D,
              spcoef2D = spcoef2D,
              spcoef2D_mat = spcoef2D_mat,
              t_knots = t_knots,
              K_t = K_t,
              BIC = BIC,
              pred = pred,
              deg_t = deg_t,
              t_limit = t_limit,
              max_vector_length = max_vector_length))
}
