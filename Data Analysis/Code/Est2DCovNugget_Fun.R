Est2DCovNugget_Fun <- function(data, K_t = 8, deg_t = 3, t_limit = c(0, 1), 
                          max_vector_length = 300000){
  
  library(splines)
  
  t_knots = seq(t_limit[1], t_limit[2], length.out = K_t+2)[-c(1,K_t+2)]
  
  index_grid <- expand.grid(unique(data$location_id),unique(data$location_id))
  index_grid <- subset(index_grid, Var2==Var1)
  index_grid$dist <- earth_dist(location[index_grid$Var1,], location[index_grid$Var2,])
  #index_grid <- subset(index_grid, dist<delta)
  #index_grid$dist[index_grid$dist == 0] <- 0.01
  
  index <- list()
  for(i in 1:nrow(index_grid)){
    index[[i]] <- expand.grid(which(data$location_id == index_grid[i,1]),
                              which(data$location_id == index_grid[i,2]))
    index[[i]]$dist <- index_grid$dist[i]
  }
  
  index_grid <- ldply (index, data.frame)
  rm(index)
  
  pair_data <- data.frame(location_id_1 = data$location_id[index_grid$Var1],
                          location_id_2 = data$location_id[index_grid$Var2],
                          dist = index_grid$dist,
                          t1 = data$t_coord[index_grid$Var1],
                          t2 = data$t_coord[index_grid$Var2],
                          pair_value = data$value_extmean[index_grid$Var1]*
                            data$value_extmean[index_grid$Var2])
  pair_data <- subset(pair_data, t1>t2|t1<t2)
  
  #mean(pair_data$pair_value)
  
  #sdesMat <- bs(x = pair_data$dist, knots = s_knots, intercept=TRUE, degree = deg_s,
  #              Boundary.knots = c(0, delta))
  tdesMat1 <- bs(x = pair_data$t1, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  tdesMat2 <- bs(x = pair_data$t2, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  
  index2d_func <- function(index){
    index_t1 <- ceiling(index/((K_t+deg_t+1)^1))
    #index_t1 <- ceiling((index-(index_s-1)*(K_t+deg_t+1)^2)/(K_t+deg_t+1))
    index_t2 <- index - (index_t1-1)*(K_t+deg_t+1)
    return(c(index_t1, index_t2))
  }
  
  # if(save_memory){
  #   sep_space <- floor(nrow(sdesMat)/max_vector_length)
  #   ind_list <- matrix(NA, sep_space +1 , 2)
  #   for(i in 1:(sep_space+1)){
  #     ind_list[i,] <- c(1+(i-1)*max_vector_length, i*max_vector_length)
  #   }
  #   if(ind_list[sep_space+1,2] == nrow(sdesMat))
  #     ind_list <- ind_list[1:sep_space,]
  #   if(ind_list[sep_space+1,2] > nrow(sdesMat))
  #     ind_list[sep_space+1,2] = nrow(sdesMat)
  #   
  #   XTX <- matrix(0, (K_t+deg_t+1)^2*(K_s+deg_s+1), (K_t+deg_t+1)^2*(K_s+deg_s+1))
  #   X <- matrix(NA, max_vector_length, ncol(sdesMat)*ncol(tdesMat1)*ncol(tdesMat2))
  #   for(loop in 1:nrow(ind_list)){
  #     print(paste0("Calculating the ", loop, " / ", nrow(ind_list), " step..."))
  #     if(loop == nrow(ind_list))
  #       X <- matrix(NA, ind_list[loop,2] - ind_list[loop,1] + 1,
  #                   ncol(sdesMat)*ncol(tdesMat1)*ncol(tdesMat2))
  #     for(i in 1:ncol(X)){
  #       index2d <- index2d_func(i)
  #       X[,i] <- sdesMat[ind_list[loop,1]:ind_list[loop,2],index2d[1]]*
  #         tdesMat1[ind_list[loop,1]:ind_list[loop,2],index2d[2]]*
  #         tdesMat2[ind_list[loop,1]:ind_list[loop,2],index2d[3]]
  #     }
  #     XTX <- XTX + t(X)%*%X
  #   }
  # }
  # 
  # 
  #  if(!save_memory){
  X <- matrix(NA, nrow(tdesMat1), ncol(tdesMat1)*ncol(tdesMat2))
  for(i in 1:(ncol(tdesMat1)*ncol(tdesMat2))){
    index2d <- index2d_func(i)
    X[,i] <- tdesMat1[,index2d[1]]*tdesMat2[,index2d[2]]
  }
  XTX <- t(X)%*%X
  #  }
  
  # XTX <- matrix(0,((K_t+deg_t+1)^2*(K_s+deg_s+1)), ((K_t+deg_t+1)^2*(K_s+deg_s+1)))
  # for(i in 1:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
  #   index2d_i <- index2d_func(i)
  #   Xcol_i <- sdesMat[,index2d_i[1]]*
  #     tdesMat1[,index2d_i[2]]*
  #     tdesMat2[,index2d_i[3]]
  #   for(j in i:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
  #     index2d_j <- index2d_func(j)
  #     if(abs(index2d_i[1]-index2d_j[1])<=deg_s&
  #        abs(index2d_i[2]-index2d_j[2])<=deg_t&
  #        abs(index2d_i[3]-index2d_j[3])<=deg_t){
  # 
  #       Xcol_j <- sdesMat[,index2d_j[1]]*
  #         tdesMat1[,index2d_j[2]]*
  #         tdesMat2[,index2d_j[3]]
  # 
  #       XTX[i,j] <- XTX[j,i] <- sum(Xcol_i*Xcol_j)
  # 
  #     }
  #   }
  #   print(i)
  # }
  
  # XTX <- round(XTX, 7)
  # 
  #  for(i in 1:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
  #    index2d_i <- index2d_func(i)
  #    for(j in i:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
  #      index2d_j <- index2d_func(j)
  #           if(abs(index2d_i[1]-index2d_j[1])>deg_s|
  #              abs(index2d_i[2]-index2d_j[2])>deg_t|
  #              abs(index2d_i[3]-index2d_j[3])>deg_t){
  #             XTX[i,j] <- XTX[j,i] <- 0
  #           }
  #    }
  #  }
  
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
  
  spcoef2D_mat <- matrix(0, K_t+deg_t+1, K_t+deg_t+1)
  for(i in 1:nrow(spcoef2D)){
    spcoef2D_mat[index2d_func(i)[1], index2d_func(i)[2]] <- spcoef2D[i]
  }

  return(list(cov2D = cov2D,
              spcoef2D = spcoef2D,
              spcoef2D_mat = spcoef2D_mat,
              t_knots = t_knots,
              K_t = K_t,
              deg_t = deg_t,
              t_limit = t_limit,
              max_vector_length = max_vector_length))
  #return(list(cov2D = cov2D, spcoef2D = spcoef2D, t_knots = t_knots, 
  #             K_t = K_t, deg_t = deg_t,
  #            t_limit = t_limit, max_vector_length = max_vector_length))
}



# for(i in 1:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
#   index2d_i <- index2d_func(i)
#   Xcol_i <- sdesMat[,index2d_i[1]]*
#     tdesMat1[,index2d_i[2]]*
#     tdesMat2[,index2d_i[3]]
#   for(j in i:((K_t+deg_t+1)^2*(K_s+deg_s+1))){
#     index2d_j <- index2d_func(j)
#     if(abs(index2d_i[1]-index2d_j[1])<=deg_s&
#        abs(index2d_i[2]-index2d_j[2])<=deg_t&
#        abs(index2d_i[3]-index2d_j[3])<=deg_t){
#       
#       Xcol_j <- sdesMat[,index2d_j[1]]*
#         tdesMat1[,index2d_j[2]]*
#         tdesMat2[,index2d_j[3]]
#       
#       XTX[i,j] <- XTX[j,i] <- sum(Xcol_i*Xcol_j)
#       
#     }
#   }
#   print(Sys.time() - a)
#   print(i)
# }
# Sys.time() - b
# cov2D(0.1,0.5,0.4)
# plot(seq(0,1,0.005),2*exp(-seq(0,1,0.005)/0.5)+
#        5*exp(-seq(0,1,0.005)/0.2)*phi2_fun(0.9)*phi2_fun(0.9), type = "l",
#      xlab="dist", ylab="cross-cov", xlim=c(0,1))
# lines(seq(0.01,0.99,0.01),sapply(seq(0.01,0.99,0.01), cov2D, 0.9, 0.9), col="red")



# cov2D(c(0.4,0.5))
# integrate(Vectorize(function(t){cov2D(t,0.4,0.5)}), 0, 1)
# cov2D(c(0.1,0.2))
# x_grid <- seq(0.01,0.99,0.01)
# y_grid <- seq(0.01,0.99,0.01)
# #apply(expand.grid(x_grid, y_grid),1,cov2D)
# fitted_cov2D <- matrix(NA, length(x_grid), length(y_grid))
# true_cov2D <- matrix(NA, length(x_grid), length(y_grid))
# for(i in 1:length(x_grid)){
#   for(j in 1:length(y_grid)){
#     fitted_cov2D[i,j] <- cov2D(c(x_grid[i], y_grid[j]))
#     true_cov2D[i,j] <- 0.5 + 0.5*cos(2*pi*x_grid[i])*cos(2*pi*y_grid[j])
#   }
#   print(i)
# }
# image(x_grid, y_grid, fitted_cov2D)
# plot_ly(x = x_grid, y = y_grid, z = fitted_cov2D) %>% add_surface( opacity = 0.98) %>%add_surface(z = true_cov2D)
# 
# true_cov2D <- function(t){
#   2*(1+cos(pi*t[1])*cos(pi*t[2]))
# }

#setwd("~/Dropbox/Projects/SPCA/sfpcaSim")
#source("EstMean_Fun.R")
#data <- read.csv("data_extmean.csv", header = TRUE)
#  
# pair_fun <- function(index){
#   index1 <- index[1]
#   index2 <- index[2]
#   c(sqrt((data$x_coord[index1] - data$x_coord[index2])^2+
#            (data$y_coord[index1] - data$y_coord[index2])^2),
#     data$t_coord[index1], 
#     data$t_coord[index2],
#     data$value_extmean[index1]*data$value_extmean[index2],
#     as.integer(data$location_id[index1]),
#     as.integer(data$location_id[index2]))
# }
# 
# index_grid <- expand.grid(1:nrow(data), 1:nrow(data))
# index_grid <- subset(index_grid, Var2>Var1)
# rownames(index_grid) <- NULL
# pair_data <- t(apply(index_grid, 1, pair_fun))
# colnames(pair_data) <- c("dist", "t1", "t2", "pair_value", "location_id_1", "location_id_2")
# pair_data <- data.table::data.table(pair_data)
# pair_data <- subset(pair_data, dist<1)
# dim(pair_data)
