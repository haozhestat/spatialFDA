Est2DCov_Fun <- function(Est3DCov_output,
                         delta_lower = 0,
                         delta_upper = NULL){
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

  if(is.null(delta_upper))
    delta_upper <- Delta

  index3d_func <- function(index){
    index_s <- ceiling(index/((K_t+deg_t+1)^2))
    index_t1 <- ceiling((index-(index_s-1)*(K_t+deg_t+1)^2)/(K_t+deg_t+1))
    index_t2 <- index-(index_s-1)*(K_t+deg_t+1)^2-(index_t1-1)*(K_t+deg_t+1)
    return(c(index_s, index_t1, index_t2))
  }

  spcoef2D_mat <- matrix(0, K_t+deg_t+1, K_t+deg_t+1)
  sbs <- ibs(x = delta_upper, knots = s_knots, intercept=TRUE, degree = deg_s,Boundary.knots =c(0, Delta)) -
    ibs(x = delta_lower, knots = s_knots, intercept=TRUE, degree = deg_s,Boundary.knots =c(0, Delta))

  for(i in 1:nrow(spcoef3D)){
    spcoef2D_mat[index3d_func(i)[2], index3d_func(i)[3]] <-
      spcoef2D_mat[index3d_func(i)[2], index3d_func(i)[3]] +
      sbs[index3d_func(i)[1]]*spcoef3D[i]
  }

  # cov2D <- function(t){
  #   t1 = t[,1]
  #   t2 = t[,2]
  #   t1bs <- bs(x = t1, knots = t_knots, intercept=TRUE, degree = deg_t,
  #              Boundary.knots = t_limit)
  #   t2bs <- bs(x = t2, knots = t_knots, intercept=TRUE, degree = deg_t,
  #              Boundary.knots = t_limit)
  #   cov2D_value <- NULL
  #   for(i in 1:nrow(t)){
  #     cov2D_value[i] <- (t1bs[i,]%*%spcoef2D_mat)%*%t(t2bs)[,i]
  #   }
  #   return(cov2D_value)
  # }

  return(list(#cov2D = cov2D,
              spcoef2D_mat =spcoef2D_mat,
              spcoef3D = spcoef3D,
              t_limit = t_limit,
              t_knots = t_knots,
              s_knots = s_knots,
              index3d_func = index3d_func,
              Delta = Delta,
              K_t = K_t,
              K_s = K_s,
              deg_t = deg_t,
              deg_s = deg_s))
}
