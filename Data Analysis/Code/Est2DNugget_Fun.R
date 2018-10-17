Est2DNugget_Fun <- function(Est2DCovYao_output, Est3DCov_output){

  t_knots <- Est3DCov_output$t_knots
  s_knots <- Est3DCov_output$s_knots
  Delta <- Est3DCov_output$Delta
  K_t <- Est3DCov_output$K_t
  K_s <- Est3DCov_output$K_s
  deg_t <- Est3DCov_output$deg_t
  deg_s <- Est3DCov_output$deg_s
  t_limit <- Est3DCov_output$t_limit

  cov2D <- function(t){
    Est2DCovYao_output$cov2D(t) - Est3DCov_output$cov3D(c(0,t[1],t[2]))
  }

  index3d_func <- function(index){
    index_s <- ceiling(index/((K_t+deg_t+1)^2))
    index_t1 <- ceiling((index-(index_s-1)*(K_t+deg_t+1)^2)/(K_t+deg_t+1))
    index_t2 <- index-(index_s-1)*(K_t+deg_t+1)^2-(index_t1-1)*(K_t+deg_t+1)
    return(c(index_s, index_t1, index_t2))
  }

  spcoef2D_mat <- matrix(0, K_t+deg_t+1, K_t+deg_t+1)
  sdesMat <- bs(x = 0, knots = s_knots, intercept=TRUE, degree = deg_s,
                Boundary.knots = c(0, Delta))
  for(i in 1:nrow(Est3DCov_output$spcoef3D)){
    spcoef2D_mat[index3d_func(i)[2], index3d_func(i)[3]] <-
      spcoef2D_mat[index3d_func(i)[2], index3d_func(i)[3]] +
      sdesMat[index3d_func(i)[1]]*Est3DCov_output$spcoef3D[i]
  }

  spcoef2D_mat <- Est2DCovYao_output$spcoef2D_mat - spcoef2D_mat


  return(list(cov2D = cov2D,
              spcoef2D_mat = spcoef2D_mat,
              t_knots = t_knots,
              K_t = K_t,
              deg_t = deg_t,
              t_limit = t_limit))
}
