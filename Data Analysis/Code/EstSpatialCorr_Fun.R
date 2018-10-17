EstSpatialCorr_Fun <- function(Est3DCov_output, fpca_output){
  
  Delta <- Est3DCov_output$Delta
  s_knots <- Est3DCov_output$s_knots
  t_knots <- Est3DCov_output$t_knots
  deg_s <- Est3DCov_output$deg_s
  deg_t <- Est3DCov_output$deg_t
  K_s <- Est3DCov_output$K_s
  K_t <- Est3DCov_output$K_t
  spcoef3D <- Est3DCov_output$spcoef3D
  t_limit <- Est3DCov_output$t_limit
  
  basisfd=create.bspline.basis(t_limit,K_t+deg_t+1)
  jmat <- inprod(basisfd, basisfd)
  
  SpatialCorr_Fun <- function(s, k){
    sbs <- bs(x = s, knots = s_knots, intercept=TRUE, degree = deg_s, Boundary.knots = c(0, Delta))
    t1bs <- c(t(fpca_output$EstPhi_beta[,k])%*%jmat)
    t2bs <- c(t(fpca_output$EstPhi_beta[,k])%*%jmat)
    tbs <- rep(t1bs, each = (K_t+deg_t+1))*rep(t2bs, (K_t+deg_t+1))
    basis3D <- rep(sbs, each = (K_t+deg_t+1)^2)*rep(tbs, (K_s+deg_s+1))
    return(sum(basis3D*spcoef3D))
  }
  
  return(list(SpatialCorr_Fun = SpatialCorr_Fun))
}