EstFPC_Fun <- function(Est2DCov_output){
  library(fda)

  spcoef2D_mat = as.matrix(Est2DCov_output$spcoef2D_mat)
  t_limit <- Est2DCov_output$t_limit
  t_knots <- Est2DCov_output$t_knots
  K_t <- Est2DCov_output$K_t
  deg_t <- Est2DCov_output$deg_t
  
  basisfd=create.bspline.basis(t_limit, nbasis = K_t+deg_t+1, breaks=c(t_limit[1], t_knots, t_limit[2]))
  jmat <- inprod(basisfd, basisfd)
  Smat=chol(jmat)
  eigen_output <- eigen(Smat%*%Est2DCov_output$spcoef2D_mat%*%t(Smat))
  
  pd_spcoef2D_mat <-
    solve(Smat)%*%(eigen_output$vectors)%*%
    diag(sapply(eigen_output$values,max,0))%*%t(eigen_output$vectors)%*%solve(t(Smat))
  
  EstPhi_beta <- solve(Smat)%*%eigen_output$vectors
  
  EstPhi_Fun <- function(t, k){
    bs(x = t, knots = t_knots, intercept=TRUE, degree = deg_t, Boundary.knots = t_limit)%*%EstPhi_beta[,k]
  }
  
  return(list(EstOmega = eigen_output$values,
              EstPhi_Fun = EstPhi_Fun,
              EstPhi_beta = EstPhi_beta,
              spcoef2D_mat = spcoef2D_mat,
              pd_spcoef2D_mat = pd_spcoef2D_mat))
}