EstMeasurementError_Fun <- function(Est2DCovYao_output, 
                                    data,
                                    lower_bound = 0.1,
                                    upper_bound = 0.9){
  
  library(splines)
  library(splines2)
  
  t_knots <- Est2DCovYao_output$t_knots
  deg_t <- Est2DCovYao_output$deg_t
  K_t <- Est2DCovYao_output$K_t
  t_limit <- Est2DCovYao_output$t_limit
  
  t_knots = quantile(data$t_coord, seq(0,1,length.out = K_t+2)[-c(1,K_t+2)])
  
  tdesMat <- bs(x = data$t_coord, knots = t_knots, intercept=TRUE, degree = deg_t,
                 Boundary.knots = t_limit)
  
  totalvariation_spline_coef <- 
    solve(t(tdesMat)%*%tdesMat)%*%t(tdesMat)%*%data$value_extmean^2
  
  TotalVariation <- function(t){
    c(bs(x = t, knots = t_knots, intercept=TRUE, Boundary.knots = t_limit)%*%totalvariation_spline_coef)
  }
  
  return(sapply(seq(lower_bound,upper_bound,0.01), 
                     function(t) TotalVariation(t)-Est2DCovYao_output$cov2D(c(t,t))))
}
