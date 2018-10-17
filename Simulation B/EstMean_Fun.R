EstMean_Fun <- function(sim_data = sim_data, K_m = 10, t_limit = c(0, 1)){
  library(splines)
  
  m_knots = quantile(sim_data$t_coord, seq(0,1,length.out = K_m+2)[-c(1,K_m+2)])
  mdesMat <- bs(x = sim_data$t_coord, knots = m_knots, intercept=TRUE, Boundary.knots = t_limit)
  mean_spline_coef <- solve(t(mdesMat)%*%mdesMat)%*%t(mdesMat)%*%sim_data$value
  
  EstMean <- function(t){
    c(bs(x = t, knots = m_knots, intercept=TRUE, Boundary.knots = t_limit)%*%mean_spline_coef)
  }
  
  EstMean <- Vectorize(EstMean)
  
  pred <- EstMean(sim_data$t_coord)
  
  BIC <- log(sum((sim_data$value - pred)^2))*length(pred) + 
    (K_m + 4)*log(length(pred))
  
  AIC <- log(sum((sim_data$value - pred)^2))*length(pred) + (K_m + 4)
  
  return(list(EstMean = EstMean,
              mean_spline_coef = mean_spline_coef,
              BIC = BIC,
              AIC = AIC,
              pred = pred,
              K_m = K_m,
              m_knots = m_knots))
}