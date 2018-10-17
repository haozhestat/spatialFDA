lambda_spatial = 10
lambda_temporal = 10
x_limit = c(0, 10)
y_limit = c(0, 10)
t_limit = c(0, 1)
omega1 = 3
omega2 = 2
omega3 = 1
nugget_omega1 = 2
nugget_omega2 = 1
nu1 = 5.5
nu2 = 3.5
nu3 = 1.5
scale1 = 1
scale2 = 0.5
scale3 = 0.5
res_sigma = 0.25
mean_fun = function(t) {2*t*sin(2*pi*t)}
phi1_fun = function(t) {sqrt(2)*cos(2*pi*t)}
phi2_fun = function(t) {sqrt(2)*sin(2*pi*t)}
phi3_fun = function(t) {sqrt(2)*cos(4*pi*t)}
nugget_phi1_fun = function(t) {besselJ(2.4048*t,0)*sqrt(t)*sqrt(2)/abs(besselJ(2.4048,1))}
nugget_phi2_fun = function(t) {besselJ(5.5201*t,0)*sqrt(t)*sqrt(2)/abs(besselJ(5.5201,1))}

num_rep = 200
save_memory = TRUE
max_vector_length = 450000

Delta = 2.5
delta_lower = 0
delta_upper = 1.5
K_t_range <- 1:10
K_s_range <- 1:10
K_m_range <- 1:10
K_m = 3
K_t = 3
K_s = 3
deg_t = 3
deg_s = 3

