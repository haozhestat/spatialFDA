lambda_spatial = 10
lambda_temporal = 10
x_limit = c(0, 10)
y_limit = c(0, 10)
t_limit = c(0, 1)
omega1 = 3
omega2 = 2
omega3 = 1
nu1 = 5.5
nu2 = 3.5
nu3 = 1.5
res_sigma = 0.25

mean_fun = function(t) {2*t*sin(2*pi*t)}
phi1_fun = function(t) {sqrt(2)*cos(2*pi*t)}
phi2_fun = function(t) {sqrt(2)*sin(2*pi*t)}
phi3_fun = function(t) {sqrt(2)*cos(4*pi*t)}

num_rep = 200
save_memory = TRUE
max_vector_length = 450000

Delta = 2.5
delta_lower = 0
delta_upper = 1.5
K_t_range <- c(4,3,2,5,6)
K_s_range <- 1:6
K_m_range <- 1:10
K_m = 5
K_t = 4
K_s = 3
deg_t = 3
deg_s = 3


