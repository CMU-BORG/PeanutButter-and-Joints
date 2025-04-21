function theta = theta_model(mu_s,mu_k,mu_d,b,t_data,theta0,dtheta0, A, I, r_com, r_ring)
t_end = t_data(end);

dt = 0.01*mean(diff(t_data));
T = 0:dt:t_end;
thetas = zeros(size(T));
theta_i = theta0;
thetas(1) = theta0;
dtheta_i = dtheta0;
ddtheta_i = 0;

for i = 2:length(thetas)
    
    ddtheta = dynamics(theta_i,dtheta_i,ddtheta_i,mu_s,mu_k,mu_d,b,A,I,r_ring,r_com);
    dtheta_i = dtheta_i + ddtheta*dt;
    theta_i = theta_i + dtheta_i*dt;
    
    ddtheta_i = ddtheta;
    
    thetas(i) = theta_i;

end

theta = interp1(T,thetas,t_data);

end
