function theta = theta_model(mu_s,mu_k,mu_d,b,t_data,theta0,dtheta0, A, I, r_com, r_ring)
%{
    Function to solve the dynamical model presented for the pendulum
    experiments and reinterpolating at the corresponding experimental time
    points.

    Inputs:
     - mu_s, mu_k, mu_d, b: model parameters (static friction (not used),
     kinetic friction, viscous friction, viscous damping coefficients)
     - t_data: time points from the experimental data that will be used to
     reinterpolate the model
     - theta0, dtheta0: initial conditions for the model, calculated from
     the experimental data
     - A: collective inertial contributions (sum m_i r_i)
     - I: collective moment of interia
     - r_com: center of mass
     - r_ring: radius of the joint ring
%}

%% Setting up the simulations
% length of the simulation
t_end = t_data(end);

% preallocating vectors and setting up time vector
dt = 0.01*mean(diff(t_data));   % simulate 100x the resolution of the data
T = 0:dt:t_end;
thetas = zeros(size(T));
theta_i = theta0;
thetas(1) = theta0;
dtheta_i = dtheta0;
ddtheta_i = 0;

%% Solving the model
for i = 2:length(thetas)
    % forward Euler integration of the dynamics

    ddtheta = dynamics(theta_i,dtheta_i,ddtheta_i,mu_s,mu_k,mu_d,b,A,I,r_ring,r_com);
    dtheta_i = dtheta_i + ddtheta*dt;
    theta_i = theta_i + dtheta_i*dt;
    
    ddtheta_i = ddtheta;
    
    thetas(i) = theta_i;

end

%% Reinterpolating the model data
theta = interp1(T,thetas,t_data);

end
