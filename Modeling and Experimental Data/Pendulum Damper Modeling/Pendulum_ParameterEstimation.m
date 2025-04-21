%% Pendulum Tests
% Last Updated: MB 04/14/2025
clear all;close all;clc

if_save_friction = 0;
if_save_damping = 0;

outfolder = "20250414_Parameters/";
%% Component Properties
% units: g-mm-s
%global A I r_com r_ring

% Bar:
m_bar = 3.9642; % [g] mass
r_com_bar = 59.9; % [mm] radius to center of mass
L_bar = 80; % [mm] length of bar
w_bar = 15; % [mm] width of bar
I_bar = (1/12)*m_bar*(L_bar^2 + w_bar^2) + m_bar*r_com_bar^2; % [g mm^2] moment of inertia

% Finger:
m_finger_calc = 8.13; % [g] mass (from solidworks based on 100% infill)
m_finger = 4.9849; % [g] measured mass of finger
r_com_finger = 17.5; % [mm] radius to center of mass
I_finger = 1159.4*(m_finger / m_finger_calc) + m_finger*r_com_finger^2; % [g mm^2] moment of inertia
r_ring = 11;

% Bolts:
m_B1 = 0; % [g] mass of bolt 1
r_B1 = 24.048; % [mm] radius to center of mass
I_B1 = m_B1*r_B1^2; % [g mm^2] moment of inertia

m_B2 = 0; % [g] mass of bolt 1
r_B2 = 30.9731; % [mm] radius to center of mass
I_B2 = m_B2*r_B2^2; % [g mm^2] moment of inertia

m_B3 = 11.2; % [g] mass of bolt 1
r_B3 = 95.0676; % [mm] radius to center of mass
I_B3 = m_B3*r_B3^2; % [g mm^2] moment of inertia

%% Inertia Properties

g = 9806.65; % [mm s^-2] gravitational acceleration

A = m_bar*r_com_bar + m_finger*r_com_finger + m_B1*r_B1 + m_B2*r_B2 + m_B3*r_B3; % [g mm] total center of mass x mass
r_com = A / (m_bar + m_finger + m_B1 + m_B2 + m_B3); 
I = I_bar + I_finger + I_B1 + I_B2 + I_B3; % [g mm^2] total moment of inertia

%% Simulation parameters

N = 300; % number of bootstrap resamples
K = 5;  % number of samples per boostrap resample

%% Boot Strap Parameter Estimation for Undamped Case

MU_K = zeros(N,1);
MU_D = zeros(N,1);

folder = "No Damper Pendulum Tests\";
file = "No Damper Pendulum Data.xlsx";

% pulling meta data from data file
[~,sheets] = xlsfinfo(folder+file);
num_trials = length(sheets);

if if_save_friction
    tic
    parfor (kk=1:N,10)
        % sample initial condition
        mu_k0 = unifrnd(0,1)*(1e-2);
        mu_d0 = unifrnd(0,1)*(1e-5);

        % randomly select (with replacement) 5 trials
        trials = randi([1,num_trials],K,1);
        fprintf("\nBootstrap Sample %i: Datasets:\n",kk);
    
        % optimize the parameters
        error = @(params) AllTrialError([0,params(1),params(2),0],trials,0,A, I, r_com, r_ring,folder,file);
        params = fminsearch(error,[mu_k0,mu_d0]);%,optimset("Display","iter"));
            
        % save the optimized parameters
        MU_K(kk) = params(1);
        MU_D(kk) = params(2);

    end
    toc
    save(outfolder+"MU_K_vals.mat","MU_K")
    save(outfolder+"MU_D_vals.mat","MU_D")
else
    % load pre-saved parameters
    MU_K = load(outfolder+"MU_K_vals.mat").MU_K;
    MU_D = load(outfolder+"MU_D_vals.mat").MU_D;
end


figure;
subplot(1,2,1)
histogram(MU_K);
title("Coefficient of Kinetic Friction")
subplot(1,2,2)
histogram(MU_D);
title("Coefficient of Viscous Friction")

trials = 1:15; 
fprintf("\nAll Datasets:")
fprintf("\nOptimal Params: mu_k: %f +/- %f, mu_d: %f +/- %f\n",mean(MU_K),std(MU_K),mean(MU_D),std(MU_D))
errorfinal = AllTrialError([0,mean(MU_K),mean(MU_D),0],trials,1,A, I, r_com, r_ring,folder,file);

%% Fitting Damping Parameter for Mid Damping Case

folder_damp_mid = "Peanut Butter Damper Pendulum Tests\";
file_damp_mid = "FullyInjectedPeanutButter.xlsx";

[~,sheets] = xlsfinfo(folder_damp_mid+file_damp_mid);
num_trials = length(sheets);

if if_save_damping

    B_mid = zeros(10*N,1); % more samples because we also need to sample the other param distribution
    mu_k_resamp = zeros(10*N,1); % store the values that were used to fit the data to maintain correlations
    mu_d_resamp = zeros(10*N,1);

    trials = 1:num_trials;
    error0 = AllTrialError([0,mean(MU_K),mean(MU_D),b_0],trials,1,A, I, r_com, r_ring,folder_damp_mid,file_damp_mid);
    sgtitle("Initial Guess")
    
    tic
    parfor (kk=1:length(B_mid),10)
        % randomly sample initial condition
        b0 = unifrnd(0,1)*(1e5);
        trials = randi([1,num_trials],K,1);
        fprintf("\nBootstrap Sample %i: Datasets:\n",kk);
        
        % resample from previous boot strapped distributions
        mu_k = randsample(MU_K,1);
        mu_d = randsample(MU_D,1);
        
        % optimize the parameter
        error = @(param) AllTrialError([0,mu_k,mu_d,param],trials,0,A, I, r_com, r_ring,folder_damp_mid,file_damp_mid);
        param = fminsearch(error,[b_0]);
            
        % save the joint parameter set
        B_mid(kk) = param(1);
        mu_k_resamp(kk) = mu_k;
        mu_d_resamp(kk) = mu_d;
        
    
    end
    toc

    save(outfolder+"B_vals_mid.mat","B_mid")
    save(outfolder+"MU_K_vals_resampled.mat","mu_k_resamp")
    save(outfolder+"MU_D_vals_resampled.mat","mu_d_resamp")

else
    % reload the joint parameter set if not reoptimizing
    B_mid = load(outfolder+"B_vals_mid.mat").B_mid;
    mu_k_resamp = load(outfolder+"MU_K_vals_resampled.mat").mu_k_resamp;
    mu_d_resamp = load(outfolder+"MU_D_vals_resampled.mat").mu_d_resamp;
    
end

trials = 1:num_trials; %randi([1,15],10,1);
fprintf("\nAll Datasets:")
fprintf("\nOptimal Params: b: %f +/- %f\n",mean(B_mid / (1e9)),std(B_mid / (1e9)))
errorfinal = AllTrialError([0,mean(MU_K),mean(MU_D),mean(B_mid)],trials,1,A, I, r_com, r_ring,folder_damp_mid,file_damp_mid);
sgtitle("Optimized Values")

%% Compiling the final bootstrapped joint distribution
MU_K_final = mu_k_resamp;
MU_D_final = mu_d_resamp;
B_final = B_mid;

%% Mean Data
fprintf("\nFrictional Params: \n\tmu_k: %f [%f,%f] x10^-3,\n\tmu_d: %f [%f,%f] x10^-3\n",1000*mean(MU_K_final),1000*quantile(MU_K_final,0.025),1000*quantile(MU_K_final,0.975),...
                                                                                          1000*mean(MU_D_final),1000*quantile(MU_D_final,0.025),1000*quantile(MU_D_final,0.975))
fprintf("\tmu_d median: %.3f x10^-10\n",median(MU_D_final)/(1e-10))

% diving by 10^9 because model is in g-mm-s and we want the parameter in N
% (1000*g * 1000*mm) * m (1000*m) *s/rad
fprintf("\nDamping Param: b = %f [%f,%f] x10^-3 N m s / rad\n ",1000*mean(B_final)/(1e9),1000*quantile(B_final,0.025)/(1e9),1000*quantile(B_final,0.975)/(1e9))

%% Running Monte Carlo sampling on the parameter sets to obtain average model dynamics
theta0 = pi*120/180;    % initial conditions for all simulations
t_data = 0:0.001:6;     % time vector
N = 20;                 % number of resamplings

% preallocating
thetas_undamped = zeros(N,length(t_data));
thetas_damped_mid = zeros(N,length(t_data));
thetas_human = zeros(N,length(t_data));

rng(1000); % for reproducibility of the figure
for i=1:N
    ind = randi(10*N,1); % get random index into bootstrapped joint distribution
    mu_k = MU_K_final(ind); 
    mu_d = MU_D_final(ind); 
    b_mid = B_final(ind); 
    b_human = unifrnd(8.1e6, 14.2e6);

    % run simulation for undamped, PB damped, and human damped
    thetas_undamped(i,:) = theta_model(0,mu_k,mu_d,0,t_data,theta0,0, A, I, r_com, r_ring);
    thetas_damped_mid(i,:) = theta_model(0,mu_k,mu_d,b_mid,t_data,theta0,0,A,I,r_com,r_ring);
    thetas_human(i,:) = theta_model(0,mu_k,mu_d,b_human,t_data,theta0,0,A,I,r_com,r_ring);

end

% Calculate the point-wise means and standard deviations for each group
theta_undamped_mean = mean(thetas_undamped,1);
theta_undamped_std = std(thetas_undamped,[],1);

theta_damped_mid_mean = mean(thetas_damped_mid,1);
theta_damped_mid_std = std(thetas_damped_mid,[],1);

theta_human_mean = mean(thetas_human,1);
theta_human_std = std(thetas_human,[],1);

%% Plotting the average dynamics and confidence intervals
figure('Position',[10,10,1000,400],'Color',[1,1,1]); 
subplot(6,4,[1,2,5,6,9,10,13,14,17,18,21,22]);hold all
plot(t_data,180*theta_undamped_mean/pi,'-','Color',1.2*[0.2,0.2,0.8],'DisplayName','No Damper','LineWidth',1.5)
plot(t_data,180*(theta_undamped_mean+1.96*theta_undamped_std)/pi,'--','Color',0.7*[0.2,0.2,0.8],'HandleVisibility','off','LineWidth',1)
plot(t_data,180*(theta_undamped_mean-1.96*theta_undamped_std)/pi,'--','Color',0.7*[0.2,0.2,0.8],'HandleVisibility','off','LineWidth',1)

plot(t_data,180*theta_damped_mid_mean/pi,'-','Color',1.2*[0.8,0.2,0.2],'DisplayName','With Mid Damper','Linewidth',1.5)
plot(t_data,180*(theta_damped_mid_mean+1.96*theta_damped_mid_std)/pi,'--','Color',0.7*[0.8,0.2,0.2],'HandleVisibility','off','LineWidth',1)
plot(t_data,180*(theta_damped_mid_mean-1.96*theta_damped_mid_std)/pi,'--','Color',0.7*[0.8,0.2,0.2],'HandleVisibility','off','LineWidth',1)

plot(t_data,180*theta_human_mean/pi,'-','Color',1.2*[0.2,0.8,0.2],'DisplayName','Human Finger Damping','Linewidth',1.5)
plot(t_data,180*(theta_human_mean+1.96*theta_human_std)/pi,'--','Color',0.7*[0.2,0.8,0.2],'HandleVisibility','off','LineWidth',1)
plot(t_data,180*(theta_human_mean-1.96*theta_human_std)/pi,'--','Color',0.7*[0.2,0.8,0.2],'HandleVisibility','off','LineWidth',1)

xlabel("Time [s]")
ylabel("Pendulum Angle [deg]")
legend('Location','SE');

set(gca,"FontSize",10)

% Plotting example fits for undamped data
subplot(6,4,[3,7,11]); hold all
folder = "No Damper Pendulum Tests\";
file = "No Damper Pendulum Data.xlsx";
data = readmatrix(folder + file,"Sheet","Trial "+num2str(3));
nan_vals = isnan(data(:,1));
data(nan_vals,:) = [];

t = data(:,1);
theta = data(:,end-2);
theta0 = theta(1);
dtheta0 = (-3*theta(1) + 4*theta(2) - theta(3))/(2*mean(diff(t)));

% running model with mean parameters
theta_data_undamped = theta_model(0,mean(MU_K),mean(MU_D),0,t,theta0,dtheta0,A, I, r_com, r_ring);

plot(t,180*theta/pi,'--','Color',0.6*[0.2,0.2,0.8],'DisplayName','Data','LineWidth',1.5) 
plot(t,180*theta_data_undamped/pi,'-','Color',1.2*[0.2,0.2,0.8],'DisplayName','Model Fit','LineWidth',1.5)
xlabel("Time [s]")
ylabel("Angle [deg]")
legend('Location','SE')
set(gca,"FontSize",10)

% Plotting example fits for damped case
subplot(6,4,[15,19,23]); hold all
data = readmatrix(folder_damp_mid + file_damp_mid,"Sheet","Trial "+num2str(3));
nan_vals = isnan(data(:,1));
data(nan_vals,:) = [];

t = data(:,1);
theta = data(:,end-2);
theta0 = theta(1);
dtheta0 = (-3*theta(1) + 4*theta(2) - theta(3))/(2*mean(diff(t)));

theta_data_damped = theta_model(0,mean(MU_K),mean(MU_D),mean(B_mid),t,theta0,dtheta0,A, I, r_com, r_ring);
plot(t,180*theta/pi,'--','Color',0.6*[0.8,0.2,0.2],'DisplayName','Data','LineWidth',1.5) 
plot(t,180*theta_data_damped/pi,'-','Color',1.2*[0.8,0.2,0.2],'DisplayName','Model Fit','LineWidth',1.5)
xlabel("Time [s]")
ylabel("Angle [deg]")
legend('Location','SE')
set(gca,"FontSize",10)


% Plotting histograms of parameters
subplot(6,4,[4,8]);hold all
yyaxis left
histogram(MU_K,30,'Normalization','probability');
xlabel("\mu_K")
set(gca,"FontSize",10)

subplot(6,4,[12,16]);
histogram(MU_D,30,'Normalization','probability');
xlabel("\mu_D")
set(gca,"FontSize",10)


subplot(6,4,[20,24]);
histogram(B_mid,30,'Normalization','probability');
yticks([])
xlabel("b")
set(gca,"FontSize",10)

%% Calculating gross parameters from the difference tests

folder = "No Damper Pendulum Tests\";
file = "No Damper Pendulum Data.xlsx";

% getting meta data
[~,sheets] = xlsfinfo(folder+file);
num_trials = length(sheets);

% preallocating
T = zeros(num_trials,1);        
N_osc = zeros(num_trials,1);
P = zeros(size(T));

for i=1:num_trials
    % reading in data
    data = readmatrix(folder + file,"Sheet","Trial "+num2str(i));
    nan_vals = isnan(data(:,1));
    data(nan_vals,:) = [];
    
    t = data(:,1);                  % time vector
    theta = data(:,end-2);          % angle vector

    % finding the peaks in the data to count the cycles
    [pks,locs] = findpeaks(180*theta/pi,t,'MinPeakHeight',1.005*180*mean(theta(end-40:end))/pi,'MinPeakDistance',0.4);
    
    N_osc(i) = length(pks);         % number of oscillations
    period = mean(diff(locs));      % average period of oscillations
    P(i) = period;                  
    T(i) = period*(length(pks));    % average length of oscillations
end

fprintf("\nUndamped Case:")
fprintf("\n    Number of oscillations: %f +/- %f",mean(N_osc),std(N_osc));
fprintf("\n    Length of oscillations: %f +/- %f s\n",mean(T),std(T));

% Damped Trials
folder= "Peanut Butter Damper Pendulum Tests\";
file = "FullyInjectedPeanutButter.xlsx";

% getting meta data
[~,sheets] = xlsfinfo(folder+file);
num_trials = length(sheets);

% preallocating
T = zeros(num_trials,1);
N_osc = zeros(num_trials,1);

for i=1:num_trials
    % reading in data
    data = readmatrix(folder + file,"Sheet","Trial "+num2str(i));
    nan_vals = isnan(data(:,1));
    data(nan_vals,:) = [];
    
    t = data(:,1);
    theta = data(:,end-2);hold all
    
    % finding the peaks
    [pks,locs] = findpeaks(180*theta/pi,t,'MinPeakHeight',1.005*180*mean(theta(end-40:end))/pi,'MinPeakDistance',0.4);
    
    N_osc(i) = length(pks);             % number of oscillations
    period = mean(diff(locs));          % average period
    T(i) = 2*locs(1)*(length(pks));     % different calculation because only every one peak 
end

fprintf("\nDamped Case:")
fprintf("\n    Number of oscillations: %f +/- %f",mean(N_osc),std(N_osc));
fprintf("\n    Length of oscillations: %f +/- %f s\n",mean(T),std(T));


