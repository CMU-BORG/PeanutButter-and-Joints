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

mu_s0 = 0.001;
mu_k0 = 0.9*mu_s0;
mu_d0 = 0.0000;
b_0 = 0;

N = 300; % number of bootstrap resamples
K = 5;  % number of samples per boostrap resample

%% Boot Strap Parameter Estimation for Undamped Case

MU_K = zeros(N,1);
MU_D = zeros(N,1);

folder = "No Damper Day 2\";
file = "No Damper Case Day 2.xlsx";

[~,sheets] = xlsfinfo(folder+file);
num_trials = length(sheets);

if if_save_friction
    tic
    parfor (kk=1:N,10)
        mu_k0 = unifrnd(0,1)*(1e-2);
        mu_d0 = unifrnd(0,1)*(1e-5);

        trials = randi([1,num_trials],K,1);
        fprintf("\nBootstrap Sample %i: Datasets:\n",kk);
    
        error = @(params) AllTrialError([0,params(1),params(2),0],trials,0,A, I, r_com, r_ring,folder,file);
        params = fminsearch(error,[mu_k0,mu_d0]);%,optimset("Display","iter"));
            
        MU_K(kk) = params(1);
        MU_D(kk) = params(2);

    end
    toc
    save(outfolder+"MU_K_vals.mat","MU_K")
    save(outfolder+"MU_D_vals.mat","MU_D")
else
    MU_K = load(outfolder+"MU_K_vals.mat").MU_K;
    MU_D = load(outfolder+"MU_D_vals.mat").MU_D;
end

mu_k_mean = mean(MU_K);
mu_k_std = std(MU_K);

mu_d_mean = 0;
mu_d_std = std(MU_D);



figure;
subplot(1,2,1)
histogram(MU_K);
title("Coefficient of Kinetic Friction")
subplot(1,2,2)
histogram(MU_D);
title("Coefficient of Viscous Friction")

trials = 1:15; %randi([1,15],10,1);
fprintf("\nAll Datasets:")
fprintf("\nOptimal Params: mu_k: %f +/- %f, mu_d: %f +/- %f\n",mean(MU_K),std(MU_K),mean(MU_D),std(MU_D))
errorfinal = AllTrialError([0,mean(MU_K),mean(MU_D),0],trials,1,A, I, r_com, r_ring,folder,file);

%% Fitting Damping Parameter for Mid Damping Case

folder_damp_mid = "Concentric Damper\Fully Injected Peanut Butter\";
file_damp_mid = "FullyInjectedPeanutButter.xlsx";

[~,sheets] = xlsfinfo(folder_damp_mid+file_damp_mid);
num_trials = length(sheets);

if if_save_damping

    b_0 = 50000;
    B_mid = zeros(10*N,1); % more samples because we also need to sample the other param distribution
    
    mu_k_resamp = zeros(10*N,1); % store the values that were used to fit the data to maintain correlations
    mu_d_resamp = zeros(10*N,1);

    trials = 1:num_trials;
    error0 = AllTrialError([0,mean(MU_K),mean(MU_D),b_0],trials,1,A, I, r_com, r_ring,folder_damp_mid,file_damp_mid);
    sgtitle("Initial Guess")
    
    tic
    parfor (kk=1:length(B_mid),10)
        b0 = unifrnd(0,1)*(1e5);%
        trials = randi([1,num_trials],K,1);
        fprintf("\nBootstrap Sample %i: Datasets:\n",kk);
        
        % resample from previous boot strapped distributions
        mu_k = randsample(MU_K,1);%normrnd(mu_k_mean,mu_k_std);
        mu_d = randsample(MU_D,1);%abs(normrnd(mu_d_mean,mu_d_std));
    
        error = @(param) AllTrialError([0,mu_k,mu_d,param],trials,0,A, I, r_com, r_ring,folder_damp_mid,file_damp_mid);
        param = fminsearch(error,[b_0]);%,optimset("Display","iter"));
            
        B_mid(kk) = param(1);
        mu_k_resamp(kk) = mu_k;
        mu_d_resamp(kk) = mu_d;
        
    
    end
    toc

    save(outfolder+"B_vals_mid.mat","B_mid")
    save(outfolder+"MU_K_vals_resampled.mat","mu_k_resamp")
    save(outfolder+"MU_D_vals_resampled.mat","mu_d_resamp")

else
    B_mid = load(outfolder+"B_vals_mid.mat").B_mid;
    mu_k_resamp = load(outfolder+"MU_K_vals_resampled.mat").mu_k_resamp;
    mu_d_resamp = load(outfolder+"MU_D_vals_resampled.mat").mu_d_resamp;
    
end

trials = 1:num_trials; %randi([1,15],10,1);
fprintf("\nAll Datasets:")
fprintf("\nOptimal Params: b: %f +/- %f\n",mean(B_mid / (1e9)),std(B_mid / (1e9)))
errorfinal = AllTrialError([0,mean(MU_K),mean(MU_D),mean(B_mid)],trials,1,A, I, r_com, r_ring,folder_damp_mid,file_damp_mid);
sgtitle("Optimized Values")
b_mean_mid = mean(B_mid);
b_std_mid = std(B_mid);

%% Compiling the final bootstrapped joint distribution
MU_K_final = mu_k_resamp;
MU_D_final = mu_d_resamp;
B_final = B_mid;

%% Mean Data
% MU_K(MU_K > mean(MU_K) + 4*std(MU_K)) = []; % removing outliers
% MU_D(MU_D > 4*std(MU_D)) = []; % removing outliers

% fprintf("\nFrictional Params: mu_k: %f +/- %f, mu_d: %f +/- %f\n",mean(MU_K),std(MU_K),mean(MU_D),std(MU_D))
% fprintf("\nDamping Param: b = %f +/- %f\n ",mean(B_mid)/(1e9),std(B_mid)/(1e9))
fprintf("\nFrictional Params: \n\tmu_k: %f [%f,%f] x10^-3,\n\tmu_d: %f [%f,%f] x10^-3\n",1000*mean(MU_K_final),1000*quantile(MU_K_final,0.025),1000*quantile(MU_K_final,0.975),...
                                                                                          1000*mean(MU_D_final),1000*quantile(MU_D_final,0.025),1000*quantile(MU_D_final,0.975))
fprintf("\tmu_d median: %.3f x10^-10\n",median(MU_D_final)/(1e-10))
% diving by 10^9 because model is in g-mm-s and we want the parameter in N
% (1000*g * 1000*mm) * m (1000*m) *s/rad
fprintf("\nDamping Param: b = %f [%f,%f] x10^-3 N m s / rad\n ",1000*mean(B_final)/(1e9),1000*quantile(B_final,0.025)/(1e9),1000*quantile(B_final,0.975)/(1e9))


theta0 = pi*120/180;
t_data = 0:0.001:6;
N = 20;
thetas_undamped = zeros(N,length(t_data));
thetas_damped_mid = zeros(N,length(t_data));
thetas_human = zeros(N,length(t_data));

rng(1000); % for reproducibility of the figure
for i=1:N
    ind = randi(10*N,1); % get random index into bootstrapped joint distribution
    mu_k = MU_K_final(ind); %normrnd(mu_k_mean,mu_k_std);
    mu_d = MU_D_final(ind); %abs(normrnd(mu_d_mean,mu_d_std));
    b_mid = B_final(ind); %normrnd(b_mean_mid,b_std_mid);
    b_human = unifrnd(8.1e6, 14.2e6);
    thetas_undamped(i,:) = theta_model(0,mu_k,mu_d,0,t_data,theta0,0, A, I, r_com, r_ring);
    thetas_damped_mid(i,:) = theta_model(0,mu_k,mu_d,b_mid,t_data,theta0,0,A,I,r_com,r_ring);
    thetas_human(i,:) = theta_model(0,mu_k,mu_d,b_human,t_data,theta0,0,A,I,r_com,r_ring);

end
theta_undamped_mean = mean(thetas_undamped,1);
theta_undamped_std = std(thetas_undamped,[],1);

theta_damped_mid_mean = mean(thetas_damped_mid,1);
theta_damped_mid_std = std(thetas_damped_mid,[],1);

theta_human_mean = mean(thetas_human,1);
theta_human_std = std(thetas_human,[],1);

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

subplot(6,4,[3,7,11]); hold all
folder = "No Damper Day 2\";
file = "No Damper Case Day 2.xlsx";
data = readmatrix(folder + file,"Sheet","Trial "+num2str(3));
nan_vals = isnan(data(:,1));
data(nan_vals,:) = [];

t = data(:,1);
theta = data(:,end-2);
theta0 = theta(1);
dtheta0 = (-3*theta(1) + 4*theta(2) - theta(3))/(2*mean(diff(t)));

theta_data_undamped = theta_model(0,mean(MU_K),mean(MU_D),0,t,theta0,dtheta0,A, I, r_com, r_ring);
plot(t,180*theta/pi,'--','Color',0.6*[0.2,0.2,0.8],'DisplayName','Data','LineWidth',1.5) 
plot(t,180*theta_data_undamped/pi,'-','Color',1.2*[0.2,0.2,0.8],'DisplayName','Model Fit','LineWidth',1.5)
xlabel("Time [s]")
ylabel("Angle [deg]")
legend('Location','SE')
set(gca,"FontSize",10)

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

subplot(6,4,[4,8]);hold all
yyaxis left
histogram(MU_K,30,'Normalization','probability');
xrange = linspace(min(MU_K),max(MU_K),100);
pd = fitdist(MU_K,'normal');
mu_k_dist = pdf(pd,xrange);
yyaxis right
plot(xrange,mu_k_dist,'-r','LineWidth',2)
yticks([])

xlabel("\mu_K")
set(gca,"FontSize",10)

subplot(6,4,[12,16]);
histogram(MU_D,30,'Normalization','probability');
xrange = linspace(min(MU_D),max(MU_D),100);
pd = fitdist(MU_D,'HalfNormal','mu',0);
mu_d_dist = pdf(pd,xrange);
yyaxis right
plot(xrange,mu_d_dist,'-r','LineWidth',2)
yticks([])
xlabel("\mu_D")
set(gca,"FontSize",10)


subplot(6,4,[20,24]);
histogram(B_mid,30,'Normalization','probability');
xrange = linspace(min(B_mid),max(B_mid),100);
pd = fitdist(B_mid,'Burr');
b_dist = pdf(pd,xrange);
yyaxis right
plot(xrange,b_dist,'-r','LineWidth',2)
yticks([])

xlabel("b")
set(gca,"FontSize",10)

%% Values from the difference tests

folder = "No Damper Day 2\";
file = "No Damper Case Day 2.xlsx";

[~,sheets] = xlsfinfo(folder+file);
num_trials = length(sheets);
T = zeros(num_trials,1);
N_osc = zeros(num_trials,1);
P = zeros(size(T));
figure();
for i=1:num_trials
    data = readmatrix(folder + file,"Sheet","Trial "+num2str(i));
    nan_vals = isnan(data(:,1));
    data(nan_vals,:) = [];
    
    t = data(:,1);
    theta = data(:,end-2);hold all
    subplot(ceil(length(trials)/3),3,i); hold all
    plot(t,180*theta/pi,'DisplayName','Data')
    
    v =  (diff(theta))./diff(t);
    [pks,locs] = findpeaks(180*theta/pi,t,'MinPeakHeight',1.005*180*mean(theta(end-40:end))/pi,'MinPeakDistance',0.4);
    plot(locs,pks,'or');
    yyaxis right
    plot(t(1:end-1),v)
    N_osc(i) = length(pks);
    period = mean(diff(locs));
    P(i) = period;
    T(i) = period*(length(pks));
end

fprintf("\nUndamped Case:")
fprintf("\n    Number of oscillations: %f +/- %f",mean(N_osc),std(N_osc));
fprintf("\n    Length of oscillations: %f +/- %f s\n",mean(T),std(T));


folder = "Concentric Damper\Fully Injected Peanut Butter\";
file = "FullyInjectedPeanutButter.xlsx";

[~,sheets] = xlsfinfo(folder+file);
num_trials = length(sheets);
T = zeros(num_trials,1);
N_osc = zeros(num_trials,1);
figure();
for i=1:num_trials
    data = readmatrix(folder + file,"Sheet","Trial "+num2str(i));
    nan_vals = isnan(data(:,1));
    data(nan_vals,:) = [];
    
    t = data(:,1);
    theta = data(:,end-2);hold all
    subplot(ceil(length(trials)/3),3,i); hold all
    plot(t,180*theta/pi,'DisplayName','Data')
    
    v =  (diff(theta))./diff(t);
    [pks,locs] = findpeaks(180*theta/pi,t,'MinPeakHeight',1.005*180*mean(theta(end-40:end))/pi,'MinPeakDistance',0.4);
    plot(locs,pks,'or');
    yyaxis right
    plot(t(1:end-1),v)
    N_osc(i) = length(pks);
    period = mean(diff(locs));
    T(i) = 2*locs(1)*(length(pks));
end

fprintf("\nDamped Case:")
fprintf("\n    Number of oscillations: %f +/- %f",mean(N_osc),std(N_osc));
fprintf("\n    Length of oscillations: %f +/- %f s\n",mean(T),std(T));


function t_end = find_steady_state_time(t,x)
    dt = mean(diff(t));
    t_guess = mean(t);
    ss_not_found = 1;

    x_cutoff = abs(max(x(end-10:end)) - min(x(end-10:end))); %std(x(end-10:end));

    while ss_not_found
        i = find(abs(t-t_guess) < 0.6*dt);
        
        range_vals = abs(max(x(i:end)) - min(x(i:end)));

        if range_vals > x_cutoff
            % steady state not found
            t_guess = mean([t_guess,t(end)]);
        else
            % steady state region found. back track to find the beginning
            while range_vals < x_cutoff
                t_guess = t_guess - 10*dt;
                i = find(abs(t-t_guess) < 0.6*dt);
                range_vals = abs(max(x(i:end)) - min(x(i:end)));
            end
            ss_not_found = 0;
        end

    end
    t_end = t_guess;


end
