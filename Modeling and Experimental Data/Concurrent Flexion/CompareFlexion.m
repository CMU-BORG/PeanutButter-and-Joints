%% Comparison of Concurrent Flexion with and without a Parallel Elastic Element
%{
    This script compares the joint kinematics of our robotic finger with
    and without the parallel elastic element in place. The pair-wise
    correlation coefficients between normalized joint kinematics is also
    calculated.
%}

clear all; close all; clc
warning('off') % suppressing warning messages about xlsx column names not working with MATLAB conventions

%% Reading in data from Tracker spreadsheet
Elastic = readtable("Tracking Data Flexion.xlsx",'Sheet','With Elastic Stiffness');         % data with an elastic element
NonElastic = readtable("Tracking Data Flexion.xlsx",'Sheet','Without Elastic Stiffness');   % data without an elastic element

% pulling out steady state values for each motor position 
idx=315+130.*(0:27).';

q1E=Elastic{idx,"q1"};  % palm to phalange 1 angle
q2E=Elastic{idx,"q2"};  % phalange 1 to phalange 2 angle
q3E=Elastic{idx,"q3"};  % phalange 2 to phalange 3 angle

% same as above but for non-elastic element test
idx2=245+130.*(0:27).';
q1N=NonElastic{idx,"q1"};
q2N=NonElastic{idx,"q2"};
q3N=NonElastic{idx,"q3"};

% position of the motor for each step
motorAngle=[0:27]*10;

%% Calculating the normalized kinematics
T_nom = linspace(0,1,length(q1N));              % normalized motor position

% range-of-motion normalized non-elastic kinematics
J1N = (q1N - min(q1N))/(max(q1N) - min(q1N));
J2N = (q2N - min(q2N))/(max(q2N) - min(q2N));
J3N = (q3N - min(q3N))/(max(q3N) - min(q3N));

% range-of-motion normalized elastic kinematics
J1E = (q1E - min(q1E))/(max(q1E) - min(q1E));
J2E = (q2E - min(q2E))/(max(q2E) - min(q2E));
J3E = (q3E - min(q3E))/(max(q3E) - min(q3E));



%% Plotting the normalized kinematics of the fingers
ms = 6; % marker size

figure('Position',[100,100,600,300],'Color',[1,1,1])
tiledlayout(2,1)

ax1=nexttile;
hold on
plot(T_nom,J1N,'ro','MarkerFaceColor','r','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 1')
plot(T_nom,J2N,'b^','MarkerFaceColor','b','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 2')
plot(T_nom,J3N,'s','MarkerFaceColor',[0 200 0]./255,'MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 3')
xticks([])
ylabel({"Normalized","Joint Angle"})
set(gca,'FontSize',12)


ax2=nexttile;
hold on
plot(T_nom,J1E,'ro','MarkerFaceColor','r','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 1')
plot(T_nom,J2E,'b^','MarkerFaceColor','b','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 2')
plot(T_nom,J3E,'s','MarkerFaceColor',[0 200 0]./255,'MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 3')
xlabel("Normalized Motor Position")
ylabel({"Normalized","Joint Angle"})

set(gca,'FontSize',12)
legend('NumColumns',3,'Location','NW');


%% Calculating the pair-wise correlation coefficients

% concatenating the three joints' kinematics into a single matrix
JN = [J1N,J2N,J3N];
JE = [J1E,J2E,J3E];

R_E = corrcoef(JE);
corrs_E = [R_E(2,1),R_E(3,1),R_E(3,2)];
fprintf("Elastic Joint average correlation: %.3f (range: [%.3f, %.3f])\n", mean(corrs_E),min(corrs_E),max(corrs_E))

R_N = corrcoef(JN);
corrs_N = [R_N(2,1),R_N(3,1),R_N(3,2)];
fprintf("Non-Elastic Joint average correlation: %.3f (range: [%.3f, %.3f])\n", mean(corrs_N),min(corrs_N),max(corrs_N))