%% Parameter Sweep of Effective Damping Coefficient v2
%{
    This script produces heat maps of the achieveable normalized damping
    coefficient for the interlaced fin damper design

    T_joint(dtheta) = -     mu      * b_tilde * dtheta
         [N m]      =   [N s / m^2]    [m^3]    [1/s]

%}

close all; clear all; clc

scale = 1;

D = scale*0.0117; % [m] total diameter of the damper
L_tilde = scale*0.00725*2; % [m] total length of the damper

res = scale*0.0003; % [m] resolution of the printer/minimal wall thickness
tol = scale*0.0004; % [m] tolerance of the printer/minimal channel width

% number of fins for which to generate heatmaps
Ns = [2,3,4,5];
numW = 500;

valid_params = @(w,delta,N) (w>=res).*(delta>=tol);

damp_norm = zeros(numW,numW,length(Ns)); % preallocating matrix for normalized damping coeffs
ws = linspace(0,0.15*D,numW);
deltas = linspace(0,0.15*D,numW);
[WS,DELTAS] = meshgrid(ws,deltas);

for i=1:length(Ns)
    for j=1:numW
        for k=1:numW
            if valid_params(WS(j,k),DELTAS(j,k),Ns(i))
                % valid widths
                damp_norm(j,k,i) = Effective_Damping_Coeff_v2(Ns(i),WS(j,k),DELTAS(j,k),D,L_tilde,0)/(1e-6);
            else
                % invalid widths
                damp_norm(j,k,i) = 0;
            end
        end
    end
end

%% Final Design
final_w = 0.5/1000; % [m] final fin width 0.5mm
final_delta = 0.4/1000; % [m] final gap 0.4mm
final_N_fins = 5;

final_design = Effective_Damping_Coeff_v2(final_N_fins,final_w,final_delta,D,L_tilde,1)/(1e-6);
fprintf("Final normalize damping estimate: %.3f x10^-6 [m^3]\n",final_design)
fprintf("Final damping coefficient estimate: [%.3f, %.3f] x10^-3[N m / s]\n",150000*final_design*(1e-6),250000*final_design*(1e-6))
fprintf("Final design required viscosity range: [%.3f, %.3f] cP\n", 1000*0.0081 ./ (final_design*1e-6) , 1000*0.0142 ./ (final_design*1e-6))


%%
close all
max_damps = zeros(length(Ns),1);

for i=1:length(Ns)
    figure('Position',[10,10,2*300,2*250],'Color',[1,1,1]);

    % subplot(1,length(Ns),i)
    imagesc(ws*1000,deltas*1000,damp_norm(:,:,i));
    max_damps(i) = max(damp_norm(:,:,i),[],'all');

    if Ns(i)==final_N_fins
        % plot the location of the final design
        hold all
        plot(1000*final_w,1000*final_delta,'pw','MarkerSize',15,"MarkerFaceColor",'k')
    end

    set(gca, 'ydir', 'normal','FontName','Arial','FontSize',25)
    colormap turbo
    colorbar
    xlabel('Wall Width [mm]')
    ylabel('Channel Width [mm]')
    % title(['N=',num2str(Ns(i))])
    axis equal

    clim([0,80])
end
% sgtitle('Normalized Damping Coefficient (1e-6 [m^3])')

figure('Position',[10,10,400,250],'Color',[1,1,1])
yyaxis left
plot(Ns,max_damps,'LineWidth',1.5,'HandleVisibility','off')
xlabel('Number of Inner Fins')
ylabel({'Maximum Normalized Damping','Coefficient [1e-6 m^3]'})

C = colororder;

yyaxis right
hold all
mult_fact = 1;
plot(Ns,1000*mult_fact*0.0142 ./ (max_damps*1e-6),'--','Color',C(2,:),'LineWidth',1.5,'DisplayName','Upper Limit')
plot(Ns,1000*mult_fact*mean([0.0142,0.0081]) ./ (max_damps*1e-6),'-','Color',C(2,:),'LineWidth',1.5,'DisplayName','Mean')
plot(Ns,1000*mult_fact*0.0081 ./ (max_damps*1e-6),'--','Color',C(2,:),'LineWidth',1.5,'DisplayName','Lower Limit')
ylabel({'Required Viscosicty for',num2str(mult_fact)+"x Human Damping [cP]"})
%legend()

yline(150e3,'--k','HandleVisibility','off')
yline(200e3,'-k','DisplayName','Peanut Butter Viscocity')
yline(250e3,'--k','HandleVisibility','off')

ylim([0.8,2.6]*(1e5))
xticks(Ns)
%title("Design 2: Maximally Distal Fins")
title("Concentric Cylinder Damper")
