function PlotFins(outer_fins, inner_fins, D, w_in, w_out, L, dL)
%{
    Function to plot the cross-sectional geometry of the concentric ring
    dampers, showing the interdigitated structures.
%}

cap_y = [];
cap_x = [];

%% Get all of the outer finger points
for i=1:length(outer_fins)
    cap_y = [cap_y, outer_fins(i) + [-w_out(i)/2,-w_out(i)/2,w_out(i)/2,w_out(i)/2]];
    cap_x = [cap_x, [0,L,L,0]];
end

% adding end cap
cap_y = [cap_y, [0.5*D,-0.5*D,-0.5*D]];
cap_x = w_out(1) + [cap_x, [-w_out(1),-w_out(1),0]];

%% Get all of the inner finger points
fin_y = [];
fin_x = [];
for i=1:length(inner_fins)
    fin_y = [fin_y, inner_fins(i) + [-w_in(i)/2,-w_in(i)/2,w_in(i)/2,w_in(i)/2]];
    fin_x = [fin_x, [dL+L,dL,dL,dL+L]];
end

% adding end cap
fin_y = [fin_y, [0.5*D,0.5*D,-0.5*D,-0.5*D,fin_y(1)]];
fin_x = w_out(1) + [fin_x, [dL+L,dL+L+w_in(1),dL+L+w_in(1),dL+L,dL+L]];

%% Plotting the two interdigitating structures
figure; hold all
plot(cap_x,cap_y);
plot(fin_x,fin_y);
ylim([-1.1*(D/2),1.1*(D/2)])
axis equal

end