Elastic=readtable("Tracking Data Flexion.xlsx",'Sheet','With Elastic Stiffness');
NonElastic=readtable("Tracking Data Flexion.xlsx",'Sheet','Without Elastic Stiffness');

idx=315+130.*(0:27).';
q1E=Elastic{idx,"q1"};
q2E=Elastic{idx,"q2"};
q3E=Elastic{idx,"q3"};
x3E = Elastic{idx,"x_2"};
y3E = Elastic{idx,"y_2"};

idx2=245+130.*(0:27).';
q1N=NonElastic{idx,"q1"};
q2N=NonElastic{idx,"q2"};
q3N=NonElastic{idx,"q3"};
x3N = NonElastic{idx,"x_2"};
y3N = NonElastic{idx,"y_2"};

motorAngle=[0:27]*10;


%% Plot
figure()
tiledlayout(3,1)
ax1=nexttile;
hold on
plot(q1E,'ro')
plot(q1N,'r^')
nexttile
hold on
plot(q2E,'bo');
plot(q2N,'b^');
nexttile
hold on
plot(q3E,'mo');
plot(q3N,'m^');


%% Plot Concurrent flexion

ms = 6;

figure('Position',[100,100,600,300],'Color',[1,1,1])
tiledlayout(2,1)
ax1=nexttile;
hold on
%legend();
T_nom = linspace(0,1,length(q1N));
J1N = (q1N - min(q1N))/(max(q1N) - min(q1N));
J2N = (q2N - min(q2N))/(max(q2N) - min(q2N));
J3N = (q3N - min(q3N))/(max(q3N) - min(q3N));

JN = [J1N,J2N,J3N];

plot(T_nom,J1N,'ro','MarkerFaceColor','r','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 1')
plot(T_nom,J2N,'b^','MarkerFaceColor','b','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 2')
plot(T_nom,J3N,'s','MarkerFaceColor',[0 200 0]./255,'MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 3')
% xlabel("normalized time")
% ylabel("Normalized angle")
% title("Flexion angle without elastic element")
xticks([])
ylabel({"Normalized","Joint Angle"})
set(gca,'FontSize',12)


ax2=nexttile;
hold on
J1E = (q1E - min(q1E))/(max(q1E) - min(q1E));
J2E = (q2E - min(q2E))/(max(q2E) - min(q2E));
J3E = (q3E - min(q3E))/(max(q3E) - min(q3E));

JE = [J1E,J2E,J3E];

plot(T_nom,J1E,'ro','MarkerFaceColor','r','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 1')
plot(T_nom,J2E,'b^','MarkerFaceColor','b','MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 2')
plot(T_nom,J3E,'s','MarkerFaceColor',[0 200 0]./255,'MarkerSize',ms,'MarkerEdgeColor','k','DisplayName','Joint 3')
xlabel("Normalized Motor Position")
ylabel({"Normalized","Joint Angle"})
set(gca,'FontSize',12)
%title("Flexion angle with elastic element")
legend('NumColumns',3,'Location','NW');

R_E = corrcoef(JE);
corrs_E = [R_E(2,1),R_E(3,1),R_E(3,2)];
fprintf("Elastic Joint average correlation: %.3f (range: [%.3f, %.3f])\n", mean(corrs_E),min(corrs_E),max(corrs_E))

R_N = corrcoef(JN);
corrs_N = [R_N(2,1),R_N(3,1),R_N(3,2)];
fprintf("Non-Elastic Joint average correlation: %.3f (range: [%.3f, %.3f])\n", mean(corrs_N),min(corrs_N),max(corrs_N))