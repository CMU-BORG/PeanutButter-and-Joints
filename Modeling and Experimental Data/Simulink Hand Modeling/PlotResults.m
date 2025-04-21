%% Plot kinematics, torque and contact info for one finger

figure("Position",[100,100,500,400],"Color","w")
subplot(2,1,1); hold all

plot(out.Finger1Data.Segment1_Kinematics.Angle.Time, ...
    rad2deg(out.Finger1Data.Segment1_Kinematics.Angle.Data),'r-','DisplayName',"\theta_0^1","LineWidth",2);
plot(out.Finger1Data.Segment2_Kinematics.Angle.Time, ...
    rad2deg(out.Finger1Data.Segment2_Kinematics.Angle.Data)-rad2deg(out.Finger1Data.Segment1_Kinematics.Angle.Data),'b--','DisplayName',"\theta_1^2","LineWidth",2);
plot(out.Finger1Data.Segment2_Kinematics.Angle.Time, ...
    rad2deg(out.Finger1Data.Segment3_Kinematics.Angle.Data)-rad2deg(out.Finger1Data.Segment2_Kinematics.Angle.Data),'m.-','DisplayName',"\theta_2^3","LineWidth",2);
legend("NumColumns",3,"Location","nw");
xlabel('Time [s]');
ylabel('Angle [deg]');
set(gca,"FontName","Arial","FontSize",15)

subplot(2,1,2); hold all
plot(out.ballInfo.ball_y,'b-','DisplayName','Ball y',"LineWidth",2);
xlabel('Time [s]');
ylabel('Ball Position [m]');
set(gca,"FontName","Arial","FontSize",15)
