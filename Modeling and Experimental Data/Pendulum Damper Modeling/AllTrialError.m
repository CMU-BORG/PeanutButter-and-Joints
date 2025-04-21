function error = AllTrialError(params, trials, if_plot,A, I, r_com, r_ring,folder,file)

error = 0;

if if_plot
    figure('Position',[10,10,1200,600]);
end

for i=1:length(trials)
    data = readmatrix(folder + file,"Sheet","Trial "+num2str(trials(i)));
    nan_vals = isnan(data(:,1));
    data(nan_vals,:) = [];
    
    t = data(:,1);
    theta = data(:,end-2);
    theta0 = theta(1);
    dtheta0 = (-3*theta(1) + 4*theta(2) - theta(3))/(2*mean(diff(t)));
    model_data = theta_model(params(1),params(2),params(3),params(4),t,theta0,dtheta0,A, I, r_com, r_ring);

    if if_plot
        subplot(ceil(length(trials)/3),3,i); hold all
        plot(t,180*theta/pi,'DisplayName','Data')
        plot(t,180*model_data/pi,'DisplayName','Initial Guess Model')
        legend()
        title("Trial" + num2str(trials(i)))
    end
    a = 180/pi;
    error = error + sum((a*theta - a*model_data).^2);

end

negative_penalty = 0;
for i=1:length(params)
    if params(i)<0
        negative_penalty = negative_penalty + 1e9;
    end
end

error = 1*error + negative_penalty;% + 1000*(params(1)<params(2)) ;


end
