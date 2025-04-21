function norm_damp = Effective_Damping_Coeff_v2(N,w,delta,D,L_tilde,to_plot)
%{
    This function is used to calculate the normalized effective damping 
    coefficient for an interlaced cyclindrical damper with N inner fins 
    with inner wall widths w_in and outer wall widths w_out. The returned
    normalized damping coefficient will need to be multiplied by the
    working fluid's viscocity to get the true damping coefficient.

    Fixed design parameters:
     - D: diameter of the total damper = 10 mm
     - dL: width of the channel between inner and outer fins at the end cap
          = 2 mm
     - L_tilde: total length of the damper = 10 mm
    
    Input Parameters:
     - N: number of inner fins [count]
     - w_in: wall width of inner fins [m]
     - w_out: wall width of outer fins [m]
%}

D_eff = D - w; % [m] effective total damper diameter

dL = delta; % [m] end cap channel width
L = L_tilde - w - dL; % [m] overlap between the fins

%{
    Variable names:
     - outer_fins, inner_fins: radius to the mid wall of the different fins
     (this is equivalent to rho in the model).
     - w_cent, ws_i, ws_out: wall widths of the different fin structures
%}

if N==1
    % if there is only going to be one inner fin
    outer_fins = [0.5*D_eff];
    inner_fins = 0;
    w_cent = 2*(outer_fins(1) - delta - 0.5*w);
    ws_in = [w_cent];
    ws_out = w*ones(size(outer_fins));
    
    outer_fins = [-outer_fins(end:-1:1), outer_fins];
    ws_out = [ws_out, ws_out];
else
    outer_fins = linspace(0.5*D_eff- 2*(floor(N/2 + 1.5)-1)*(w + delta), 0.5*D_eff, floor(N/2 + 1.5));
    if mod(N,2)==1
        % if odd number of inner fins
        outer_fins(1) = [];
        inner_fins = (outer_fins(2:end) + outer_fins(1:end-1))/2;
        ws_in = w*ones(size(inner_fins));
        ws_out = w*ones(size(outer_fins));
        
        outer_fins = [-outer_fins(end:-1:1), outer_fins];
        ws_out = [ws_out, ws_out];
        w_cent = 2*(inner_fins(1) - 2*(w+delta) + 0.5*w);
        ws_in = [ws_in, w_cent,ws_in];
        inner_fins = [-inner_fins(end:-1:1), 0, inner_fins];
        
    else
        % if even number of outer fins
        inner_fins = (outer_fins(2:end) + outer_fins(1:end-1))/2;
        outer_fins(1) = [];
        ws_in = w*ones(size(inner_fins));
        ws_out = w*ones(size(outer_fins));
        
        w_cent = 2*(outer_fins(1) - 2*(w+delta) + 0.5*w);
        ws_out = [ws_out, w_cent, ws_out];
        ws_in = [ws_in, ws_in];
        outer_fins = [-outer_fins(end:-1:1),0, outer_fins];
        inner_fins = [-inner_fins(end:-1:1), inner_fins];
    end
end

if w_cent < w
    % check if the center pin is less than the wall thickness. this would
    % be an unprintable geometry and thus should not be assigned a value
    norm_damp = 0;
    return;
end


if to_plot
    PlotFins(outer_fins, inner_fins, D, ws_in, ws_out, L, dL);
end
%% Calculating the effective damping

% contribution from one edge of the fin
t_bar = @(rho,w) ( ((rho + (w/2) + delta).^2) .* ((rho + (w/2)).^2) .* L ) ./ ( delta * (2*(rho + (w/2)) + delta) );
% contribution from the end cap of the fin
t_end = @(rho,w) ( (rho + (w/2)).^4 - (rho - (w/2)).^4 ) ./ ( 8*dL );

% total contribution from the inner most pillar
t_bar_center = @(w) t_bar(0,w) + ((w/2).^4)./(8*dL);
% total contribution from internal fin
t_bar_internal = @(rho,w) t_bar(rho,w) - t_bar(-rho,w) + t_end(rho,w);
% total contribution from outer most fin
t_bar_outer = @(rho,w) -t_bar(-rho,w) + t_end(rho,w);

% initializing the effective damping coefficient for both inner and outer
% fins to make sure they are equal
T_in = 0;
T_out = 0;

if mod(N,2)==0
    % if there are an even number of inner fins (meaning the inner pillar
    % is part of the outer portion)

    % center pin portion
    T_in = T_in + t_bar_center(w_cent); 
    T_out = T_out + ((w_cent/2).^4)./(8*dL);%t_bar_center(w_cent);
    
else
    % if there are an odd number of inner fins (meaning the inner pillar is
    % part of the inner portion)
    
    T_in = T_in + ((w_cent/2).^4)./(8*dL);%t_bar_center(w_in);
    T_out = T_out + t_bar_center(w_cent); % adding reaction force
end

Rs_in = inner_fins(inner_fins > 0); % get only the positive portions (negative portions are connected)

for i = 1:length(Rs_in)
    T_in = T_in + t_bar_internal(Rs_in(i),w); % adding a full fin contribution
    T_out = T_out + t_end(Rs_in(i),w); % adding reaction force from end cap
end

Rs_out = outer_fins(outer_fins > 0);

for i=1:length(Rs_out)-1
    % -1 because last fin is the outer most
    T_out = T_out + t_bar_internal(Rs_out(i),w); % adding a full fin contribution
    T_in = T_in + t_end(Rs_out(i),w); % adding reaction force from end cap
end
T_out = T_out + t_bar_outer(Rs_out(end),w);
T_in = T_in + t_end(Rs_out(end),w);


T_in = 4*pi*T_in;
T_out = 4*pi*T_out;

norm_damp = mean([T_in,T_out]);
end