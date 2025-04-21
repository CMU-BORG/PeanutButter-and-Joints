function ddtheta = dynamics(theta,dtheta,ddtheta_i,mu_s,mu_k,mu_d,b,A,I,r,r_com)

g = 9806.65; % [ mm / s^2 ]
N = A.*dtheta.^2 + A.*g.*cos(pi - theta); % A = sum ( m_j r_j )

tau_fric = abs(N).*r.*sign(dtheta).*(mu_k + mu_d.*r.*dtheta);

ddtheta = (1./I).*(A.*g.*sin(theta) - tau_fric - b.*dtheta);

end
