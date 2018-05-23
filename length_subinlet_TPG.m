function [l_c, l_d] = length_subinlet_TPG(A_u, alpha_d, alpha_t, theta_c, theta_d)

% alpha_d = A_d/A_u;
% alpha_t = A_u/A_t;

% Length of the inlet:
A_d = alpha_d*A_u;
A_t = alpha_t*A_u; %throat area
l_c = (A_u-A_t)/tand(theta_c);
l_d = (A_d-A_t)/tand(theta_d);

end

