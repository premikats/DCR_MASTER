function Length = length_nozzle_TPG(T_u, P_u, x_u, M_u, A_u, T_d, P_d, x_d, M_d, alpha_d)
% Length of Nozzle:

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
gamma_u = cp_mass(gas)/cv_mass(gas);
set(gas, 'T', T_d, 'P', P_d, 'X', x_d);
gamma_d = cp_mass(gas)/cv_mass(gas);
gamma = average(gamma_u, gamma_d);

nu_d = Prandlt_Meyer_CPG(M_d, gamma);
nu_u = Prandlt_Meyer_CPG(M_u, gamma);
theta = (nu_d-nu_u)/2;
% theta = 2;
Length = A_u*(alpha_d-1)/tand(theta);

end

