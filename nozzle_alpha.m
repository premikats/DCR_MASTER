function alpha_d = nozzle_alpha(T_u, P_u, x_u, M_u, P_d, eta)
%nozzle Gives the exit properties across the nozzle

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
h_u = enthalpy_mass(gas);
s = entropy_mass(gas);
v_u = M_u*soundspeed(gas);

set(gas, 'S', s, 'P', P_d);
h_di = enthalpy_mass(gas);

h_d = h_u + eta*(h_di-h_u);
v_d = sqrt(v_u^2 + 2*(h_u-h_d));
set(gas, 'H', h_d, 'P', P_d);
T_d = temperature(gas);

alpha_d = (v_u/v_d)*(P_u/P_d)*(T_d/T_u);


end

