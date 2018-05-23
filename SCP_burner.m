function [T_d, P_d, x_d, M_d, alpha_dr] = SCP_burner(T_ur, P_ur, x_ur, M_ur, T_us, P_us, x_us, M_us, ra, f)

gas = GRI30;
set(gas, 'T', T_ur, 'P', P_ur, 'X', x_ur)
rho_ur = density(gas);
h_ur = enthalpy_mass(gas);
v_ur = M_ur*soundspeed(gas);

set(gas, 'T', T_us, 'P', P_us, 'X', x_us)
h_us = enthalpy_mass(gas);
v_us = M_us*soundspeed(gas);

v_d = ((1-ra+f)*v_ur+ra*v_us)/(1+f);
h_d = (h_ur+h_us)+0.5*(v_ur^2+v_us^2-v_d^2);

x_u = x_ur+x_us;
P_d = P_ur;

set(gas, 'H', h_d, 'P', P_d, 'X', x_u);
equilibrate(gas, 'HP');

M_d = v_d/soundspeed(gas);

if M_d>1
    T_d = temperature(gas);
    rho_d = density(gas);
    x_d = moleFractions(gas);
    alpha_dr = ((1+f)/(1-ra+f))*(rho_ur/rho_d)*(v_ur/v_d);
else
    T_d = 0;
    rho_d = 0;
    x_d = x_ur;
    alpha_dr = 0;
end