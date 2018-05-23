function [T_d, P_d, x_d, M_d, alpha_d, alpha_dc, R_d] = Iso(T_u, P_u, x_u, M_u, k)
%Isolator Gives the exit properties of after the isolator given the factor
%           of decrease of the exit Mach number

M_d = k*M_u;

gas = GRI30;
set(gas, 'P', P_u, 'T', T_u, 'X', x_u);
h_u = enthalpy_mass(gas);
gamma_u = cp_mass(gas)/cv_mass(gas);
a_u = soundspeed(gas);
v_u = M_u*a_u;

X0 = [h_u;
      T_u];
% sol = fsolve(@(x)E_eqn(x, gas, M_d, v_u, h_u), X0);
sol = newton(@E_eqn, X0, gas, M_d, v_u, h_u, 1e-8);

% h_d = sol(1);
T_d = sol(2);

setTemperature(gas, T_d);
gamma_d = cp_mass(gas)/cv_mass(gas);

tau = T_d/T_u;

pi = (1 + gamma_u*M_u^2 - M_u*M_d*sqrt(gamma_u*gamma_d)*sqrt(tau));
P_d = pi*P_u;

alpha_dc = (1/pi)*(M_u/M_d)*sqrt((gamma_u*tau)/(gamma_d));
alpha_d = 1;

x_d = moleFractions(gas);
R_d = gasconstant/meanMolecularWeight(gas);


end

