function [T_d, P_d, x_d, M_d, alpha_d] = Comp_tau(T_u, P_u, x_u, M_u, tau, eta)
%Comp_tau Given the upstream conditions, the temperature ratio and the 
%isentropic efficiency, the downstream conditions are calculated across any compression 

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u)
s = entropy_mass(gas);
h_u = enthalpy_mass(gas);
gamma_u = cp_mass(gas)/cv_mass(gas);
a_u = soundspeed(gas);

T_d = tau*T_u;
setTemperature(gas, T_d);
h_d = enthalpy_mass(gas);

h_di = eta*(h_d-h_u)+h_u;
set(gas, 'P', P_u, 'H', h_di);
T_di = temperature(gas);

P_d = set_ST(moleFractions(gas), s, T_di, T_u, P_u);

pi = P_d/P_u;

set(gas, 'T', T_d, 'P', P_d);
gamma_d = cp_mass(gas)/cv_mass(gas);
a_d = soundspeed(gas);

M_d = sqrt((a_u^2*M_u^2)/a_d^2 - (2/a_d^2)*(h_d-h_u));

alpha_d = (sqrt(tau)/pi)*(M_u/M_d)*sqrt(gamma_u/gamma_d);

x_d = moleFractions(gas);

end

