function [T_d, P_d, x_d, alpha_d] = Comp_M(T_u, P_u, x_u, M_u, M_d, eta)
%Comp_M Gives the exit properties of after adiabatic compression given the
%       exit Mach number

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u)
s = entropy_mass(gas);
h_u = enthalpy_mass(gas);
gamma_u = cp_mass(gas)/cv_mass(gas);
a_u = soundspeed(gas);
v_u = M_u*a_u;

T(1) = T_u;
count = 0;

for i=1:1:10
setTemperature(gas, T(i));
gamma(i) = cp_mass(gas)/cv_mass(gas);
a(i) = soundspeed(gas);
v(i) = M_d*a(i);
h(i) = h_u + (1/2)*(v_u^2-v(i)^2);
set(gas, 'H', h(i), 'P', P_u);
T(i+1) = temperature(gas);
T_f = T(i);
h_f = h(i);
v_f = v(i);
tol = T(i+1)-T(i);
count = count + 1;
if abs(tol)<1e-3
    break;
end
end

T_d = T_f;
tau = T_d/T_u;
h_d = h_f;

h_di = eta*(h_d - h_u) + h_u;
set(gas, 'P', P_u, 'H', h_di);
T_di = temperature(gas);

P_d = set_ST(moleFractions(gas), s, T_di, T_u, P_u);

pi = P_d/P_u;

set(gas, 'T', T_d, 'P', P_d);
gamma_d = cp_mass(gas)/cv_mass(gas);

alpha_d = (sqrt(tau)/pi)*(M_u/M_d)*sqrt(gamma_u/gamma_d);

x_d = moleFractions(gas);

end

