function [F,J] = Throat_isen_newt(x, gas, M_u, T_u, P_u, s_u, gamma_u, h0, alpha_ux)
%Throat_isen Isentropic relations across a diffuser with exit area ratio
%given
% Set of equations to be solved to obtain the downstream conditions
% The downstream conditions are given in the following order:
% x = [M_d, tau_d, pi_d]

set(gas, 'S', s_u, 'P', x(3)*P_u);
gamma_d = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);

f1 = (x(3)*(x(1)/M_u)*sqrt((gamma_d*T_u)/(gamma_u*x(2)*T_u)))-alpha_ux; % continuity equation
f2 = x(2)*T_u-temperature(gas); % Isentropic equation
f3 = enthalpy_mass(gas)+((x(1)^2*gamma_d*R*x(2)*T_u)/2)-h0; % Energy equation
% f4 = x(4)-enthalpy_mass(gas);

F = [f1;
    f2;
    f3];

j11 = x(3)*(1/M_u)*sqrt((gamma_d*T_u)/(gamma_u*x(2)*T_u));
j12 = -0.5*(x(3)*(x(1)/M_u)*sqrt((gamma_d*T_u)/(gamma_u*T_u)))*x(2)^(-3/2);
j13 = (x(1)/M_u)*sqrt((gamma_d*T_u)/(gamma_u*x(2)*T_u));
j21 = 0;
j22 = T_u;
j23 = 0;
j31 = x(1)*gamma_d*R*x(2)*T_u;
j32 = (x(1)^2*gamma_d*R*T_u)/2;
j33 = 0;

J = [j11 j12 j13;
    j21 j22 j23;
    j31 j32 j33];

end