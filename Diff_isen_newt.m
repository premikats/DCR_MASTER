function [F,J] = Diff_isen_newt(x, gas, M_u, T_u, P_u, s_u, gamma_u, h0, M_d)
%Diff_isen Isentropic relations across a diffuser with exit Mach number
%given
% Set of equations to be solved to obtain the downstream conditions
% The downstream conditions are given in the following order:
% x = [tau_d, pi_d, alpha_d]

set(gas, 'S', s_u, 'P', x(2)*P_u);
gamma_d = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);

f1 = (x(2)*x(3)*sqrt((gamma_d*T_u)/(gamma_u*x(1)*T_u)))-(M_u/M_d); % continuity equation
f2 = x(1)*T_u-temperature(gas); % Isentropic equation
f3 = enthalpy_mass(gas)+((M_d^2*gamma_d*R*x(1)*T_u)/2)-h0; % Energy equation
% f4 = x(4)-enthalpy_mass(gas);

F = [f1;
    f2;
    f3];

j11 = -0.5*(x(2)*x(3)*sqrt((gamma_d*T_u)/(gamma_u*T_u)))*(x(1)^(-3/2));
j12 = x(3)*sqrt((gamma_d*T_u)/(gamma_u*x(1)*T_u));
j13 = x(2)*sqrt((gamma_d*T_u)/(gamma_u*x(1)*T_u));
j21 = T_u;
j22 = 0;
j23 = 0;
j31 = (M_d^2*gamma_d*R*T_u)/2;
j32 = cp_mass(gas);
j33 = 0;

J = [j11 j12 j13;
    j21 j22 j23;
    j31 j32 j33];

end
