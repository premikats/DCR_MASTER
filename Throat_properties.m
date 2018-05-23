function [F, J] = Throat_properties(x, gas, M_u, T_u, P_u, h_u, gamma_u, alpha_ux)
% Throat_properties
% Set of equations to be solved to obtain the downstream conditions
% The upstream conditions are mentioned in the following order:
% u = [M_u, T_u, P_u, h_u]
% The downstream conditions are given in the following order:
% x = [M_d, T_d, P_d, h_d]

set(gas, 'T', x(2), 'P', x(3));
gamma_d = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);

f1 = ((x(3)/P_u)*(x(1)/M_u)*sqrt((gamma_d*T_u)/(gamma_u*x(2))))-alpha_ux; % continuity equation
f2 = ((x(3)/P_u)*((1+gamma_d*x(1)^2)/(1+gamma_u*M_u^2)))-alpha_ux; % Momentum equation
f3 = x(4)+((x(1)^2*gamma_d*R*x(2))/2)-h_u-((M_u^2*gamma_u*R*T_u)/2); % Energy equation
f4 = x(4)-enthalpy_mass(gas);

F = [f1;
    f2;
    f3;
    f4];

if nargout>1
    j11 = (1/M_u)*(x(3)/P_u)*sqrt((gamma_d*T_u)/(gamma_u*x(2)));
    j12 = -0.5*(x(3)/P_u)*(x(1)/M_u)*sqrt(gamma_d*T_u*x(2)/gamma_u);
    j13 = (1/P_u)*(x(1)/M_u)*sqrt((gamma_d*T_u)/(gamma_u*x(2)));
    j14 = 0;
    j21 = (x(3)/P_u)*((2*gamma_d*x(1))/(1+gamma_u*M_u^2));
    j22 = 0;
    j23 = (1/P_u)*((1+gamma_d*x(1)^2)/(1+gamma_u*M_u^2));
    j24 = 0;
    j31 = x(1)*gamma_d*R*x(2);
    j32 = (gamma_d*R*x(1)^2)/2;
    j33 = 0;
    j34 = 1;
    j41 = 0;
    j42 = -cp_mass(gas);
    j43 = 0;
    j44 = 1;
    
    J = [j11 j12 j13 j14;
        j21 j22 j23 j24;
        j31 j32 j33 j34;
        j41 j42 j43 j44];
end
end

