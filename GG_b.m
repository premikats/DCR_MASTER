function [T_d, P_d, x_d, M_d, alpha_d] = GG_b(T_u, P_u, x_u, M_u, fuel, stfu, f, phi, ra)
% GG_b Solves the exit conditions of the Gas Generator Burner

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
h_u = enthalpy_mass(gas);
R_u = gasconstant/meanMolecularWeight(gas);
v_u = M_u*soundspeed(gas);
v_d = ((1-ra)/(1-ra+f))*v_u;
h_d = ((1-ra)/(1-ra+f))*(h_u+(0.5*v_u^2)) - 0.5*v_d^2;

i_f = speciesIndex(gas, fuel);
x_u(i_f) = phi*stfu;

set(gas,'H', h_d,'P', P_u,'MoleFractions', x_u);
equilibrate(gas, 'HP');

P_d = P_u;
T_d = temperature(gas);
x_d = moleFractions(gas);
MW = meanMolecularWeight(gas);
R_d = gasconstant/MW;
gamma_d = cp_mass(gas)/cv_mass(gas);

alpha_d = ((1-ra+f)/(1-ra))*(v_u/v_d)*(R_d/R_u)*(T_d/T_u);

M_d = v_d/soundspeed(gas);


end

