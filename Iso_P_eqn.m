function [F,J] = Iso_P_eqn(x, gas, P_d, P_u, M_u, T_u, gamma_u, R, h0)

set(gas, 'T', x(2), 'P', P_d);
gamma_d = cp_mass(gas)/cv_mass(gas);

F = [M_u*x(1)*sqrt(gamma_u*gamma_d)*sqrt(x(2)/T_u)-1-gamma_u*M_u^2+(P_d/P_u);
    enthalpy_mass(gas)+(gamma_d*R*x(2)*x(1)^2)/2-h0];

J = [M_u*sqrt(gamma_u*gamma_d)*sqrt(x(2)/T_u) (0.5*M_u*x(1)*sqrt(gamma_u*gamma_d))/sqrt(x(2)*T_u);
    gamma_d*R*x(2)*x(1) cp_mass(gas)+(gamma_d*R*x(1)^2)/2];
end