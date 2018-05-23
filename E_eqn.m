function [F, J] = E_eqn(X, gas, M_d, v_u, h_u)

set(gas, 'T', X(2));

gamma = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);

F = [X(1) - enthalpy_mass(gas);
    (X(1)-h_u) + (1/2)*(M_d^2*gamma*R*X(2)-v_u^2)];

J = [1 -cp_mass(gas);
    1 0.5*M_d^2*gamma*R];

end

