function [F] = mix_eqn(X, gas, term, h_ur, h_us, v_ur, v_us)
%mix_eqn energy equation and momentum equation to be solved

set(gas, 'T', X(2));
gamma = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);

F = [X(1)-enthalpy_mass(gas);
    X(1)-(h_ur+h_us)+0.5*(X(3)^2*(gamma*R*X(2))-(v_ur^2+v_us^2));
    X(3)^2*(gamma*R*X(2))+term*X(3)*sqrt(gamma*R*X(2))+R*X(2)];

end

