function [F] = Comb_eqn(X, gas, v_u, h_u, T_u)
%Comb_eqn momentum and energy equations to be solved

set(gas, 'T', X(2));
R = gasconstant/meanMolecularWeight(gas);

F = [X(1)-enthalpy_mass(gas);
    X(1)-h_u + 0.5*((X(3)*soundspeed(gas))^2-v_u^2);
    (X(3)*soundspeed(gas))^2-(v_u+((R*T_u)/v_u))*X(3)*soundspeed(gas)+R*X(2)];

end

