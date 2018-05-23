function [F] = Throat_eqn(x, gas, P_u, T_u, v_u, h0, alpha_ux)


% x = [P_d; T_d; v_d; h_d];
set(gas, 'T', x(2), 'P', x(1));
R = gasconstant/meanMolecularWeight(gas);
F = [(x(1)*x(3))/x(2) - (P_u*v_u*alpha_ux)/T_u;
    x(1)*(1+(x(3)^2)/(R*x(2))) - P_u*(1+(v_u^2)/(R*T_u))*alpha_ux;
    x(4)+(x(3)^2)/2 - h0;
    x(4)-enthalpy_mass(gas)];

end

