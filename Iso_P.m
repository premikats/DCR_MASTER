function [T_d, x_d, M_d, alpha_d, alpha_dc, R_d] = Iso_P(T_u, P_u, x_u, M_u, P_d)

pi = P_d/P_u;

if pi>1
    gas = GRI30;
    set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
    gamma_u = cp_mass(gas)/cv_mass(gas);
    R = gasconstant/meanMolecularWeight(gas);
    h_u = enthalpy_mass(gas);
    h0 = h_u+(gamma_u*R*T_u*M_u^2)/2;
    
    x0 = [M_u;
        T_u];
    sol = newton(@Iso_P_eqn, x0, gas, P_d, P_u, M_u, T_u, gamma_u, R, h0, 1e-8);
    
    M_d = sol(1);
    T_d = sol(2);
    
    if M_d>1
        set(gas, 'T', T_d, 'P', P_d);
        gamma_d = cp_mass(gas)/cv_mass(gas);
        
        tau = T_d/T_u;
        alpha_dc = (1/pi)*(M_u/M_d)*sqrt((gamma_u*tau)/(gamma_d));
        alpha_d = 1;
        
        x_d = moleFractions(gas);
        R_d = gasconstant/meanMolecularWeight(gas);
    else
        T_d = 0;
        M_d = 0;
        alpha_dc = 0;
        alpha_d = 0;
        x_d = x_u;
        R_d = 0;
    end
    
else
    T_d = 0;
    M_d = 0;
    alpha_dc = 0;
    alpha_d = 0;
    x_d = x_u;
    R_d = 0;
end
