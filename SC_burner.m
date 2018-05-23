function [T_d, P_d, x_d, M_d, alpha_dr] = SC_burner(T_ur, P_ur, x_ur, M_ur, R_ur, T_us, P_us, x_us, M_us, alpha_usc, R_us, ra, f)
%SC_burner Gives the exit properties across the supersonic combustor

gas_ur = GRI30;
set(gas_ur, 'T', T_ur, 'P', P_ur, 'X', x_ur)
u_ur = intEnergy_mass(gas_ur);
v_ur = M_ur*soundspeed(gas_ur);
rho_ur = density(gas_ur);

gas_us = GRI30;
set(gas_us, 'T', T_us, 'P', P_us, 'X', x_us)
u_us = intEnergy_mass(gas_us);
v_us = M_us*soundspeed(gas_us);
rho_us = density(gas_us);

alpha_dr = 1 + ((ra/(1-ra+f))*(rho_ur/rho_us)*(v_ur/v_us))/alpha_usc;
B = ((1-ra+f)/(1+f))*value(R_ur, T_ur, v_ur) + (ra/(1+f))*value(R_us, T_us, v_us);
C = ((1-ra+f)/(1+f))*(u_ur+R_ur*T_ur+v_ur^2/2) + (ra/(1+f))*(u_us+R_us*T_us+v_us^2/2);

X0 = [average(u_ur, u_us);
    average(T_ur, T_us);
    average(v_ur, v_us)];
x_d = x_ur+x_us;
gas_n = GRI30;
setMoleFractions(gas_n, x_d);
i = 1;
R(i) = gasconstant/meanMolecularWeight(gas_n);

while i>0
%     sol = fsolve(@(x)SC_eqn(x, gas_n, B, C, R(i)), X0);
    sol = newton(@SC_eqn, X0, gas_n, B, C, R(i), 1e-3);
    u_n = sol(1);
    T_n = sol(2);
    v_n = sol(3);
    rho_n = ((rho_ur*v_ur)/alpha_dr + rho_us*v_us - (rho_us*v_us)/alpha_dr)/v_n;
    set(gas_n, 'U', u_n, 'Rho', rho_n)
    equilibrate(gas_n, 'UV');
    
    i = i+1;
    R(i) = gasconstant/meanMolecularWeight(gas_n);
    tol = R(i)-R(i-1);
    if abs(tol)<1e-5
        break;
    end
end

gas_d = GRI30;
set(gas_d, 'U', u_n, 'Rho', rho_n, 'X', x_d)
equilibrate(gas_d, 'UV');

T_d = temperature(gas_d);
P_d = pressure(gas_d);
x_d = moleFractions(gas_d);
M_d = v_n/soundspeed(gas_d);

end

