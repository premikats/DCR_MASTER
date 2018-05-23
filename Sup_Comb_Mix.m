function data_d = Sup_Comb_Mix(data_ur, data_us, ra, f, h, T, M)
%Sup_Comb Gives the exit properties after the supersonic combustor

T_us = data_us.T;
P_us = data_us.P;
x_us = data_us.x;
h_us = data_us.h;
gamma_us = data_us.gamma;
R_us = data_us.R;
M_us = data_us.M;
v_us = M_us*sqrt(gamma_us*R_us*T_us);

T_ur = data_ur.T;
P_ur = data_ur.P;
x_ur = data_ur.x;
h_ur = data_ur.h;
gamma_ur = data_ur.gamma;
R_ur = data_ur.R;
M_ur = data_ur.M;
v_ur = M_ur*sqrt(gamma_ur*R_ur*T_ur);

x_d = x_us + x_ur;

alpha_dr = 1 + (ra/(1-ra+f))*(P_ur/P_us)*(R_us/R_ur)*(T_us/T_ur)*(v_ur/v_us);

term = -((ra/(1+f))*(v_us+((R_us*T_us)/v_us))+((1-ra+f)/(1+f))*(v_ur+((R_ur*T_ur)/v_ur)));

if nargin == 4
    X0=[(h_ur+h_us)/2;
        (T_ur+T_us)/2;
        (M_ur+M_us)/2];
else
    X0 = [h;
        T;
        M];
end

gas = GRI30;
set(gas, 'X', x_d)

sol = fsolve(@(x)mix_eqn(x, gas, term, h_ur, h_us, v_ur, v_us), X0);
h_d = sol(1);
T_d = sol(2);
M_d = sol(3);

set(gas, 'T', T_d, 'X', x_d)
v_d = M_d*soundspeed(gas);
R_d = gasconstant/meanMolecularWeight(gas);
P_d = ((R_d*T_d)/v_d)*(((P_ur*v_ur)/(R_ur*T_ur))-((P_us*v_us)/(R_us*T_us)))*(1/alpha_dr)+((R_d*T_d)/v_d)*((P_us*v_us)/(R_us*T_us));
set(gas, 'T', T_d, 'P', P_d, 'X', x_d)

data_d.T = T_d;
data_d.P = P_d;
data_d.x = x_d;
data_d.h = h_d;
data_d.M = M_d;
data_d.alpha = alpha_dr;

end

