function data_d = Sup_Comb(data_u, h, T, M)
%Sup_Comb exit properties after combustion

gas = GRI30;
set(gas, 'T', data_u.T, 'P', data_u.P, 'X', data_u.x);
v_u = data_u.M*soundspeed(gas);

equilibrate(gas, 'UV');

T_u = temperature(gas);
P_u = pressure(gas);
h_u = enthalpy_mass(gas);
M_u = v_u/soundspeed(gas);

if nargin == 1
    X0=[h_u;
        T_u;
        M_u];
else
    X0 = [h;
        T;
        M];
end
sol = fsolve(@(x)Comb_eqn(x, gas, v_u, h_u, T_u), X0);
T_d = sol(2);
M_d = sol(3);

setTemperature(gas, T_d)
v_d = M_d*soundspeed(gas);
P_d = (T_d/T_u)*(v_u/v_d)*P_u;

data_d.T = T_d;
data_d.P = P_d;
data_d.x = moleFractions(gas);
data_d.h = enthalpy_mass(gas);
data_d.M = M_d;
data_d.alpha = 1;

end

