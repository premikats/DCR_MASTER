function [T_d, P_d, x_d, M_d, alpha_d, alpha_t] = Inlet_sub(T_u, P_u, x_u, M_u)
%Inlet_sub Exit properties after subsonic compression

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
h_u = enthalpy_mass(gas);
s_u = entropy_mass(gas);
gamma_u = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);
hx0 = h_u+(gamma_u*R*T_u*M_u^2)/2;
alphastar_u = A_Astar_TPG(x_u, M_u, T_u);
PI_yu = NormalShock_TPG(x_u, M_u, T_u, 'Pt');
alpha_ux = alphastar_u*PI_yu;

alphastar_t = 1/PI_yu;

%Isentropic Relations between x and u:
Mx_CPG = fzero(@(M)(A_Astar_CPG(M,gamma_u)-alphastar_t),M_u);
taux_CPG = Isen_CPG(M_u, gamma_u, 'T')/Isen_CPG(Mx_CPG, gamma_u, 'T');
pix_CPG = Isen_CPG(M_u, gamma_u, 'P')/Isen_CPG(Mx_CPG, gamma_u, 'P');

% CPG initial conditions for faster convergence:
xx0 = [Mx_CPG;
    taux_CPG;
    pix_CPG];
% solx = fsolve(@(x)Throat_isen(x, gas, M_u, T_u, P_u, s_u, gamma_u, hx0, alpha_ux), xx0, options);
solx = newton(@Throat_isen_newt, xx0, gas, M_u, T_u, P_u, s_u, gamma_u, hx0, alpha_ux);

M_x = solx(1);
tau_xu = solx(2);
pi_xu = solx(3);

% Normal Shock relations between x and y:
T_x = tau_xu*T_u;
P_x = pi_xu*P_u;
M_y = NormalShock_TPG(x_u, M_x, T_x,'M');
pi_yx = NormalShock_TPG(x_u, M_x, T_x,'P'); % Py/Px
tau_yx = NormalShock_TPG(x_u, M_x, T_x,'T'); % Ty/Tx
alpha_yx = 1; % Ay/Ax

T_y = tau_yx*T_x;
P_y = pi_yx*P_x;
set(gas, 'T', T_y, 'P', P_y, 'X', x_u);
gamma_y = cp_mass(gas)/cv_mass(gas);
s_y = entropy_mass(gas);
h_y = enthalpy_mass(gas);
R = gasconstant/meanMolecularWeight(gas);
hy0 = h_y+(gamma_y*R*T_y*M_y^2)/2;

% Isentropic relations between y and 4r:
M_d = 0.7*M_y;
alphastar_y = A_Astar_CPG(M_y,gamma_y);
alphastar_d = A_Astar_CPG(M_d,gamma_y);
alphad_CPG = alphastar_d/alphastar_y; % Ad/Ay
taud_CPG = Isen_CPG(M_y,gamma_y,'T')/Isen_CPG(M_d,gamma_y,'T');
pid_CPG = Isen_CPG(M_y,gamma_y,'P')/Isen_CPG(M_d,gamma_y,'P');

xd0 = [taud_CPG;
       pid_CPG;
       alphad_CPG];
options = optimoptions('fsolve', 'MaxFunEvals', 100000, 'MaxIter', 100000);
sold = fsolve(@(x)Diff_isen(x, gas, M_y, T_y, P_y, s_y, gamma_y, hy0, M_d), xd0, options);

tau_dy = sold(1);
pi_dy = sold(2);
alpha_dy = sold(3);

% Final Ratios:
alpha_du = (alpha_dy*alpha_yx)/alpha_ux;
tau_du = tau_dy*tau_yx*tau_xu;
pi_du = pi_dy*pi_yx*pi_xu;

T_d = tau_du*T_u;
P_d = pi_du*P_u;
x_d = x_u;
alpha_d = alpha_du;
alpha_t = 1/alpha_ux;

end

