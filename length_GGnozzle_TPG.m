function [l_c, l_d, theta_d, alpha_su] = length_GGnozzle_TPG(T_u, P_u, x_u, M_u, A_u, T_d, P_d, x_d, M_d, A_d, theta_c, eta)
%Non-isentropic nozzle sizing: (Ref:Zucrow)

gas = GRI30;
set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
gamma_u = cp_mass(gas)/cv_mass(gas);
set(gas, 'T', T_d, 'P', P_d, 'X', x_d);
gamma_d = cp_mass(gas)/cv_mass(gas);
gamma = average(gamma_u, gamma_d);

Pt = P_u*Formula(gamma,M_u)^(gamma/(gamma-1));
Tt = T_u*Formula(gamma,M_u);

eta_n = (1-1/Formula(gamma,M_d))/(1-1/Formula(gamma,M_u)-(1/Formula(gamma,M_d)-1/Formula(gamma,M_u))/eta); %Equivalent Kinetic energy efficiency
A = eta_n-1;
B = (3*eta_n-((3*gamma-1)/(gamma-1)))/(gamma-1);
C = (2*eta_n)/(gamma-1)^2;
x1 = sqrt((-B-sqrt(B^2-4*A*C))/(2*A));
x2 = sqrt((-B+sqrt(B^2-4*A*C))/(2*A));
if imag(x1)==0
    Mstar = x1;
elseif imag(x2)==0
    Mstar = x2;
else
    error('No minimum throat area');
end

Pstar = Pt*(1-((Formula(gamma,Mstar)-1)/(eta_n*Formula(gamma,Mstar))));
Tstar = Tt/Formula(gamma,Mstar);
alpha_su = (P_u/Pstar)*(M_u/Mstar)*sqrt(Tstar/T_u); % A/Astar for the converging portion
alpha_sd = (Pstar/P_d)*(Mstar/M_d)*sqrt(T_d/Tstar); % A/Astar for the diverging portion

l_c = A_u*(1-alpha_su)/tand(theta_c);
nu_d = Prandlt_Meyer_CPG(M_d, gamma);
theta_d = nu_d/2; %(Ref: High Speed Wind Tunnel Testing, Alan Pope)
% theta_d = 2;
l_d = A_d*(1-1/alpha_sd)/tand(theta_d);

end

