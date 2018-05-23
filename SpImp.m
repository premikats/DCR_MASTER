function [st,Isp] = SpImp(M_e, M_0, a_e, a_0, gamma_0, Alpha_e, Pi_e, f, g)
%SpImp Specific Impulse

%   Alpha_e = A_e/A_0
%   Pi_e = P_e/P_0

%   M - Mach Number
%   T - temperature
%   A - Area
%   P - Pressure
%   gamma - specific heat ratio
%   f - fuel air ratio
%   g - acceleration due to gravity

%   ()_e - Exit conditions
%   ()_0 - Freestream conditions

st=(1+f)*M_e*a_e - M_0*a_0 + Alpha_e*(1/M_0)*(a_0/gamma_0)*(Pi_e-1); % (m/s)
Isp=st/(f*g); % (s)

end

