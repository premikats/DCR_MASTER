function Length = length_GGburner_TPG(P_u, T_u, A_u, phi, V_0)
%Length of the Combustor: (Ref: Tiago C. Rolin and Frank K. Lu)

% For H2-Air
p = 9.86923e-6*P_u; % Pressure in atm
t_ign = ((8e-9)*exp(9600/T_u))/p;
t_reac = (0.000105*exp((-1.12*T_u)/1000))/(p^1.7);
t_comb = t_ign+t_reac;
l_comb = V_0*t_comb;
Cm = 60; % Approx (Ref: Heiser and Pratt)
if phi>1
    l_m = 3.333*Cm*exp(-1.204*phi)*A_u;
else
    l_m = 0.179*Cm*exp(1.72*phi)*A_u;
end
Length = l_comb+l_m;

end

