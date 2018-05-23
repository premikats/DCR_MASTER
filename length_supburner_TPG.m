function Length = length_supburner_TPG(P_ur, P_us, T_ur, T_us, A_d, phi, V_0)
% Length of Combustor: (Ref: Tiago C. Rolin and Frank K. Lu)

% For H2-Air
% Mixture Velocity and Temperature: (Length decreases with increase in
% temperature. Exit Temperature higher than the average temperature. Avg
% Temp assumed for safety)
p = 9.86923e-6*((P_ur+P_us)/2); % Avg Pressure in atm
T = (T_ur+T_us)/2; % Avg T in K
% T = max(T_ur,T_us);
t_ign = (8e-9)*(exp(9600/T)/p);
t_reac = 0.000105*(exp((-1.12*T)/1000)/p^1.7);
t_comb = t_ign+t_reac;
l_comb = V_0*t_comb;
Cm = 25; % Approx (Ref: Heiser and Pratt)
l_m = 0.179*Cm*exp(1.72*phi)*A_d; % for phi<=1
Length = l_comb+l_m;

end

