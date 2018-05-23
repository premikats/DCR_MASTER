function [T_d, P_d, x_d, M_d, R_d, a_d] = Exit_Nozzle(T_u, P_u, x_u, M_u, alpha, eta)
%Exit_nozzle exit paramaeters of the nozzle for a given area ratio

P_d = Bisection(@nozz_eqn, 0.01*oneatm, oneatm, T_u, P_u, x_u, M_u, eta, alpha, 1e-5);
[T_d, x_d, M_d, alpha_d, R_d, a_d] = nozzle(T_u, P_u, x_u, M_u, P_d, eta);
end

