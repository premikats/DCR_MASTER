function F = nozz_eqn(x, T_u, P_u, x_u, M_u, eta, alphae)
%nozz_eqn to obtain the Pressure corresponding to the given area ratio

F = nozzle_alpha(T_u, P_u, x_u, M_u, x, eta)-alphae;

end

