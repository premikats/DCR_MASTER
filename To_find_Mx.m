function F = To_find_Mx(x,gas,alphastar_t)
%T_find_Mx
% Equation to be solved to find the corresponding Mx to alphastar_t

F = [x(2)-temperature(gas);
     A_Astar_TPG(gas,x(1))-alphastar_t];


end

