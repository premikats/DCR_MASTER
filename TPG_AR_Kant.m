function [AR_Kant] = TPG_AR_Kant(Mx,Tx,R)
%The following function evealuates the Kantrowitz's limiting contraction
%ratio for Thermally perfect gas assumption. It assumes normal shock at Mx
%and an isentropic flow downstream of the shock.
% I/P - Mx=Mach number at start of internal compression, Tx=Static
% temperature at the point, R=gas constant=287.06 kJ/kg-K;
if Mx<1
    error('AR_Kant cannot be evaluated for M<1')
    return
end
[Md,Td,Temp_ratio,Den_ratio,PR] = TPG_NSW(Mx,Tx,R); % Calculate the downstream quantities behind the normal shock
[H_static_d,Cp_air_d,gamma_d,A_air_d] = TPG_Cp(Td);
[Ht_d,Tt_d,u,Pt_R,Tt_R,T_star,u_star,AR_R_star,Pt_star_R] = TPG_HtandTt(Md,Td,R); % Use isentropic relations to calculate other quantities
[H_static_star,Cp_air_star,gamma_star,A_air_star] = TPG_Cp(T_star);
B=Pt_star_R/Pt_R; %Pd/P*
C=sqrt(gamma_d/gamma_star);
D=sqrt(Td/T_star);
AR_Kant=Md*B*C*D;
end

