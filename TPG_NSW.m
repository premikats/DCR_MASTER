function [Md,Td,Temp_ratio,Den_ratio,PR] = TPG_NSW(Mu,Tu,R)
%Function to evaluate the downstream quantities of a normal shock for
%thermally perfect gas conditions (TPG).
% Inputs for the function, all the upstream quantities,Tu=Upstream Temperature,
% Mu = Upstream Mach number, R-gas constant, 287.06 kJ/kg-K 
% O/P for the function are all the downstream quantities, Md=Downstream
% Mach number, Td=downstream temp, Temp_ratio=Td/Tu, Den_ratio=rho_d/rho_u,
% PR=Pd/Pu
if Mu<1
    error('Mach number cannot be lesser than 1 for Normal shocks')
    return
end
if Tu<200
    error('TPG code cannot calculate for Tu<200')
    return
end
[hu,Cp_u_air,gamma_u,A_u_air] = TPG_Cp(Tu);
[Ht,Tt,Vu] = TPG_HtandTt(Mu,Tu,R);
Q=(R*Tu/Vu)+Vu; % Mass-momentum conservation applied to the upstream quaantities
Td_max=Tt; % Search for Temperature which is always higher than upstream temperature
Td_min=Tu;
i=1;
t=1;
while(abs(t)>1e-5)
    Td_mean=(Td_max+Td_min)/2;
    p=Td_max-Td_min;
       if i>1000 && abs(p)>1e-5
          disp('The solution could not be found in the required number of iterations')
          Td_mean=0;
          break
%           return
       end
         if p<1e-5 && i<1000
             disp('The solution could be erroneous')
             Td=Td_mean;
             break
%              return
         end
    [hd,Cp_air_mean,gamma_mean,A_mean_air] = TPG_Cp(Td_mean);
    Vd=sqrt(2*(Ht-hd));
    t=(R*Td_mean)/Vd+Vd-Q;
                    if t>0
                        Td_max=(Td_mean+Td_max)/2;
                    else
                         Td_min=(Td_mean+Td_min)/2;
                    end
                        i=i+1;
end
Td=Td_mean;
Md=Vd/sqrt(gamma_mean*R*Td_mean);
Temp_ratio=Td/Tu;
Den_ratio=Vu/Vd; %rho_d/rho_u
PR=Den_ratio*Temp_ratio;
end
