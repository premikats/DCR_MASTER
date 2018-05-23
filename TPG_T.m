function T = TPG_T(Ht,Tt,c)
%This function calculates iterativelt the static temperature, when
%stagnation temperature and flow velocity are given foor TPG assumption.
% I/P - Ht - Total enthalpy, Tt-Total temperature, c - flow velocity
T_min=220;
T_max=Tt;
i=1;
t=1;
while (abs(t)>1e-6)
    T_mean=(T_max+T_min)/2;
    p=(T_max-T_min)/2;
    if (abs(p)<1e-7)
        disp('The static temperature could be erroneous')
        break
    end
    if i>500
        disp('The solutions could not be found in the required number of iterations')
        T_mean=0;
        break
    end
    [H_mean,Cp_air_mean,gamma_mean,A_mean_air] = TPG_Cp(T_mean);
    B=H_mean+(c^2/2);
    t=B-Ht;
                    if t>0
                        T_max=(T_mean+T_max)/2;
                    else
                         T_min=(T_mean+T_min)/2;
                    end
                        i=i+1;
end
T=T_mean;
end

