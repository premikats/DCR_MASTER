function [Ht,Tt,u,Pt_R,Tt_R,T_star,u_star,AR_R_star,Pt_star_R] = TPG_HtandTt(M,T,R)
%This function calculates all the stagnation and static quantities and also
%sonic quantities for isentropic flow assumtions for TPG condition.
% NASA air-polynomials are used to evaluate the Enthalpies at temperatures
% greater than 200K and stagnation temperatures lower than 6000K. I/P -
% M=Mach number, T-Static Temp,R-gas constant=287.06 kJ/kg-K; O/P -
% Ht-Tot.enthalpy, Tt-Total Temp,u-velocity,Pt_R-Pt/P
% Tt_R=Tt/T, T_star-Sonic Temp for given T(T*),u_star-critical
% velocity(u*),Ar_R_star-A*/A,A*-critical area, Pt_star_R=Pt/P*
[h,Cp_air,gamma,A_air] = TPG_Cp(T);
u=M*sqrt(gamma*R*T);
Ht= h+u^2/2; % Calculation of stagnation enthalpy
Tt_max=T*(1+((gamma-1)/2)*M^2);
Tt_min=T; %0.5*Tt_max;
t=1; % Dummy Tolerance
i=1;
while(abs(t)>1e-5)
       Tt_mean=(Tt_max+Tt_min)/2;
       p=Tt_max-Tt_min;
       if i>1000 && abs(t)>1e-5
          disp('The solution could not be found in the required number of iterations');
          Tt_mean=0;
          break
       end
         if p<1e-7 && abs(t)>1e-5
             disp('The solution could be erroneous');
             break
         end
         [H_mean,Cp_air_mean,gamma_mean,A_air_mean] = TPG_Cp(Tt_mean);
         t=H_mean-Ht;
                    if t>0
                        Tt_max=(Tt_mean+Tt_max)/2;
                    else
                         Tt_min=(Tt_mean+Tt_min)/2;
                    end
                        i=i+1;
end
Tt=Tt_mean;
%% Expression to calculate Pt/P and rho_t/rho for a given static 
[h_t,Cp_mean_t,gamma_t,A_t_air] = TPG_Cp(Tt);
if Tt<=1000
    B=sum([-0.5*(Tt^-2-T^-2) -(Tt^-1-T^-1) log(Tt/T) (Tt-T) 0.5*(Tt^2-T^2) (Tt^3-T^3)/3 (Tt^4-T^4)/4 0].*A_t_air);
else
    if T<1000
        B=sum([-0.5*(Tt^-2-1000^-2) -(Tt^-1-1000^-1) log(Tt/1000) (Tt-1000) 0.5*(Tt^2-1000^2) (Tt^3-1000^3)/3 (Tt^4-1000^4)/4 0].*A_t_air)+sum([-0.5*(1000^-2-T^-2) -(1000^-1-T^-1) log(1000/T) (1000-T) 0.5*(1000^2-T^2) (1000^3-T^3)/3 (1000^4-T^4)/4 0].*A_air);
    else
        B=sum([-0.5*(Tt^-2-T^-2) -(Tt^-1-T^-1) log(Tt/T) (Tt-T) 0.5*(Tt^2-T^2) (Tt^3-T^3)/3 (Tt^4-T^4)/4 0].*A_t_air);
    end
end
Pt_R=exp(B/R); %Pt/P
Tt_R=Tt/T; % Ratio of total temp to static temp
Rho_R=Tt_R/Pt_R;
%% Evaluation of Sonic Temperature -
if M>1
    T_star_max=Tt; 
    T_star_min=T;
else
    T_star_max=T; 
    T_star_min=200;
end
r=1; % dummy tolerance
q=1; % Iteration counter 
while(abs(r)>1e-5)
    T_star_mean=(T_star_max+T_star_min)/2;
    s=T_star_max-T_star_min;
         if q>1000 && abs(r)>1e-5
          disp('The solution could not be found in the required number of iterations')
          T_star_mean=0;
          break
         end
         if s<1e-5 && abs(r)>1e-5
             disp('The solution could be erroneous')
             break
         end
         [H_static_star_mean,Cp_air_star_mean,gamma_star_mean,A_air_star_mean] = TPG_Cp(T_star_mean);
         r=gamma_star_mean*R*T_star_mean-2*(Ht-H_static_star_mean);   
                    if r>0
                        T_star_max=(T_star_mean+T_star_max)/2;
                    else
                         T_star_min=(T_star_mean+T_star_min)/2;
                    end
                        q=q+1;
end
   T_star=T_star_mean;
%% Expression to calculate Pt/P and rho_t/rho* for the sonic state 
[h_star,Cp_star,gamma_star,A_star_air] = TPG_Cp(T_star);
C=sum([-0.5*(Tt^-2-T_star^-2) -(Tt^-1-T_star^-1) log(Tt/T_star) (Tt-T_star) 0.5*(Tt^2-T_star^2) (Tt^3-T_star^3)/3 (Tt^4-T_star^4)/4 0].*A_star_air);
Pt_star_R=exp(C/R); %Pt/P*
Tt_star_R=Tt/T_star; 
u_star=sqrt(gamma_star*R*T_star);
Rho_star_R=Tt_star_R/Pt_star_R;
%% Evaluation of AR_R_star=A*/A  
AR_R_star = (Rho_R*u)/(Rho_star_R*sqrt(gamma_star*R*T_star));
end

