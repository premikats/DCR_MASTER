function [Md,Td,beta,PR,TPR,TR,AR] = TPG_OSW(Mu,Tu,theta,R)
% This function evaluates the oblique shock angle when the inlet Mach number, inlet Static temperature and flow deflection angle are given
%for a thermally perfect gas.
% M1, T1 - Upstream Mach number, Static temperature
% theta-Flow deflection angle
% R=287.06 kJ/(Kg-K)
% M2,T2 - Downstream Mach number, Static Temperature
% beta-shock angle,TPR-Pt,u/Pt,d, TR-Td/Tu, AR-area Ratio
theta=theta*pi/180; % Convert degrees to radians
if theta<=1e-5
   disp('Solution cannot be evaluated for theta<1e-5 (in radians)')
   beta=0; 
   Md = 0;
   TPR = 0;
   PR=0; 
   TR=0;
   AR=0;
   Td=0;
end
if Mu<=1
    disp('Solution cannot be evaluated since Mu<=1')
    beta=0; 
    Md = 0;
    TPR = 0;
    PR=0; 
    AR=0;
    TR=0;
    Td=0;    
end
%% Evaluation of Cp and Tt in the  upstream - 
[H_static,Cp_air,gammau,A_air] = TPG_Cp(Tu);
[Ht,Tt,Vu,Pt_u_R,Tt_R,T_star,u_star,AR_R_star,Pt_star_R] = TPG_HtandTt(Mu,Tu,R);
[Ht,Cp_air_t,gamma_t,A_t_air] = TPG_Cp(Tt);
betamax=89*pi/180;
betamin=theta;
t=1; % Setting a dummy tolerance 
i=1;
p=1;
%% while loop for Td and beta calculation
while (abs(t)>1e-6)
    p=betamax-betamin;
%     if (abs(p)<1e-6)
%         disp('Shock is detached');
%         google=0;
%         c=0;
%         break
%     end    
    if i>200
       disp('Shock is detached'); 
       google=0;
       c=0;
       break
    end
    google=(betamax+betamin)/2;
    c=Vu^2*(cos(google)^2)*(cos(google-theta))^-2;
    d=sqrt(c);
    Td = TPG_T(Ht,Tt,d);
    t=R*(Tu/(Vu*sin(google))-Td/(sqrt(c)*sin(google-theta)))+Vu*sin(google)-sqrt(c)*sin(google-theta);
                    if t>0
                    betamax=(betamax+google)/2;
                    else
                        betamin=(betamin+google)/2;
                    end
                    i=i+1;
end
    beta=google;   %Computed Shock angle
    Vd=sqrt(c);
    if beta==0
       Md=0;
       Td=0;
       PR=0;
       TPR=0;
       TR=0;
       AR=0;
    else
        [H_static,Cp_air,gamma2,A_air] = TPG_Cp(Td);
        Md=Vd/sqrt(gamma2*R*Td);
        Vnd=Vd*sin(beta-theta);
        Vnu=Vu*sin(beta);
        PR=1+gammau*Mu^2*sin(beta)^2*(1-(Vnd/Vnu));
        [Ht_d,Tt_d,ud,Pt_d_R,Tt_d_R,T_d_star,u_d_star,AR_R_d_star,Pt_star_d_R] = TPG_HtandTt(Md,Td,R);
        TPR=(PR*Pt_d_R)/Pt_u_R; %Pt,d/Pt,u
        TR=Td/Tu;
        AR=sin(beta-theta)/sin(beta);
    end
    beta=beta*180/pi;
end


