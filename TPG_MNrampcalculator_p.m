% function [Me,Te,Md_e,Md_i,Mx,PR_e,PR_i,TR_e,TR_i,beta_e,beta_i,OTPR,OPR,OTR,OAR_Kant,OAR_geo,delOAR,OAR_Isent,f] = TPG_MNrampcalculator(M,N,T0,M0,Theta_e,Theta_i,R)   
function [Te, Pe, xe, Me, OAR_geo, beta_e, beta_i, f] =  TPG_MNrampcalculator_p(T0,P0,x0,M0,M,N,Theta_e,Theta_i)
%This is a general 2D ramp calculator function. This function evaluates all
%the exit quantities.The program is designed to achieve
%I/P - M=no. external compression ramps, N=no.internal compression ramps,T0
%- freestream static temperature,M0=freestream Mach number, Theta_e- array of flow deflection angles where
%every element is the magnitude of flow deflection angle at each external
%compression ramp, Theta_i=array of flow deflection angles for internal
%ramps
%O/P - Me=Exit Mach number, OTPR-Overall Tot Pressure Ratio, OAR_Kant,
%OAR_geo, Mx- Mach number at start of internal compression, beta_e, beta_i,
%OTR,, PR_e, PR_i;
if M0<1
    error('This program cannot calculate inlet quantities for M0<1')
end
if M~=length(Theta_e)
    error('Number of elements in the Theta_e does not match M')
end
if N~=length(Theta_i)
    error('Number of elements in Theta_i does not match N')
end
%theta_i=(M*Theta_e)/N; % decides magnitude of flow deflection in internal ramp based on external ramp angle
gas = GRI30;
set(gas, 'T', T0, 'P', P0, 'X', x0);
R = gasconstant/meanMolecularWeight(gas);

[Ht_0,Tt_0,u0,Pt_R_0,Tt_R_0,T_star_0,u_star_0,AR_R_star_0,Pt_star_R_0] = TPG_HtandTt(M0,T0,R);
OAR_Isent=AR_R_star_0;
OTPR=1;
OPR=1;
OTR=1;
Ar_ext=1;
Ar_int=1;
i=1; % Initializing to run the Theta_e in the loop
j=1; % Initializing to run the Theta_i in the loop
theta_c=0;
p=0;
r=0;
Mu=M0;
Tu=T0;
while(p<M && Mu>1)
    if p==0
       [Md_e(i),Td_e(i),beta_e(i),PR_e(i),TPR_e(i),TR_e(i),AR_e(i)] = TPG_OSW(Mu,Tu,Theta_e(i),R);% No BL in the first ramp
    else
        [Md_e(i),Td_e(i),beta_e(i),PR_e(i),TPR_e(i),TR_e(i),AR_e(i)] = TPG_OSW_BLS(Mu,Tu,Theta_e(i),R);
    end
    Ar_ext=Ar_ext*AR_e(i);
    Mu=Md_e(i);
    Tu=Td_e(i);
    OTPR=OTPR*TPR_e(i);
    OPR=OPR*PR_e(i);
    OTR=OTR*TR_e(i);
    p=p+1;
    i=i+1;
end
Mx=Md_e(i-1);
Tx=Td_e(i-1);
if Mx<=1
    OAR_Kant=0;
    beta_i=zeros(1,N);
    Md_i=zeros(1,N);
    PR_i=zeros(1,N);
    TR_i=zeros(1,N);
    Me=0;
    Te=0;
    OPR=0;
    OTPR=0;
    OTR=0;
    Ar_int=0;
else
    AR_Kant = TPG_AR_Kant(Mx,Tx,R);
    OAR_Kant = AR_Kant*Ar_ext;
    Mi=Mx;
    Ti=Tx;
    while(r<N && Mi>1)
        [Md_i(j),Td_i(j),beta_i(j),PR_i(j),TPR_i(j),TR_i(j),AR_i(j)] = TPG_OSW_BLS(Mi,Ti,Theta_i(j),R);
        if beta_i(j)==0
            if j==1
                beta_i=[beta_i zeros(1,N-j)];
                Md_i=[zeros(1,N)];
                PR_i=[zeros(1,N)];
                OPR=0;
                OTPR=0;
                OTR=0;
                Ar_int=0;
                Me=0;
                Te=0;
                break
            else
                beta_i=[beta_i zeros(1,N-j)];
                Md_i=[Md_i zeros(1,N-j)];
                PR_i=[PR_i zeros(1,N-j)];
                OPR=0;
                OTPR=0;
                OTR=0;
                Ar_int=0;
                Me=0;
                Te=0;
                break
            end
        end
        Mi=Md_i(j);
        Ti=Td_i(j);
        OTPR=OTPR*TPR_i(j);
        OPR=OPR*PR_i(j);
        OTR=OTR*TR_i(j);
        if Mi~=0
            if r==N-1
               Ar_int=Ar_int*sind(beta_i(j)-Theta_i(j))/sind(sum(Theta_e)+beta_i(j)-sum(Theta_i));
            else
               Ar_int=Ar_int*AR_i(j);
            end
        else
            Ar_int=0;
        end
    r=r+1;
    j=j+1;
end
if OTPR~=0
    Me=Mi;
    Te=Ti;
end
end
OAR_geo=Ar_ext*Ar_int;
delOAR=OAR_geo-OAR_Kant;
if OTPR==0
    f=1;
else
    f=0;
end

Pe = OPR*P0;
set(gas, 'T', Te, 'P', Pe);
xe = moleFractions(gas);
end




