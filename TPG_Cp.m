function [H_static,Cp_air,gamma,A_air] = TPG_Cp(T)
%Thermally perfect gas modelling of Air for given altitude and Static
%Temperature. The mass fraction of all the species have been assumed to be
%constant across all altitudes till 100km. The mass fraction of all the species will
%be varied in the next code iteration. Cp of the individual species have
%been implemented as a 4th degree 7 term polynomial fit. 8th term is the 
%integration constant, which will be used for calculation of static enthalpy at a given temperature.  
% Static Enthalpy is calculated from reference temperature of 298.15K.
% The polynomials are valid for temperatures above 200 K. For temperatures
% below 200K, Cp(200K)is assumed, and all the O/P parameters are calculated
% using CPG assumptions.
%% Calculation of R for air (in J/mol-K)-
R=287.06;
if T<200
   Cp_air=999.0548;
   gamma=Cp_air/(Cp_air-R);
   H_static=Cp_air*T-302110.96;
   A_air=[0];
else
%% Universal Gas constant
    R_univ=8.31447;
%% Mass Fractions of each species - (Source-not available,Has to be verified)
    mf_O2=0.2320;
    mf_N2=0.7547;
    mf_CO2=4.6e-4;
    mf_H2=0; %(Approximated due to very low molecular weights)
    mf_Ar=0.0128;
    mf_Ne=1.2e-5;
    mf_He=7e-7;
    mf_Kr=3e-6;
    mf_Xe=4e-7;
%% Molecular weight of Species in Kg/mol
    m_O2=31.9988e-3;
    m_N2=28.1348e-3;
    m_CO2=44.0098e-3;
    m_H2=2.01588e-3; 
    m_Ar=39.48e-3;
    m_Ne=20.1797e-3;
    m_He=2.01588e-3;
    m_Kr=83.80e-3;
    m_Xe=131.29e-3;
% R=8.31447*((mf_O2/m_O2)+(mf_N2/m_N2)+(mf_CO2/m_CO2)+(mf_H2/m_H2)+(mf_Ar/m_Ar)+(mf_Ne/m_Ne)+(mf_He/m_He)+(mf_Kr/m_Kr)+(mf_Xe/m_Xe));
%% Gas constant for various gas species (in J/kg-K)-
    R_O2=R_univ/m_O2;
    R_N2=R_univ/m_N2;
    R_CO2=R_univ/m_CO2;
    R_H2=R_univ/m_H2; 
    R_Ar=R_univ/m_Ar;
    R_Ne=R_univ/m_Ne;
    R_He=R_univ/m_He;
    R_Kr=R_univ/m_Kr;
    R_Xe=R_univ/m_Xe;
%% Temperature Range for calculation
%% Temperature Matrix
        if T>=200&&T<=1000
%% Temperature Coefficients
            A_H2=[4.07832281E4 -8.00918545E2 8.21470167 -1.26971436E-2 1.75360493E-5 -1.20286016E-8 3.36809316E-12 2.682482438E3];
            A_N2=[2.21037122E4 -3.81846145E2 6.08273815 -8.53091381E-3 1.3846461E-5 -9.6257929293E-9 2.5197056E-12 7.10845911E2];
            A_CO2=[4.94365054E4 -6.26411601E2 5.30172524 2.503813816E-3 -2.127308728E-7 -7.689988780E-10 2.849677801E-13 -4.52819846E4];
            A_O2=[-3.42556269E4 4.8469986E2 1.11901159 4.29388743E-3 -6.83627313E-7 -2.02337478E-9 1.03904064E-12 -3.39145434E3];
            A_He=[0 0 2.5 0 0 0 0 -745.375];
            A_Ne=[0 0 2.5 0 0 0 0 -745.375];
            A_Ar=[0 0 2.5 0 0 0 0 -745.375];
            A_Kr=[0 0 2.5 0 0 0 0 -745.375];
            A_Xe=[6.60802392E-3 -9.53610408E-5 2.50000053 -1.49716621E-9 2.21314503E-12 -1.64711078E-3 4.84969606E-19 -745.374544];
        else T>1000&&T<6000
            A_H2=[5.6018338E5 -8.37149134E2 2.97536304 1.25224993E-3 -3.74071842E-7 5.93662825E-11 3.36809316E-12 5.33981585E3];
            A_N2=[5.87709908E5 -2.23924255E3 6.06694267 -6.13965296E-4 1.49179819E-7 -1.92309442E-11 1.06194871E-15 1.28320618E4];
            A_CO2=[1.176962419E5 -1.788791477E3 8.29152319 -9.22315678E-5 4.86367688E-9 -1.891053312E-12 6.33003659E-16 -3.90835059E4];
            A_O2=[-1.03793994E6 2.34483275E3 1.81972949 1.26784887E-3 -2.18807142E-7 2.05372411E-11 -8.19349062E-16 -1.68901253E4];
            A_He=[0 0 2.5 0 0 0 0 -745.375];
            A_Ne=[0 0 2.5 0 0 0 0 -745.375];
            A_Ar=[0 0 2.5 0 0 0 0 -745.375];
            A_Kr=[0 0 2.5 0 0 0 0 -745.375];
            A_Xe=[1104.19906 -3.31828724 2.50387762 -2.25191164E-6 6.86935273E-10 -1.0484943E-13 6.29438941E-18 -724.304574];
        end
%% Temperature Coefficients for Air
    A_air=(R_H2*mf_H2*A_H2)+(R_CO2*mf_CO2*A_CO2)+(R_N2*mf_N2*A_N2)+(R_O2*mf_O2*A_O2)+(R_He*mf_He*A_He)+(R_Ar*mf_Ar*A_Ar)+(R_Kr*mf_Kr*A_Kr)+(R_Ne*mf_Ne*A_Ne);%+(R_Xe*mf_Xe*A_Xe);
    B=[T.^-2 T.^-1 1 T T.^2 T.^3 T.^4 0];
    C=[-T.^-1 log(T) T T.^2/2 T.^3/3 T.^4/4 T.^5/5 1];
    Cp_air=sum(A_air.*B);
    gamma=Cp_air/(Cp_air-R);
    H_static=sum(A_air.*C);
%% Polynomial fit for air (Source-not verifiable)
    Cp_poly = 1.9327E-10*T^4 - 7.9999E-07*T^3 + 1.1407E-03*T^2 - 4.4890E-01*T + 1.0575E+03;
    gamma_poly=Cp_poly/(Cp_poly-R);
end
end
