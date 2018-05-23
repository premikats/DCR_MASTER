function [P,T,rho,gamma,R] = Stdatm(h)
%Standard atmosphere model
%   This model provides the value of ambient temperature and pressure for
%   the given altitude in m. This model is valid from (0-100)Km altitude.
% I/P - h-altitude in metres
% O/P - P-Pressure in Pa, Temperature in K, density(rho) in kg/m^3, mu -
% Absolute viscosity in Ns/m^2, k - thermal conductivity in J/(s.m.K)
a1=-6.5e-3;
a2=3e-3;
a3=-4.5e-3;
a4=4e-3;
if  h>100000
    disp('This atmospheric model cannot work for more than 100 Km')
    return
end
gas = GRI30;
io2 = speciesIndex(gas, 'O2');
in2 = speciesIndex(gas, 'N2');
iar = speciesIndex(gas, 'AR');
x = zeros(nSpecies(gas), 1);
x(io2, 1) = 0.21;
x(in2, 1) = 0.78;
x(iar, 1) = 0.01;
setMoleFractions(gas, x);
R = gasconstant/meanMolecularWeight(gas);

if h>=0&&h<11000
    T=288.16+a1*h;
    P=1.01325e5*((T/288.16)^(-9.81/(a1*R)));
    rho=1.2250*(T/288.16)^(-(9.81/(a1*R)+1));
elseif h>=11000&&h<25000
    T=216.66;
    P=2.2623e4*exp(-((9.81/(R*T))*(h-11000)));
    rho=0.3638*exp(-((9.81/(R*T))*(h-11000)));
elseif h>=25000&&h<47000
    T=216.66+a2*(h-25000);
    P=2.4861e3*((T/216.66)^(-9.81/(a2*R)));
    rho=0.040*(T/216.66)^(-(9.81/(a2*R)+1));
elseif h>=47000&&h<53000
    T=282.66;
    P=1.202214e2*exp(-((9.81/(R*T))*(h-47000)));
    rho=0.0015*exp(-((9.81/(R*T))*(h-47000)));
elseif h>=53000&&h<79000
    T=282.66+a3*(h-53000);
    P=58.2023*((T/282.66)^(-9.81/(a3*R)));
    rho=7.2619e-4*(T/282.66)^(-(9.81/(a3*R)+1));
elseif h>=79000&&h<=90000
    T=165.66;
    P=1.0063*exp(-((9.81/(R*T))*(h-79000)));
    rho=2.1423e-5*exp(-((9.81/(R*T))*(h-79000)));
elseif h>90000&&h<=100000
    T=165.66+a4*(h-90000);
    P=0.1040*(T/165.66)^(-9.81/(a4*R));
    rho=2.2150e-6*(T/165.66)^(-(9.81/(a3*R)+1));
end
set(gas, 'T', T, 'P', P);
gamma = cp_mass(gas)/cv_mass(gas);
% mu=1.46e-6*(T^1.5/(T+111));
% k=1.99e-3*(T^1.5/(T+112));
end



