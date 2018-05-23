function [F] = Throat(x, gas, h0, alphastar_t)
%Throat 
% temperature and Mach number at throat for given alphastar_t

setTemperature(gas, x(2));
X = moleFractions(gas);
gamma = cp_mass(gas)/cv_mass(gas);
R = gasconstant/meanMolecularWeight(gas);

F = [x(3)+(gamma*R*x(2)*x(1)^2)/2-h0;
    x(3)-enthalpy_mass(gas);
    A_Astar_TPG(X, x(1), x(2))-alphastar_t];
% 
% setTemperature(gas, x(2));
% gamma = cp_mass(gas)/cv_mass(gas);
% 
% j11 = gamma*R*x(2)*x(1);
% j12 = (gamma*R*x(1)^2)/2;
% j13 = 1;
% j21 = 0;
% j22 = -cp_mass(gas);
% j23 = 1;
% j33 = 0;
% 
% a = A_Astar_TPG(X, x(1), x(2));
% aM = A_Astar_TPG(X, (x(1)+1e-5), x(2));
% j31 = (aM-a)/1e-5;
% 
% aT = A_Astar_TPG(X, x(1), (x(2)+1e-5));
% j32 = (aT-a)/1e-5;
% 
% J = [j11 j12 j13;
%     j21 j22 j23;
%     j31 j32 j33];
end

