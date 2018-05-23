function [F,J] = SC_eqn(X, gas, B, C, R)
%SC_eqn

setTemperature(gas, X(2));

F = [X(1)-intEnergy_mass(gas);
    X(3)^2-B*X(3)+R*X(2);
    X(1)+R*X(2)+X(3)^2/2-C];

J = [1 -cv_mass(gas) 0;
    0 R 2*X(3)-B;
    1 R X(3)];
end

