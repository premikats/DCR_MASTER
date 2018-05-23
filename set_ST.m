function P = set_ST(x, S, T_f, T_ini, P_ini)
%set_ST The pressure of the gas when entropy and temperature is given

if isa(x,'double')
    gas = GRI30;
    setMoleFractions(gas, x)
elseif isa(x,'Solution')
    gas = x;
else
    error('Invalid input argument data type')
end

T = T_ini;
P = P_ini;

if T_f>=T_ini
    while (T_f-T)>0.001;
        P = 1.001*P;
        set(gas, 'S', S, 'P', P);
        T = temperature(gas);
    end
else
    while (T-T_f)>0.001;
        P = 0.999*P;
        set(gas, 'S', S, 'P', P);
        T = temperature(gas);
    end
end

end

