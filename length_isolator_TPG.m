function Length = length_isolator_TPG(T_u, P_u, x_u, M_u, A_u, P_d, phi)
% Finding the length: (Ref: Byun et al.)
    
    gas = GRI30;
    set(gas, 'T', T_u, 'P', P_u, 'X', x_u);
    a_u = soundspeed(gas);
    pi_d = P_d/P_u;
    Dh = A_u;
    w = 0.75; % For low temperatures (for high temperatures, w = 0.5)
    T_stag = Isen_TPG(x_u, M_u, T_u, 'T')*T_u;
    P_stag = Isen_TPG(x_u, M_u, T_u, 'P')*P_u;
    set(gas, 'T', T_stag, 'P', P_stag, 'X', x_u);
    a_stag = soundspeed(gas);
    nu_stag = (1e-10)*T_stag^2+(3e-8)*T_stag-3e-6; % Variation of Kinematic Viscosity with temperature (curve fit)
    nu_u = (1e-10)*T_u^2+(3e-8)*T_u-3e-6;
    l = 10;
    tol = 1;
    while tol>1e-10
        x1 = (1+(M_u^2)/5)^(-3+w);
        Re_st = (a_stag/nu_stag)*l*M_u*x1;
        theta = 0.022*((1+0.16*M_u^2)^(-0.6))*l*(Re_st)^(-1/6); % Ref: Stratford and Beavers
        Re_theta = (a_u/nu_u)*M_u*theta;
        x2 = ((M_u^2)-1)/(M_u^2);
        x3 = 1.35*(pi_d-1)+1.45*(pi_d-1)^2;
        l_n = (phi^0.52)*Dh*x2*(Re_theta^0.25)*((Dh/theta)^0.5)*x3^-1;
        tol = abs((l_n-l)/l);
        l = l_n;
    end
    Length = l;
end