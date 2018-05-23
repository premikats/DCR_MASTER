function [st, Isp, M, T, P, x, alpha, beta, x_species, f] = Specific_thrust_TPG(M0, P0, T0, rho0, gamma0, R, x_air, inlet, ra, k, fuel, phi, pi_er, Alphae, ds)
% Specific_thrust_TPG
% Properties across each components of the engine

% [P0, T0, rho0, gamma0, R] = Stdatm(alt);
gas = GRI30;
set(gas, 'T', T0, 'P', P0, 'X', x_air);
a0 = soundspeed(gas);
MW = meanMolecularWeight(gas);

eta_er = 0.9;
eta_e = 0.9;

phi_gg = phi/(1-ra);
if strcmp(fuel, 'H2')
    stfu = 2*0.21; % Stoichiometric no. of moles of fuel
    fst = (stfu*2)/MW; % Stoichiometric fuel air ratio (molecular weight of H2 = 2)
elseif strcmp(fuel, 'CH4')
    stfu = 0.21/2;
    fst = (stfu*16)/MW; % Stoichiometric fuel air ratio (molecular weight of CH4 = 16)
else
    error('Invalid fuel');
end
f= phi_gg*(1-ra)*fst;

g = 9.81; % acceleration due to gravity

Ms = inlet(1);
Ns = inlet(2);
tes = inlet(3:Ms+2);
Theta_es = transpose(tes);
tis = inlet(Ms+3:Ms+Ns+2);
Theta_is = transpose(tis);
Mr = inlet(Ms+Ns+3);
Nr = inlet(Ms+Ns+4);
ter = inlet(Ms+Ns+5:Ms+Ns+Mr+4);
Theta_er = transpose(ter);
tir = inlet(Ms+Ns+Mr+5:end);
Theta_ir = transpose(tir);

% Area Split:
alpha1r = 1-ra;
alpha1s = ra;
T1r = T0;
T1s = T0;
M1r = M0;
M1s = M0;
P1s = P0;
P1r = P0;
x1s = x_air;
x1r = x_air;


% Compression: (Supersonic stream)
[T4s, P4s, x4s, M4s, alpha4s, beta_es, beta_is] = TPG_MNrampcalculator_p(T1s,P1s,x1s,M1s,Ms,Ns,Theta_es,Theta_is); % From Nik

% External compression: (Subsonic stream)
[T2r, P2r, x2r, M2r, alpha2r, beta_er, beta_ir] = TPG_MNrampcalculator_p(T1r,P1r,x1r,M1r,Mr,Nr,Theta_er,Theta_ir); % From Nik

% Internal Compression: Subsonic part:
[T4r, P4r, x4r, M4r, alpha4r, alpha_throat] = Inlet_sub(T2r, P2r, x2r, M2r);

% Gas Generator Burner:
[T6r, P6r, x6r, M6r, alpha6r] = GG_b(T4r, P4r, x4r, M4r, fuel, stfu, f, phi_gg, ra);

% Gas Generator Nozzle:
P8r = pi_er*P0;
[T8r, x8r, M8r, alpha8r, R8r] = nozzle(T6r, P6r, x6r, M6r, P8r, eta_er);

% Isolator:
if k==0
    P8s = P8r;
    [T8s, x8s, M8s, alpha8s, alpha8sc, R8s] = Iso_P(T4s, P4s, x4s, M4s, P8s);
else
    [T8s, P8s, x8s, M8s, alpha8s, alpha8sc, R8s] = Iso(T4s, P4s, x4s, M4s, k);
end

% Supersonic Combustor:
if ds==1
    [T12, P12, x12, M12, alpha12r] = SC_burner(T8r, P8r, x8r, M8r, R8r, T8s, P8s, x8s, M8s, alpha8sc, R8s, ra, f);
elseif ds==2 && k==0
    [T12, P12, x12, M12, alpha12r] = SCP_burner(T8r, P8r, x8r, M8r, T8s, P8s, x8s, M8s, ra, f);
elseif ds==2 && k~=0
    error('Constant pressure combustion not possible in the supersonic combustor')
end

% Nozzle:
alphae = Alphae/(alpha12r*alpha8r*alpha6r*alpha4r*alpha2r*alpha1r);
[Te, Pe, xe, Me, Re, ae] = Exit_Nozzle(T12, P12, x12, M12, alphae, eta_e);
% Specific Impulse:
% Alphae = alphae*alpha12r*alpha8r*alpha6r*alpha4r*alpha2r*alpha1r;
Pi_e = Pe/P0;
[st,Isp] = SpImp(Me, M0, ae, a0, gamma0, Alphae, Pi_e, f, g);

M = [M1r, M1s, M2r, M4r, M4s, M6r, M8r, M8s, M12, Me];
T = [T1r, T1s, T2r, T4r, T4s, T6r, T8r, T8s, T12, Te];
P = [P1r, P1s, P2r, P4r, P4s, P6r, P8r, P8s, P12, Pe];
x = [x1r, x1s, x2r, x4r, x4s, x6r, x8r, x8s, x12, xe];
x_species = [cellstr(transpose(speciesNames(gas))), num2cell(x1r), num2cell(x1s), num2cell(x2r),...
    num2cell(x4r), num2cell(x4s), num2cell(x6r), num2cell(x8r),...
    num2cell(x8s), num2cell(x12), num2cell(xe)];
alpha = [alpha1r, alpha1s, alpha2r, alpha4r, alpha4s, alpha6r, alpha8r, alpha8s, alpha12r, alphae, alpha_throat];
beta = [beta_es, beta_is, beta_er, beta_ir];

end

