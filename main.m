clear
clc

tic;
Thrust_req = 5000; % Required Thrust
M0 = 4;
alt = 25000;
[P0, T0, rho0, gamma0, R] = Stdatm(alt);
% P0 = 100;
% T0 = 22;
% rho0 = 0.0105;
load x_air;
% gas = GRI30;
% set(gas, 'T', T0, 'P', P0, 'X', x_air);
% gamma0 = cp_mass(gas)/cv_mass(gas);
% R = gasconstant/meanMolecularWeight(gas);
V0 = M0*sqrt(gamma0*R*T0);
q = 0.5*gamma0*P0*M0^2;
if q>1e5
    disp('Structural limit crossed')
elseif q<0.1e5
    disp('Dynamic pressure is too small')
end
ra = 0.4;
m = 1; % m = 1: Mach number ratio should be mentioned; m = 2: Equal pressure at station 8
if m==1
    k = 0.85;
elseif m==2
    k = 0;
else
    disp('Enter a valid input for Isolator')
end
phi = 1;
phi_gg = phi/(1-ra);
fuel = 'H2';
pi_er=9;
ds = 1; % ds = 1: Const area combustion; ds = 2: Const pressure combustion
if ds~=1 && ds~=2
    disp('Enter a valid input for the type of combustion in the supersonic combustor')
end
Alphae = 1.5;

Ms = 2;
Ns = 2;
Theta_es = [8,4];
Theta_is = [6,6];
Mr = 2;
Nr = 2;
Theta_er = [10,2];
Theta_ir = [6,6];

theta4rc = 5;
theta4rd = 5; % (max possible - 7.5) Ref: High Speed Wind Tunnel Design by Alan Pope
theta8rc = 12;

theta = [theta4rc, theta4rd, theta8rc];

inlet = [Ms,Ns,Theta_es,Theta_is,Mr,Nr,Theta_er,Theta_ir];

[st, Isp, M, T, P, x, alpha, beta, x_species,f] = Specific_thrust_TPG(M0, P0, T0, rho0, gamma0, R, x_air, inlet, ra, k, fuel, phi, pi_er, Alphae, ds);

C_T = (2*st)/(V0*Alphae);
mdot0 = Thrust_req/st;
mdotf = f*mdot0;
A0 = mdot0/(rho0*V0);

[xlr,ylr,xur,yur,xls,yls,xus,yus,l,A] = Sizing_TPG(M0,P0,T0,rho0,gamma0,R,A0,inlet,phi,phi_gg,M,T,P,x,alpha,beta,theta);
if xlr(end)>9
    disp('Exceeds maximum length')
elseif yus(end)>0.65
    disp('Exceeds the maximum diameter')
end

row = {'Station','1r','1s','2r','4r','4s','6r','8r','8s','12','e'};
column = {'M';'T';'P';'alpha';'A'};
result(1,:)=cellstr(row);
result(2:length(column)+1,1)=cellstr(column);
result(2,2:end) =num2cell(M);
result(3,2:end) =num2cell(T);
result(4,2:end) =num2cell(P);
result(5,2:end) =num2cell(alpha(1:end-1));
result(6,2:end) =num2cell(A);

plot(xlr,ylr,xur,yur,xls,yls,xus,yus);
ylim([0 1])
% axis([-0.05 0.2 -0.001 0.005])
x_l = xlabel('x(m)');
y_l = ylabel('y(m)');
name1 = 'DCR Geometry for TPG';
name2 = ['M_0 = ' num2str(M0) '; r_a = ' num2str(ra) '; k = ' num2str(k) ...
    '; \phi = ' num2str(phi) '; \pi_{er} = ' num2str(pi_er) '; A_e/A_0 = ' num2str(Alphae)...
    '; Thrust = ' num2str(Thrust_req) 'N'];
name3 = ['Inlet properties => M_s = ' num2str(Ms) '; N_s = ' num2str(Ns) '; \theta_{es} = ' num2str(Theta_es)...
    '; \theta_{is} = ' num2str(Theta_is)];
name4 = ['M_r = ' num2str(Mr) '; N_r = ' num2str(Nr) '; \theta_{er} = '...
    num2str(Theta_er) '; \theta_{ir} = ' num2str(Theta_ir)];
t = title({name1,name2,name3,name4});
filename = ['geoM' num2str(M0)];
time = toc;
% save_fig(gcf, gca, x_l, y_l, t, 'pdf', filename)