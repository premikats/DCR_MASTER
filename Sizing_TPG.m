function [xlr,ylr,xur,yur,xls,yls,xus,yus,l,A] = Sizing_TPG(M0,P0,T0,rho0,gamma0,R,A0,inlet,phi,phi_gg,M,T,P,x,alpha,beta,theta)
%Sizing_TPG
%   Length and Dimensions of DCR under TPG conditions

% [P0, T0, rho0, gamma0, R] = Stdatm(alt);
V0 = M0.*sqrt(gamma0*R*T0);
eta_er=0.9;

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

theta4rc = theta(1);
theta4rd = theta(2);
theta8rc = theta(3);

beta_es = beta(1:Ms);
beta_is = beta(Ms+1:Ms+Ns);
beta_er = beta(Ms+Ns+1:Ms+Ns+Mr);
beta_ir = beta(Ms+Ns+Mr+1:end);

A1r = alpha(1)*A0;
A1s = alpha(2)*A0;

% Supersonic Compression:
A4s = alpha(5)*A1s;
[xrs,yrs,xcs,ycs,l4s] = TwoDInletcoordinates_new(alpha(5), A1s, Theta_es,Theta_is,beta_es,beta_is); % From Nik

%Subsonic Compression:
%External Compression:
A2r = alpha(3)*A1r;
[xrr,yrr,xcr,ycr,l1r] = TwoDInletcoordinates_new(alpha(3), A1r, Theta_er,Theta_ir,beta_er,beta_ir); % From Nik
%Internal Compression:
A4r = alpha(4)*A2r;
At_comp = alpha(11)*A2r;
[l4rc, l4rd] = length_subinlet_TPG(A2r, alpha(4), alpha(11), theta4rc, theta4rd);
%% Gas Generator:

A6r = alpha(6)*A4r;
l6r = length_GGburner_TPG(P(4), T(4), A4r, phi_gg, V0);

A8r = alpha(7)*A6r;
[l8rc, l8rd, theta8rd, throat_alpha8r] = length_GGnozzle_TPG(T(6), P(6), x(:,6), M(6), A6r, T(7), P(7), x(:,7), M(7), A8r, theta8rc, eta_er);
At_nozzle = A6r*throat_alpha8r;

%% Isolator:

A8s = alpha(8)*A4s;
l8s = length_isolator_TPG(T(5), P(5), x(:,5), M(5), A4s, P(8), phi);

%% Supersonic Combustor:

A12 = alpha(9)*A8r;
l12 = length_supburner_TPG(P(7), P(8), T(7), T(8), A12, phi, V0);
%% Nozzle:

Ae = alpha(10)*A12;
le = length_nozzle_TPG(T(end-1), P(end-1), x(:,end-1), M(end-1), A12, T(end), P(end), x(:,end), M(end), alpha(10));

%% Areas:
A = [A1r,A1s,A2r,A4r,A4s,A6r,A8r,A8s,A12,Ae];

%% Length:

l_mix_nozz = l12+le;
l_s = l4s+l8s;
l_r = l1r+l4rc+l4rd+l6r+l8rc+l8rd;
l_tot_s = l_s+l_mix_nozz;
l_tot_r = l_r+l_mix_nozz;
l = [l_tot_s;
          l_tot_r];
%% XY Coordinates:

% Subsonic part:
% lower coordinates:
x4rc = xrr(end)+l4rc;
x4rd = x4rc+l4rd;
x6r = x4rd+l6r;
x8rc = x6r+l8rc;
x8rd = x8rc+l8rd;
x12 = x8rd+l12;
xe = x12+le;
xlr = [xrr,x4rc,x4rd,x6r,x8rc,x8rd,x12,xe];
ylr = zeros(1,length(xlr));
ylr(1:length(yrr))=yrr;
ylr(length(yrr)+1:end)=yrr(end);
% upper coordinates:
u = length(xlr)-length(xrr)+length(xcr)+1-2;
xur = zeros(1,u);
xur(1:length(xcr)) = xcr;
xur(length(xcr)+1) = xrr(end);
xur(length(xcr)+2:end) = xlr(length(xrr)+1:length(xlr)-2);
yur = [ycr,ycr(end),yrr(end)+At_comp,yrr(end)+A4r,yrr(end)+A6r,yrr(end)+At_nozzle,yrr(end)+A8r];

% Supersonic part:
% lower coordinates:
x4s = x8rd-l8s;
x1s = x4s-l4s;
xxx = x1s*ones(1,length(xrs));
xrs = xrs+xxx;
xls = [xrs,x8rd];
yyy = yrs(end)-yrs(1);
y1s = yur(end)-yyy;
yy = y1s*ones(1,length(yrs));
yrs = yrs+yy;
yls = [yrs,yur(end)];
% upper coordinates:
xx=x1s*ones(1,length(xcs));
xcs = xcs+xx;
xus = [xcs,xrs(end),x8rd,x12,xe];
y = y1s*ones(1,length(ycs));
ycs = ycs+y;
yu = yur(end)+A8s;
yus = [ycs,ycs(end),yu,yrr(end)+A12,yrr(end)+Ae];


end

