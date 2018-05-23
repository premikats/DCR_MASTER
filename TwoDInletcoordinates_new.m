function [Xr,Yr,Xc,Yc,Length,Lx] = TwoDInletcoordinates_new(OAR_geo,A0_geo,Theta_e,Theta_i,beta_e,beta_i)
%This function will calculate the coordinates of the vertex points in a 2D
%inlet. It evaluates the coordinates of the vertex  points with Ramp
%leading edge as the origin.
%I/P - M=no. external ramps, N-no of internal ramps, OAR_geo=Overall
%geometric area ratio, Theta_e=matrix of external ramp angles,
%Theta_i=matrix of internal ramp aangles, beta_e=array of external shock angles, beta_i=array of internal shock angles;
% O/P - X=x-coordinates of all the vertex points, Y=y-coordinates of all
% the vertex points,
if OAR_geo~=0
% A0_geo=; % All the inlet coordinates are normalized to stream tube area oof 1m^2.
At_geo=OAR_geo*A0_geo;
theta_r=0;
theta_c=0;
X0=0;
Y0=0;
Y_CLE=A0_geo;
X_CLE=A0_geo/tand(beta_e(1));
RL=0;
CL=0;
%% Calculation co-ordinates of ramp vertices - 
if length(Theta_e)==1
   if length(Theta_i)==1
        i=1;
        Yr=Y_CLE-At_geo;
        Xr=Yr/tand(Theta_e(1));
        Length=Xr;
        RL=sqrt(Xr^2+Yr^2);
   else
        Ax=sqrt(X_CLE^2+Y_CLE^2)*sind(beta_e(1)-Theta_e(1));
        P=1+(tand(Theta_e(1)))^2;
        Q=-2*(X_CLE+(Y_CLE*tand(Theta_e(1))));
        R=X_CLE^2+Y_CLE^2-(Ax/sind(beta_i(1)))^2;
        S=roots([P Q R])
        if S(1)>X_CLE
           i=1;
           Xr=S(1);
           Yr=Xr*tand(Theta_e(1));
        else
           i=1;
           Xr=S(2);
           Yr=Xr*tand(Theta_e(1)); 
        end
        Length=Xr;
        RL=sqrt(Xr^2+Yr^2);
   end
else
for i=1:length(Theta_e)
    theta_r=theta_r+Theta_e(i);
    if i==1
        L(i)=(Y_CLE/sind(beta_e(1)))*sind(beta_e(2)-beta_e(1)+Theta_e(i))/sind(180-beta_e(2));
        Xr(i)=X0+L(i)*cosd(Theta_e(1));
        Yr(i)=Y0+L(i)*sind(Theta_e(1));
        RL=RL+L(i);
    else
        if i==length(Theta_e)
           L(i)=((Y_CLE-Yr(i-1))/sind(beta_e(i)-Theta_e(i)+theta_r))*sind(180-beta_e(i)-beta_i(1)+Theta_e(i))/sind(beta_i(1)); 
        else
          L(i)=((Y_CLE-Yr(i-1))/sind(beta_e(i)-Theta_e(i)+theta_r))*sind(beta_e(i+1)-beta_e(i)+Theta_e(i))/sind(180-beta_e(i+1)); 
        end
        Xr(i)=Xr(i-1)+L(i)*cosd(theta_r);
        Yr(i)=Yr(i-1)+L(i)*sind(theta_r);
        RL=RL+L(i);
    end
end
Length=Xr(i);
RL=RL+0.02;
end
% Calculation of Cowl vertex points co-ordinates - 
X_s=Xr(i);
Y_s=Yr(i);
Xc(1)=X_CLE;
Yc(1)=Y_CLE;
theta_c=0;
for j=1:length(Theta_i)-1
    d(j)=sqrt((Xc(j)-X_s)^2+(Yc(j)-Y_s)^2);
    l(j)=d(j)*sind(beta_i(j+1)-beta_i(j)+Theta_i(j))/sind(180-beta_i(j+1));
    theta_c(1)=theta_c+Theta_i(j);
    Xc(j+1)=Xc(j)+l(j)*cosd(sum(Theta_e)-theta_c);
    Yc(j+1)=Yc(j)+l(j)*sind(sum(Theta_e)-theta_c);
    CL=CL+l(j);
end
Xr=[0 Xr];
Yr=[0 Yr];
CL=CL+(Xr(i+1)-Xc(j+1))+0.02;
Lx=RL+CL;
else
    Xr=zeros(1,length(Theta_e));
    Yr=zeros(1,length(Theta_e));
    Xc=zeros(1,length(Theta_i));
    Yc=zeros(1,length(Theta_i));
    Length=100;
    Lx=100;
end
if Length<0
    Length=100;
    Lx=100;
end
end

