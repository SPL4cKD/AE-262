close  all;
clear all;
clc;

% load( 'NACA2210_XY_121pts.mat' )
L=221/12; %[ft]
rho=15e-4; %[slug/ft]
V=526; %[ft/s]
Aref=242.1; %[ft^2]
q=1/2*rho*V^2;
c=1;
X=X*c;
Y=Y*c;

for i=1:(numel(X)-1)
    S(i)=sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2);
    dydx(i)=(Y(i+1)-Y(i))/(X(i+1)-X(i));
    o(i)=atan(dydx(i));
    xi(i)=(X(i+1)+X(i))/2;
    yi(i)=(Y(i+1)+Y(i))/2;
end

for i=1:numel(xi)
    Cd(i)=-1/c*Cp(i)*S(i)*sin(o(i));
    if Y(i+1)>=0 & Y(i)>=0
        Cl(i)=-1/c*Cp(i)*S(i)*cos(o(i));
    elseif Y(i+1)<=0 & Y(i)<=0
        Cl(i)= 1/c*Cp(i)*S(i)*cos(o(i));
    end
    dA(i)=Cl(i)*(X(i+1)-X(i));
    xdA(i)=xi(i)*Cl(i)*(X(i+1)-X(i));
end

for i=1:numel(xi)/2
    Cl_si(i)=Cl(i)+Cl(numel(xi)-i);
    x_si(i)=xi(i);
    y_si(i)=yi(i);
    Cpu(i)=Cp(i);
    Cpl(i)=Cp(i+numel(xi)/2);
    Thickness(i)=abs(Y(i)-Y(numel(xi)-i));
end

CL=sum(Cl);
CD=sum(Cd);
A=sum(dA);
PC=sum(xdA)/A/c;
figure

% Cl_si=fliplr(Cl_si);
yyaxis  left
% plot(x_si,Cl_si*q)

plot(x_si,[Cpu',flip(Cpl')])
set(gca, 'Ydir' , 'reverse' )
ylabel( 'Pressure Coefficient' )
yyaxis  right

plot(X,Y)
ylim([-.4,.4])
ylabel( 'NACA 2210' )
xlim([0,1])
legend( 'Upper Surface' , 'Lower Surface' )
xlabel( ' % chord' )
title( 'Surface Pressure - Vortex Panel Method' )

xPC=PC;
yPC=0;

figure
plot(x_si,Cl_si)
ylabel( 'dC_l' )
yyaxis  right
plot(X,Y,xPC,yPC, 'ok' )
ylim([-.238,.75])
ylabel( 'NACA 2210' )
xlim([0,1])
xlabel( ' % chord' )
title( 'dC_l - Vortex Panel Method' )