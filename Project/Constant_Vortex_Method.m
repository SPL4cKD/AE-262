close all;
clear all;
clc;

%% Import Data

fid = fopen('naca2412.dat','rt');
answer = textscan(fid,'%f %f','HeaderLines',1,'Delimiter', ',');
X = answer{1,1}(:,1);
Y = answer{1,2}(:,1);
n = numel(X
figure(1)
plot(X,Y)
axis equal
%%

% for j = 1:n
%     X(j) =
%     Y(j) =
%     Phi(j) =
%
% %Converting panels to clockwise
% for i = 1:n
%     x(i,1) = endpt(n-i+1,1);
%     y(i,2) = endpt(n-i+1,2);
% end
%
% %Establish endpoints
% for i = 1:k
%     pt_1(i,1) = x(i,1);
%     pt_2(i,2) = x(i+1,1);
%     pt_1(i,1) = x(i,2);
%     pt_2(i,2) = x(i+1,2);
% end
%
% %Panel Angles
% for i = 1:k
%     x_th = pt_2(i,1) - pt_1(i,1);
%     z_th = pt_2(i,2) - pt_1(i,2);
%     theta(i) = atan(x_th,z_th);
% end
%%

% for i=1:(numel(x)-1)
%     dydx(i)=(y(i+1)-y(i))/(x(i+1)-x(i));
%     o(i)=atan(dydx(i));
%     Xi(i)=(x(i+1)+x(i))/2;
%     Yi(i)=(y(i+1)+y(i))/2;
%     S(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
%     xi(i) = (x(i+1)+x(i))/2;
%     yi(i) = (y(i+1)+y(i))/2;
%     S(i) = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
%     s(i) = (y(i+1)+y(i))/S(i);
%     c(i) = (x(i+1)+x(i))/S(i);
%     o(i)=atan(s(i) + c(i));
% end
%
% figure(2)
% plot(Xi,Yi)
% axis equal

%%

% alpha = 0;
% n = 159;
%
% for i = 1:n
%     phi(i) = -alpha + atan2((y(i+1)-y(i)),(x(i+1)-x(i)));
% %     beta(i) = phi(i)
%     X(i) = (x(i+1)-x(i))/2;
%     Y(i) = (y(i+1)-y(i))/2;
%     S(i) = sqrt((y+1)-y(i)).^2+(x(i+1)-x(i)).^2;
% end
%%

for i=1:n
    for j=1:n
    %angle of flow with tangent of panel
    phi(i)=-alfa+...
        atan2((Y(i+1)-Y(i)),(X(i+1)-X(i)));
    %angle of flow with normal of panel
    beta(i)=phi(i)+pi/2;
    midpoint_x(i)=(X(i+1)+X(i))/2;
    midpoint_y(i)=(Y(i+1)+Y(i))/2;
    S(i)=sqrt((Y(i+1)-Y(i))^2+...
        (X(i+1)-X(i))^2);%length of panel
    end
end