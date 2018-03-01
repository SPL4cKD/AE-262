close all;
clear all;
clc;
%% User Selection

choice = menu('Choose a method ','Source Panel','Vortex Panel');
%% Import Data - Source Panel Method

if choice == 1
    
    airfoil = 's1223.dat';
    spaces = ' ';
    header = 1;
    file = importdata(airfoil,spaces,header);
    pts = file.data(:,:);
    %% Calculations
    
    pnl = length(pts)-1;
    disp('Enter a value for AoA:')
    AoA = input('---> ');
    alpha = deg2rad(AoA);
    disp('Enter a value for V_infinity:')
    V_inf = input('---> ');
    %     V_inf = 1;
    
    for i = 1:pnl
        
        %         x1 = x(1:n);
        %         x2 = x(2:n+1);
        %         y1 = y(1:n);
        %         y2 = y(2:n+1);
        
        mid_x(i,1) = (pts(i+1,1) + pts(i,1))/2;
        mid_y(i,2) = (pts(i+1,2) + pts(i,2))/2;
        beta(i) = -(pi/2) + atan2((pts(i+1,2) - pts(i,2)),(pts(i+1,1) - pts(i,1)));
        phi = beta + pi/2;
        
        %     mid_x = (x1+x2)/2;
        %     mid_y = (y1+y2)/2;
        %     S(i) = sqrt((y(i+1)-y(i))^2 + (x(i+1)-x(i))^2);
        
    end
    
    for i = 1:pnl
        for j = 1:pnl
            
            x = mid_x(:,1);
            X = pts(:,1);
            y = mid_y(:,2);
            Y = pts(:,2);
            
            A = -(x(i) - X(j))*cos(phi(j))-(y(i)-Y(j))*sin(phi(j));
            B = (x(i) - X(j))^2 + (y(i)-Y(j))^2;
            C = sin(phi(i)-phi(j));
            D = (y(i)-Y(j))*cos(phi(i)) - (x(i)-X(j))*sin(phi(i));
            E = sqrt(B - A^2);
            S_j = sqrt((X(j+1)-X(j))^2 + (Y(j+1)-Y(j))^2);
            
            if i == j

                I_s(i,j) = 0;
                
            elseif j == pnl
                
                I(i,j) = 0;
                
            else
                
                I(i,j) = -(C/2)*log((S_j^2+2*A*S_j+B)/B) - ((D-A*C)/E)*(atan((S_j+A)/E)-atan(A/E));
                I_s(i,j) = -(D-(A*C))/(2*E)*log((S_j^2+2*A*S_j+B)/B) + C*(atan((S_j+A)/E)-atan(A/E));
                
            end
        end
    end
    
    lambdas = I/(2*pi);
    
    for i = 1:pnl;
        
        lambdas(i,i) = 0.5;
        
    end
    
    V_N = V_inf*cos(beta-alpha);
    V_T = -V_inf*sin(beta-alpha);
    Gamma = linsolve(lambdas,-V_N');
    V_s = (I_s*Gamma)/(2*pi);
    
    for i = 1:pnl;
        
        Vi(i,1) = V_T(i) + V_s(i);
        
    end
    
    Cps = 1-(Vi/V_inf).^2
    
%     for i = 1:pnl/2
%         
%         Cp_u(i) = Cps(i);
%         Cp_l(i) = Cps(i+numel(pnl)/2);
%         
%     end
%     
%     for i = 1:numel(pnl)
%         
%         Cds(i)= Cps(i)*S_j(i)*sin(beta(i));

%         if mid_y(i+1)>=0 && mid_y(i)>=0

%             Cl(i)=-1/c*Cps(i)*S(i)*cos(o(i));

%         elseif Y(i+1)<=0 && Y(i)<=0

%             Cl(i)= 1/c*Cps(i)*S(i)*cos(o(i));

%         end

%         dA(i)=Cl(i)*(X(i+1)-X(i));
%         xdA(i)=xi(i)*Cl(i)*(X(i+1)-X(i));

%     end
%     
%     for i=1:numel(xi)/2
%         Cl_si(i)=Cl(i)+Cl(numel(xi)-i);
%         
%     end
%     
%     Cds
    %% Plots
    
    figure(1)
    
    plot(x,y,'*');
    title('Airfoil Geometry')
    xlabel('x/c')
    ylabel('y/c')
    legend('Location','Best','Control Points')
    axis equal
    
    figure(2)
    
    plot(x(1:(pnl/2)-1),Cps(1:(pnl/2)-1),x(pnl/2:pnl,1),Cps(pnl/2:pnl));
    %     plot(mid_x,Cps);
    title('Pressure Distribution over Airfoil')
    xlabel('x/c')
    ylabel('C_p')
    legend('Location','Best','Upper Half of Airfoil','Bottom Half of Airfoil')
    axis fill
    axis ij
    
    figure(3)
    
    subplot(2,1,1)
    plot(x,y,'*')
    title('Airfoil Geometry')
    xlabel('x/c')
    ylabel('y/c')
    legend('Location','Best','Control Points')
    axis equal
    
    subplot(2,1,2)
    plot(x(1:(pnl/2)-1),Cps(1:(pnl/2)-1),x(pnl/2:pnl,1),Cps(pnl/2:pnl));
    title('Pressure Distribution over Airfoil')
    xlabel('x/c')
    ylabel('C_p')
    legend('Location','Best','Upper Half of Airfoil','Bottom Half of Airfoil')
    axis fill
    axis ij
    %% Import Data - Vortex Panel Method
    
elseif choice == 2
    
    airfoil = 'naca2412.dat';
    spaces = ' ';
    header = 1;
    file = importdata(airfoil,spaces,header);
    pts = file.data(:,:);
    %% Calculations
    
    pnl = length(pts)-1;
    disp('Enter a value for AoA:')
    AoA = input('---> ');
    alpha = deg2rad(AoA);
    disp('Enter a value for V_infinity:')
    V_inf = input('---> ');
    
    for i=1:(pnl(pts)-1)

        mid_x(i)=(pts(i+1) + pts(i))/2;
        mid_y(i)=(pts(i+1) + pts(i))/2;        
        dydx(i)=(pts(i+1) - pts(i))/(pts(i+1)-pts(i));
        phi(i)=atan(dydx(i));        
        S(i)=sqrt((pts(i+1)-pts(i))^2 + (pts(i+1)-pts(i))^2);
        
    end
    
%     for i=1:pnl(mid_x)
%         
%         Cd(i)=-1/c*Cp(i)*S(i)*sin(phi(i));
%         
%         if Y(i+1)>=0 && Y(i)>=0
%             
%             Cl(i)=-1/c*Cp(i)*S(i)*cos(phi(i));
%             
%         elseif Y(i+1)<=0 && Y(i)<=0
%             
%             Cl(i)= 1/c*Cp(i)*S(i)*cos(phi(i));
%             
%         end
%         
%         dA(i)=Cl(i)*(X(i+1)-X(i));
%         xdA(i)=mid_x(i)*Cl(i)*(X(i+1)-X(i));
%         
%     end
%     
%     for i=1:pnl/2
%         
%         Cl_si(i)=Cl(i)+Cl(numel(mid_x)-i);
%         x_si(i)=mid_x(i);
%         y_si(i)=mid_y(i);
%         Cpu(i)=Cp(i);
%         Cpl(i)=Cp(i+numel(mid_x)/2);
%         
%     end
    
    
end