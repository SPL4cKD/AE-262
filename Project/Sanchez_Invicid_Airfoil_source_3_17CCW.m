% AE 262                                                      Roger Sanchez
% Montgomery                                                due May 2, 2017

%                        Invicid Airfoil Program

%% SETUP

clear;  clc; close all;
V_inf	= 10                % Free Stream velocity [m/s]
AoAd    = 0                 % Angle of Attack [degrees]
AoA     = deg2rad(AoAd);    % Angle of Attack [radians]

filename1 = 'naca2412.dat';
filename2 =    's1223.dat';

delimiter = ' ';        % character that specifies separation of values
headerlines = 1;        % number of headerlines (to exclude from data)
DAT = importdata(filename1,delimiter,headerlines)
EP  = DAT.data(:, :)
                figure
                subplot (3, 1, 1)           % 2x1 plot organization
                plot(EP(:,1),EP(:,2))
                    title('NACA 2412 Airfoil Panels')
                    pbaspect([4 1 1])       % attempts a 4x1 aspect ratio
                    hold on                 % turn on to hold plot
n = length(EP)-1;
      
% Critical Point Coordinates
for i = 1:n
    CP(i,1) = [ EP(i+1,1)+ EP(i,1) ] / 2;   % x coordinate
    CP(i,2) = [ EP(i+1,2)+ EP(i,2) ] / 2;   % y coordinate
                scatter (CP(:,1), CP(:,2))
end
                    hold off                % not needed?
                
% Beta CP (angle of panel)
for i = 1:n
    % normal outward (Anderson figure 3.40)                   (mod for CCW)
    bCP(i) = atan2( (EP(i+1,2)-EP(i,2)) , (EP(i+1,1)-EP(i,1)) ) - pi/2; 
end
  
%% CALCULATION

% ####### Linear evaulation of integration (Subtitution) technique ####### 
% Anderson 3.17 example
            x = CP(:,1);
            y = CP(:,2);
            X = EP(:,1);
            Y = EP(:,2);
            P = bCP + pi/2;     % phi (Anderson 3.160)        (mod for CCW)
            
% Affectiveness (coefficients) matrix       
for i = 1:n
    % Panel affects on Critical Point
    for j = 1:n
        if j == i
            In(i,j)  = 0;
            Is(i,j)  = 0;
        else
            % Anderson 3.162
            A    = -(x(i) - X(j))*cos(P(j)) - (y(i) - Y(j))*sin(P(j));
            B    =  (x(i) - X(j))^2 + (y(i) - Y(j))^2;
            C    =   sin(P(i) - P(j));
            D    =  (y(i) - Y(j))*cos(P(i)) - (x(i) - X(j))*sin(P(i));
            E    =   sqrt(B - A^2);
            S(j) =   sqrt([X(j+1) - X(j)]^2 + [Y(j+1) - Y(j)]^2);
            
            % Anderson 3.163                                  (mod for CCW)
            In(i,j)  = -C/2 * log( [S(j)^2 + 2*A*S(j) + B] /B)...
                       - (D - A*C)/E * [ atan( (S(j)+A) /E) - atan(A/E) ];
                 
            % Anderson 3.165                                  (mod for CCW)
            Is(i,j) = -(D -A*C)/(2*E) * log( [S(j)^2 + 2*A*S(j) + B] /B)...
                      + C * [ atan( (S(j)+A) /E) - atan(A/E) ];
        end % if
    end % j
end % i

% Anderson 3.152 and 3.153
An = 1/(2*pi)*In;
for i = 1:n
    % [An] diagonals should be 1/2 (Anderson, 3.150, pg 251)
    An(i,i) = 1/2;              % self induction is for constant source
end

            
Tf = (bCP - AoA);
            
Vn =  V_inf * cos(Tf);          % positive is out for CW pattern
Vt = -V_inf * sin(Tf);          % positive is CW  for CW pattern
    
% Kutta Jakofski

% Matrix Calculations
    %       [An][G] + [Vn] =  [Wn]  =   [0]
    %   --> [An][G]        =         - [Vn]
G = linsolve(An, -Vn');

V0 = An*G + Vn';                % should be equal to zeros

% Anderson 3.156
Vs = 1/(2*pi) * Is * G;         % Panel-induced velocities (summed)

for i = 1:n
	Vi(i,1) = Vt(i) + Vs(i);    % Total induced tangential velocities
end

Rv = Vi/V_inf;                  % ratio of velocities
C_P = 1 - Rv.^2
                subplot (3, 1, 2)
                plot(CP(:,1), C_P)
                    title('NACA 2412 Pressure Distribution')
                    xlabel('distance along chordline')
                    ylabel('Coefficient of Pressure (C_P)')
                    set(gca,'YDir', 'reverse')
                    subplot (3, 1, 3)
                plot(CP(1:(n/2)-1), C_P(1:(n/2)-1))
                    title('NACA 2412 Pressure Distribution')
                    xlabel('distance along chordline')
                    ylabel('Coefficient of Pressure (C_P)')
                    set(gca,'YDir', 'reverse') 
  