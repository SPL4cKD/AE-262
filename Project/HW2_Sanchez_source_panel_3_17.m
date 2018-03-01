% AE 262                                          Roger Sanchez
% Montgomery                              due February 28, 2017

%                          Homework #2

%% SETUP
clear;
clc;
V_inf	= 10            % Free Stream velocity [m/s]
AoAd    = 0
AoA     = deg2rad(AoAd) % Angle of Attack [radians]

% End/Edge Point Coordinates
n2 = 6
T2 = 2*pi/n2
for i = 1:n2+1
    EP(i,1) =   cos((2*i-3)*T2/2 - pi);
    EP(i,2) =   sin((2*i-3)*T2/2);
end
    
                plot(EP(:,1),EP(:,2))
                    axis([-1 1 -1 1])   % x_min/max y_min/max
                    hold on             % turn on to hold plot
       
n = length(EP)-1;
      
% Critical Point Coordinates
for i = 1:n
    CP(i,1) = [ EP(i+1,1)+ EP(i,1) ] / 2;    % x coordinate
    CP(i,2) = [ EP(i+1,2)+ EP(i,2) ] / 2;    % y coordinate
                   scatter (CP(:,1), CP(:,2))
end

% Beta CP (angle of panel)
for i = 1:n
    bCP(i) = atan2( (EP(i+1,2)-EP(i,2)) , (EP(i+1,1)-EP(i,1)) ) + pi/2; % normal outward (Anderson figure 3.40)
end
    rad2deg(bCP)      % should be 180, -120, -60, 00, 60, 120 for CCW
    

%% CALCULATION

% Affectiveness (coefficients) matrix

h = 5 * 10^2;    % number of pieces

% ############# integration technique (not working correctly)
% for i = 1:n
% %     % Panel affects on Critical Point
% %     for j = 1:n
% %         if j == i
% %           % [An] diagonals should be 1/2 (Anderson, 3.150, pg 251)
% %             An(i,j) = 1/2;  % self induction is for constant source
% %         else            
% %             dx = [EP(j+1,1) - EP(j,1)] / h;
% %             dy = [EP(j+1,2) - EP(j,2)] / h;
% %             ds = sqrt(dx^2 + dy^2);
% %     
% %             % along the panel
% %             for k = 1:h+1
% %               % local points on influencing panel
% %                 x = EP(j,1) + (k-1)*dx;
% %                 y = EP(j,2) + (k-1)*dy;
% %               % Distance to Critical Point (to find affectiveness)
% %                 X = CP(i,1) - x;
% %                 Y = CP(i,2) - y;
% %                 r = sqrt(X^2 + Y^2);
% % % ##########################################################
% %               % Anderson 3.158 and 3.159
% %                     %  d                 (xi - xj)cos(Bi) + (yi - yj)sin(Bi)
% %                     % ---- * (ln r_ij)= -------------------------------------
% %                     % dn_i                    (xi - xj)^2 + (yi - yj)^2
% %                 Int(k) = [X *cos(bCP(i)) + Y *sin(bCP(i))] / (X^2 + Y^2) * ds;
% %             end % k
% %           % Anderson 3.152 and 3.153
% %             I(i,j) = sum(Int);
% %             An(i,j) = 1/(2*pi)*I(i,j)
% % 
% % % ##########################################################
% % 
% %         end % if
% %     end % j
% %     
% %     Tf(i) = bCP(i) - AoA;
% end % i


% ############# linear evaulation of integration (subtitution) technique
% Anderson 3.17 example
            x = CP(:,1);
            y = CP(:,2);
            X = EP(:,1);
            Y = EP(:,2);
            P = bCP - pi/2;     % phi (Anderson 3.160)
for i = 1:n
    % Panel affects on Critical Point
    for j = 1:n
        if j == i
            I(i,j)  = 0;
            I2(i,j) = 0;
        else
            % Anderson 3.162
            A    = -(x(i) - X(j))*cos(P(j)) - (y(i) - Y(j))*sin(P(j));
            B    = (x(i) - X(j))^2 + (y(i) - Y(j))^2;
            C    = sin(P(i) - P(j));
            D    =  (y(i) - Y(j))*cos(P(i)) - (x(i) - X(j))*sin(P(i));
            E    = sqrt(B - A^2);
            S(j) = sqrt([X(j+1) - X(j)]^2 + [Y(j+1) - Y(j)]^2);
            
            % Anderson 3.163
            I(i,j)  = C/2 * log( [S(j)^2 + 2*A*S(j) + B] /B)...
                     + (D - A*C)/E * [ atan( (S(j)+A) /E) - atan(A/E) ];
                 
            % Anderson 3.165
            I2(i,j) = (D - A*C)/(2*E) * log( [S(j)^2 + 2*A*S(j) + B] /B)...
                     - C * [ atan( (S(j)+A) /E) - atan(A/E) ];
            
        end % if
    end % j
end % i

            % Anderson 3.152 and 3.153
            An = 1/(2*pi)*I;
            for i = 1:n
                % [An] diagonals should be 1/2 (Anderson, 3.150, pg 251)
                An(i,i) = 1/2;  % self induction is for constant source
            end
            
    Tf = (bCP - AoA);
            
    Vn = V_inf * cos(Tf);          % positive is out for CW pattern
    Vt = V_inf * sin(Tf);          % positive is CW  for CW pattern
    
% Matrix Calculations        SEE "%----%" for SHORTER syntax
    %       [An][G] + [Vn] =  [Wn]  =   [0]
    %   --> [An][G]        =         - [Vn]
Aug = [         An          ,  -Vn' ];
R   = rref(Aug);
G   = R(:,n+1); 


% ##############################################################
% NEED TO FINISH THIS

V0 = An*G + Vn'    % should be equal to zeros


% Anderson 3.156
%   Vi = Vt.' + 1/(2*pi)*I*G;
% Vs = I * G;     %%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ISSUE HERE!!!!!!!!! d/ds not d/dn (different I)
Vs = 1/(2*pi) * I2 * G;

%------------------------------------------------
                    %-----------------------%
                    R2  = linsolve(An, -Vn');
                    %-----------------------%
                    G2  = R(:,n+1);
                    V02 = An*G2 + Vn.'
                    Vs2 = 1/(2*pi) * I2 * G2;
%------------------------------------------------


for i = 1:n
	Vi(i,1) = Vt(i) + Vs(i);
    
    %------------------------------------------------
                    Vi2(i,1) = Vt(i) + Vs2(i);
    %------------------------------------------------
end

% --------- Seems correct (correct signs)
Rv = Vi/V_inf;   % ratio of velocities
C_P = 1 - Rv.^2


%------------------------------------------------
                    Rv2  = Vi2/V_inf;   % ratio of velocities
                    C_P2 = 1 - Rv2.^2
%------------------------------------------------


% ##############################################################


%% test
  