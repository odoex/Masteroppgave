function [b] = boundaryValuesEuler(x,y,t)
% boundaryValues returns boundary values for the domain
%   Takes the discrete x and y values, and calculates the boundary
%   at the given time for the euler equations

% x and y have to be the same dimension. Can be solved by using the if else
% sentence below and letting v decide the size of b and making a zero
% vector of the same size as v. 

%     if x == 0
%         v = y;
%     else 
%         v = x;
%     end

%     % Constants
%     Ma = 0.1; % Mach number 
%     rho_inf = 1; % Freestream state
%     T_inf = 273.15; % Freestream state
%     R_g = 287.15; % Ideal gas constant
%     gamma = 1.4; % Adiabatic exponent
%     R_v = 0.1; % Radius of vortex 
%     beta = 1; % Vortex strength
%     
%     alpha = 0;
%     c = 331; % Speed of sound
%     v_inf = Ma*c; % Local speed past the boundary: Mach number times the speed of sound
%     c_p = gamma*(R_g/(gamma-1)); % Heat capacity at constant pressure for air, room conditions
% 
%     % Center of vortex
%     x_c = 0.5 + t*v_inf*cos(alpha);
%     y_c = 0.5 + t*v_inf*sin(alpha);
%     
%     % Equations for vortex
%     % Fjern v_inf fra ligningene 
%     f = @(x,y) ((x - x_c).^2 + (y - y_c).^2)./(R_v)^2; 
%     dv_1 = @(x,y) -(beta).*((y-y_c)./R_v).*exp(-f(x,y)./2); % Her regnes med |v_inf|? styrke gange fart? 
%     dv_2 = @(x,y) (beta).*((x-x_c)./R_v).*exp(-f(x,y)./2); % 
%     dT = @(x,y) 0.5*((beta).^2).*(exp(-f(x,y)))./(c_p); % Fjerne fra dT også?
%     
%     b = ones(length(x),4);
%     
%     b(:,1) = b(:,1).*rho_inf.*((T_inf-dT(x,y))./T_inf).^(1./(gamma-1)); % Endrer seg ikke med tid? Hvorfor minus?
%     b(:,2) = b(:,1).*b(:,2).*(v_inf*cos(alpha) + dv_1(x,y));
%     b(:,3) = b(:,1).*b(:,3).*(v_inf*sin(alpha) + dv_2(x,y)); % No velocity in y-direction
%     b(:,4) = b(:,1).*(R_g.*(T_inf-dT(x,y)))/(gamma-1)+(1/2).*(b(:,2).^2+b(:,3).^2)./b(:,1); % P = R*(rho)*T         
  
    u = exactSolEuler(G,t)
    
	if x == 0
        b = u()
    else 
        v = x;
    end

end