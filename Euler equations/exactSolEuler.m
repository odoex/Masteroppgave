function [u] = exactSolEuler(x0,y0,t)
% exactSolEuler Creates vector u containing solutions for grid G
%   G - grid 
%   t - the time for the exact solution
%   the exact solution u is found at time t on the given grid G 

    % Constants
    % RADIUS 
    R_v = 0.1; % Radius of vortex
    % STRENGTH
    beta = 1; % Vortex strength. Higher means more error
    R_g = 287.15; % Ideal gas constant
    gamma = 1.4; % Adiabatic exponent
    c_p = gamma*(R_g/(gamma-1)); % Heat capacity at constant pressure
                                 % for air, room conditions
    % MACH NUMBER
    Ma = 0.0; % Mach number (0.1 when not 0)
    c = 331; % Speed of sound
    
    % ANGLE
    alpha = 0; % angle of free stream velocity
    
    v_inf = Ma*c; % Local speed past the boundary: 
                  % Mach number times the speed of sound
    rho_inf = 1; % Freestream state
    T_inf = 273.15; % Freestream state
    
    % Center of vortex
    x_c = 0.5 + t*v_inf*cos(alpha);
    y_c = 0.5 + t*v_inf*sin(alpha);
    
%     x0 = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
%     y0 = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';
    [X,Y] = meshgrid(x0,y0);
    X = X';
    Y = Y';
    
    % Equations for vortex
    f = @(x,y) ((x - x_c).^2 + (y - y_c).^2)./((R_v).^2); 
    dv_1 = @(x,y) -(beta).*((y-y_c)./R_v).*exp(-f(x,y)./2);
    dv_2 = @(x,y) (beta).*((x-x_c)./R_v).*exp(-f(x,y)./2);
    dT = @(x,y) (0.5).*((beta).^2).*(exp(-f(x,y)))./(c_p);
    
    % Solution
    u = ones(length(x0),length(y0),4);
    
    u(:,:,1) = u(:,:,1).*rho_inf.*((T_inf-dT(X,Y))./T_inf).^(1./(gamma-1)); 
    u(:,:,2) = u(:,:,1).*(v_inf*cos(alpha) + dv_1(X,Y));
    u(:,:,3) = u(:,:,1).*(v_inf*sin(alpha) + dv_2(X,Y));
    u(:,:,4) = u(:,:,4).*((R_g.*u(:,:,1).*(T_inf-dT(X,Y)))./(gamma-1)+(1/2).*(u(:,:,2).^2+u(:,:,3).^2)./u(:,:,1)); 
                 
%     figure
%     quiver(X,Y,u(:,:,2),u(:,:,3));
     
end

