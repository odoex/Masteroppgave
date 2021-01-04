function [u] = initialConditionsEuler(G)
% initialConditionsEuler Creates vector U containing solutions for grid G
%   meshgrid for grid G is created by the information about start
%   position, number of nodes and space step for the grid G.

    % Center of vortex
    x_c = 0.5;
    y_c = 0.5;

    % Constants
    Ma = 0.01; % Mach number (0.1 when not 0) 
    rho_inf = 1; % Freestream state
    T_inf = 273.15; % Freestream state
    R_g = 287.15; % Ideal gas constant
    gamma = 1.4; % Adiabatic exponent
    R_v = 0.1; % Radius of vortex 
    beta = 5; % Vortex strength
    
    c = 331; % Speed of sound
    v_inf = Ma*c; % Local speed past the boundary: Mach number times the speed of sound
    c_p = gamma*(R_g/(gamma-1)); % Heat capacity at constant pressure for air, room conditions

    x0 = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
    y0 = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';
    [Y,X] = meshgrid(x0,y0);
    
    % Equations for vortex
    f = @(x,y) ((x - x_c).^2 + (y - y_c).^2)./((R_v).^2); 
    dv_1 = @(x,y) -(beta).*((y-y_c)./R_v).*exp(-f(x,y)./2);
    dv_2 = @(x,y) (beta).*((x-x_c)./R_v).*exp(-f(x,y)./2);
    dT = @(x,y) (0.5).*((beta).^2).*(exp(-f(x,y)))./(c_p);
    
    u = ones(G.m,G.m,4);
    
    u(:,:,1) = u(:,:,1).*rho_inf.*((T_inf-dT(X,Y))./T_inf).^(1./(gamma-1)); 
    u(:,:,2) = u(:,:,1).*(v_inf + dv_1(X,Y));%u(:,:,2).*u(:,:,1).*(v_inf + dv_1(X,Y)); % rho*v1
    u(:,:,3) = u(:,:,1).*(dv_2(X,Y));%u(:,:,3).*u(:,:,1).*dv_2(X,Y); % rho*v2 -Because there is no initial flow in y-direction
    u(:,:,4) = u(:,:,4).*((R_g.*u(:,:,1).*(T_inf-dT(X,Y)))./(gamma-1)+(1/2).*(u(:,:,2).^2+u(:,:,3).^2)./u(:,:,1)); 
                 
    figure
    quiver(X,Y,u(:,:,2),u(:,:,3)); 
    
end