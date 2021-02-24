function [u] = initialConditionsEulerConstant(x0,y0,t)
% initialConditionsEuler Creates vector U containing solutions for grid G
%   meshgrid for grid G is created by the information about start
%   position, number of nodes and space step for the grid G.

    gamma = 1.4;

%     x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
%     y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';
%     
%     [X,Y] = meshgrid(x,y);
    
    u = ones(length(x0),length(y0),4);
    
    A = meshgrid(x0,y0)';
%     A = randi(10,length(x0),length(y0));
    u(:,:,1) = A.*u(:,:,1); % rho
%     u(:,:,1) = 1*u(:,:,1); % rho
    u(:,:,2) = 2*u(:,:,1).*u(:,:,2); % rho*v1
    u(:,:,3) = 3*u(:,:,1).*u(:,:,3); % rho*v2
    u(:,:,4) = (4/(gamma-1)+(1/2)*(u(:,:,2).^2+u(:,:,3).^2)./u(:,:,1)).*u(:,:,4); % p = 4
    
end

