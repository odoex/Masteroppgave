function [U] = RK_4(G,t,a,b,f)
% RK_4 function for approximation with Runge-Kutta 4
%   Using Runge-Kutta 4 to calculate the given time step in the method from
%   time t to time t + k. 
    
    k1 = rhs(G.u,t,a,b,f,G);
    k2 = rhs(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
    k3 = rhs(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
    k4 = rhs(G.u + G.k*k3,t + G.k,a,b,f,G);
    
    [g_x,g_y] = boundary_t(G,f,t);
    
    U = G.u + (G.k/6)*(k1 + 2*k2 + 2*k3 + k4);

%     U(:,1) = g_x;
%     U(1,:) = g_y;
    
    % Hvordan setter man inn randbetingelser direkte/uten å regne ut flux
    
end

