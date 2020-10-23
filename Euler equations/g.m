function [g_u] = g(u,gamma)
%f right hand side function 
%   The function g on the right hand side, which is to be differentiated by
%   y, is calculated using the values form u at the given time.
    
    %gamma = 1.81*10^5;
    p = (gamma-1)*(u(:,:,4) - (1/2)*(u(:,:,2)^2 + u(:,:,3)^2)/u(:,:,1));

    g_u(1) = u(:,:,3);
    g_u(2) = (u(:,:,2)*u(:,:,3))/u(:,:,1);
    g_u(3) = (u(:,:,3)^2)/u(:,:,1) + p;
    g_u(4) = (u(:,:,4) + p)*u(:,:,3)/u(:,:,1);
    
end % Dobbeltsjekk at denne funksjonen er riktig 