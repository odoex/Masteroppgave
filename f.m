function [f_u] = f(u)
%f right hand side function 
%   The function f on the right hand side, which is to be differentiated by
%   x, is calculated using the values form u at the given time.
    
    % E = 1/2 rho(v1^2 + v2^2) + p/(gamma-1)
    % E - 1/2 rho(v1^2 + v2^2) = p/(gamma-1) % - eller +?
    
    gamma = 1.81*10^5; % Dobbelsjekk 
    p = (gamma-1)*(u(:,:,4) - (1/2)*(u(:,:,2).^2 + u(:,:,3).^2)./u(:,:,1));

    f_u = u;
    
    f_u(:,:,1) = - u(:,:,2);
    f_u(:,:,2) = - (u(:,:,2).^2)./u(:,:,1) + p; % Dette er feil tror jeg. ro er ikke i andre. Korreksjon: u(:,:,2)^2/u(:,:,1)
    f_u(:,:,3) = - (u(:,:,2).*u(:,:,3))./u(:,:,1);
    f_u(:,:,4) = - (u(:,:,4) + p).*u(:,:,2)./u(:,:,1);% Dont know if E is expressed here or elsewhere, currently expressing it here with last timesteps v1,v2
    
end

