function [f_u] = f(u)
%f right hand side function 
%   The function f on the right hand side, which is to be differentiated by
%   x, is calculated using the values form u at the given time.
    
    gamma = 1.4; 
    p = (gamma-1)*(u(:,:,4) - (1/2)*(u(:,:,2).^2 + u(:,:,3).^2)./u(:,:,1));
    
    f_u = u;
    
    f_u(:,:,1) = u(:,:,2);
    f_u(:,:,2) = ((u(:,:,2).^2)./u(:,:,1) + p); 
    f_u(:,:,3) = (u(:,:,2).*u(:,:,3))./u(:,:,1);
    f_u(:,:,4) = (u(:,:,4) + p).*u(:,:,2)./u(:,:,1);
    
end

