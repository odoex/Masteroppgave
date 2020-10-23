function [f_u] = f_adv(u)

    % PDE coefficients and exact solution
    a = 0.5;

%     f = @(x,y,s) sin(x - a*s) + sin(y - b*s);
    
    f_u = - a*u;

end