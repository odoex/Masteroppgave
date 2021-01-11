function [u] = exactSolAdv(x0,y0,t)
%EXACTSOLADV Summary of this function goes here
%   Detailed explanation goes here
    
    a = 0.5;
    b = 1;

    [X,Y] = meshgrid(x0,y0);
    X = X';
    Y = Y';
    
%     u = zeros(length(x0),length(y0),1);
%     u(:,:,1) = sin(X - a*t) + sin(Y - b*t);
    
    u = sin(X - a*t) + sin(Y - b*t);

end

