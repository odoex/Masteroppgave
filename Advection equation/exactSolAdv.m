function [u] = exactSolAdv(G,t)
%EXACTSOLADV Summary of this function goes here
%   Detailed explanation goes here
    
    a = 0.5;
    b = 1;

    [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
    
    u = sin(X - a*t) + sin(Y - b*t);

end

