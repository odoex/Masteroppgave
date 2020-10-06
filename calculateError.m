function [error] = calculateError(G,a,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
    y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';
    
    [X,Y] = meshgrid(x,y);
    
    sol = sin(X - a*G.t) + sin(Y - b*G.t);
    
    E = abs(G.u-sol');
    error = norm(E);

    
end
