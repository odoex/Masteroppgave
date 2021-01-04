function [U] = createSolutionVector(G)
% createSolutionVector Creates vector U containing solutions for grid G
%   meshgrid for grid G is created by the information about start
%   position, number of nodes and space step for the grid G.

    x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
    y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';
    
    [X,Y] = meshgrid(x,y);
    
    F = sin(X) + sin(Y);
    
    U = F;
    
end

