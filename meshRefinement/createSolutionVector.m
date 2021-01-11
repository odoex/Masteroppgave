function [U] = createSolutionVector(G)
% createSolutionVector Creates vector U containing solutions for grid G
%   meshgrid for grid G is created by the information about start
%   position, number of nodes and space step for the grid G.

    x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
    y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
    
    [X,Y] = meshgrid(x,y);
    
    F = sin(X) + sin(Y);
    
    U = F;
    
end

