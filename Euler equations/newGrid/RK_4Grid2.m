function [U] = RK_4Grid2(G,t)
% RK_4 function for approximation with Runge-Kutta 4 with stable mesh
% refinement
%   Using Runge-Kutta 4 to calculate the given time step in the method from
%   time t to time t + k. Distinguishes between the root grid and the
%   subgrids which are updated with the function rhsSubgrid instead of rhs.
    
    k1 = rhsGrid2(G.u,t,G);
    k2 = rhsGrid2(G.u + (G.k/2)*k1,t + G.k/2,G);
    k3 = rhsGrid2(G.u + (G.k/2)*k2,t + G.k/2,G);
    k4 = rhsGrid2(G.u + G.k*k3,t + G.k,G);
    
    U = G.u + (G.k/6)*(k1 + 2*k2 + 2*k3 + k4);
end
