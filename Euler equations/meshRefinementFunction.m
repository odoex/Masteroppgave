function [G] = meshRefinementFunction(G,t_0,t_n)
% meshRefinementFunction Function running the finite volume scheme.
%   This function takes a grid and a time interval as parametres. While 
%   still inside the time interval the solution is
%   approximated at the current time step on the current grid G (starting
%   at the coarsest grid) with Runge-Kutta 4. The time step for this grid
%   is updated to the time it was last calculated. If this grid has a finer
%   grid, the solution for this grid is recursively calculated with this
%   method. The method currently only works with one subgrid.
    
    t = t_0;
    n = round((t_n-t_0)/G.k);
    
	for i = 1:n

        u = RK_4(G,t);
        
        if (G.child ~= 0)
            
            % Marching in time from t to t + G.k, G.k being time step for 
            % parent grid
            G.child = meshRefinementFunction(G.child,t,t+G.k);
            
            G.u = u;
            G = projectFlux(G,G.child);
        else
            G.u = u;
        end

        t = t + G.k;
        G.t = t;
        
	end
    
end
