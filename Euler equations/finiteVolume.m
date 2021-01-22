function [G] = finiteVolume(G,t_0,t_n)
% finiteVolume Function running the finite volume scheme.
%   This function takes a grid and a time interval as parametres. While 
%   still inside the time interval the solution is
%   approximated at the current time step on the current grid G (starting
%   at the coarsest grid) with Runge-Kutta 4. The time step for this grid
%   is updated to the time it was last calculated. If this grid has a finer
%   grid, the solution for this grid is recursively calculated with this
%   method.
    
    t = t_0;
    n = round((t_n-t_0)/G.k);
    
    o=0;
     while o < n %abs(t_n-t) > 1e-08
        o=o+1;

        u = RK_4(G,t);
        
        if (G.child ~= 0)
            
            G.child = finiteVolume(G.child,t,t+G.k);
            % Stops at t_n=t+G.k, where G.k is time step at parent grid 
            
            % Her må det skje noe mer hvis G skal ha flere undergrid 
             G.u = u;
             G = projectFlux(G,G.child);
        else
             G.u = u;
        end

        t = t + G.k;
        G.t = t;
        
    end
    
end
