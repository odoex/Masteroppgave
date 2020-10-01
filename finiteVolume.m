function [G] = finiteVolume(G,a,b,t_0,t_n,f)
% finiteVolume Function running the finite volume scheme.
%   This function takes a grid, a time interval and information about the
%   advection PDE whose solution is to be approximated by the finite volume
%   method. While still inside the time interval the solution is
%   approximated at the current time step on the current grid G (starting
%   at the coarsest grid) with Runge-Kutta 4. The time step for this grid
%   is updated to the time it was last calculated. If this grid has a finer
%   grid, the solution for this grid is recursively calculated with this
%   method. 
    
    t = t_0; 

    while t <= (t_n-G.k)
        
        u = RK_4(G,t,a,b,f); 
        
        if (G.child ~= 0)
            G.child = finiteVolume(G.child,a,b,t,t+G.k,f);
            % Stops at t_n=t+G.k, where G.k is time step at parent grid 
            
            % Her må det skje noe mer hvis G skal ha flere undergrid 
            G.u = u;
            G = projectFlux(G,G.child);
        else
            G.u = u;
        end
        
        t = t + G.k; 
        G.t = t;
        
            

            [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
            %[X1,Y1] = meshgrid(G.child.location(1):G.child.h:G.child.location(1)+G.child.h*(G.child.m-1));

            mesh(X,Y,G.u) 
            hold on
            %mesh(X1,Y1,G.child.u)
            pause
        
        
    end
    
end

