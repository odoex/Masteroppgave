function [G] = finiteVolume(G,t_0,t_n) % Sjekk om det kun er parametrene som er anderledes
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

     while abs(t_n-t) > 1e-08
        
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
        
        [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1)); 
        mesh(X,Y,G.u(:,:,1))
        mesh(X,Y,G.u(:,:,2))
        mesh(X,Y,G.u(:,:,3))
        mesh(X,Y,G.u(:,:,4))
% hold on
        
%         disp(G.u(:,:,1))
%         disp(G.u(:,:,2))
%         disp(G.u(:,:,3))
%         disp(G.u(:,:,4))

        t = t + G.k;
        G.t = t;
        
    end
    
end
