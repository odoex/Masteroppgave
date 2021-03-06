function [G] = finiteVolumeStableMethod(G,t_0,t_n) % Sjekk om det kun er parametrene som er anderledes
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
    n = round((t_n-t_0)/G.k);
    
     for i = 1:n

        u = RK_4StableMethod(G,t);
        
        if (G.child ~= 0)
            
            G.child = finiteVolumeStableMethod(G.child,t,t+G.k);
            % Stops at t_n=t+G.k, where G.k is time step at parent grid 
            
             G.u = u;
             G = projectFluxStableMethod(G,G.child);
        else
             G.u = u;
        end

%         if (G.child ~= 0)
%             x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
%             y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
%             [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m_x-1));
%             X = X';
%             Y = Y';
% %             figure
%             mesh(X,Y,G.u(:,:,1))
%             pause;
%         end
        
        t = t + G.k;
        G.t = t;
    end
    
end

