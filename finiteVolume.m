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

    while floor(t*100000)/100000 < floor(t_n*100000)/100000
        
        u = RK_4(G,t,a,b,f); 
        
        %
        G_n = Node(-1,G.location,G.h,G.k,G.n,G.m);
        G_n.u = u;
        
        %
        
        if (G.child ~= 0)
            %
            G.child.sibling = G_n;
            %
            G.child = finiteVolume(G.child,a,b,t,t+G.k,f);
            % Stops at t_n=t+G.k, where G.k is time step at parent grid 
            
            % Her må det skje noe mer hvis G skal ha flere undergrid 
             G.u = u;
             G = projectFlux(G,G.child);
        else
             G.u = u;
        end
        
%         G.u = RK_4(G,t,a,b,f); 
        
        t = t + G.k; 
        G.t = t;
        
            
        %if G.parent ~= 0
            [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
            %[X1,Y1] = meshgrid(G.child.location(1):G.child.h:G.child.location(1)+G.child.h*(G.child.m-1));
            E = abs((sin(X-a*t) + sin(Y-b*t))'-G.u);
            disp(t)
            disp(E)
            
%             if G.parent == 0
%                 mesh(X,Y,(sin(X-a*t) + sin(Y-b*t))')
%             end

%             mesh(X,Y,G.u)
%             hold on
%             pause
            
            
%             mesh(X,Y,E)
%             hold on 
%             pause
            
        %end
        
    end
    
end

