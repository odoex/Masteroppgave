function [G] = finiteVolumeGrid2(G,t_0,t_n)
% finiteVolume Function running the finite volume scheme.
%   This function takes a grid and a time interval as parametres. While 
%   still inside the time interval the solution is
%   approximated at the current time step on the current grid G (starting
%   at the coarsest grid) with Runge-Kutta 4. The time step for this grid
%   is updated to the time it was last calculated. If this grid has a finer
%   grid, the solution for this grid is recursively calculated with this
%   method.
    j = 100;
    t = t_0;
    n = round((t_n-t_0)/G.k);
    figure
	for i = 1:n

        u = RK_4Grid2(G,t);
        
        if (G.child ~= 0)
            
            G.child = finiteVolumeGrid2(G.child,t,t+G.k);
            % Stops at t_n=t+G.k, where G.k is time step at parent grid 
            
            G.u = u;
            G = projectFlux(G,G.child);
        else
            G.u = u;
        end
        
%         if (i > j)
%             j = j + 100;
%             x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
%             y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
%             [X,Y] = meshgrid(x,y);
%             X = X'; Y = Y';
% 
%             mesh(X,Y,G.u(:,:,1))
%             pause
%         end

        t = t + G.k;
        G.t = t;
        
	end %sette på pause og sjekke at alle funksjonene som skal bli kalt blir kalt 
    
end