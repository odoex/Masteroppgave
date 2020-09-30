function [g_x,g_y] = boundary(G,f,t)
% boundart function that calculates the vectors placed at boundaries of
% a grid. 
%   Creates boundaries from the current solution vector from parent. (If
%   parent is updated this could give a mistake as the numbers from which
%   the boundary should be created must be from the same time step as the
%   current one) 

    x = G.location(3);
    y = G.location(4);

    if G.parent == 0
        g_x = f(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m),0,t)'; 
        g_y = f(0,linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m),t)';
    else
        
        r = G.parent.h/G.h; 
        g_x = zeros(G.m,1);
        g_y = zeros(G.m,1);
        
        for i = 0:((G.m-1)/r)
            
            g_x(i*r+1) = G.parent.u(x+i,y);
            g_y(i*2-1) = G.parent.u(x,y+i);
            
        end
        
        i=2;
        
        while i < G.m
            
            g_x(i) = (g_x(i-1)+g_x(i+1))/2;
            g_y(i) = (g_y(i-1)+g_y(i+1))/2;
            
            if mod(i,r) == 0
                i=i+1;
            end
            i=i+1;
        end
    end
    
end

