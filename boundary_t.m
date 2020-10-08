function [g_x,g_y] = boundary_t(G,f,t)
% boundart function that calculates the vectors placed at boundaries of
% a grid. 
%   Creates boundaries from the current solution vector from parent. (If
%   parent is updated this could give a mistake as the numbers from which
%   the boundary should be created must be from the same time step as the
%   current one) 

    x = G.location(3);
    y = G.location(4);

    if G.parent == 0
        g_x = f(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m),0,t)'; % Hva gjør jeg med tiden her?
        g_y = f(0,linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m),t)'; % Hva gjør jeg med tiden her?
    else
        
        r = G.parent.h/G.h; 
        g_x = zeros(G.m,1);
        g_y = zeros(G.m,1);
        
        %U = [G.parent.u,G.sibling.u];
        
        for i = 0:((G.m-1)/r)
            
            g_x(i*r+1) = G.parent.u(x+i,y);
            g_y(i*r+1) = G.parent.u(x,y+i);
            
        end
        
        i=2;
        gx=[g_x(1),g_x(1+r)];
        gy=[g_y(1),g_y(1+r)];
        
        while i < G.m
            s = mod(i-1,r);
            
            g_x(i) = ((r-s)*gx(1) + s*gx(2))/r;
            g_y(i) = ((r-s)*gy(1) + s*gy(2))/r;
            
            if (s == r-1) && (i+r+1 <= G.m)
                gx = [gx(2),g_x(i+r+1)];
                gy = [gy(2),g_y(i+r+1)];
                i=i+1;
            end
            i=i+1;
        end
        %
        %
        %g_x=-g_x;
        %g_y=-g_y;
        %
        %
    end
%     if G.parent ~= 0
%         disp(G.parent.u)
% 
%         disp("Boundary: ")
%         disp(g_x)
%         disp(g_y)
%     end
    
    
    
end

