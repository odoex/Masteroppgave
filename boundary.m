function [g_x,g_y] = boundary(G,f,t)
%UNTITLED2 Summary of this function goes here
%   Creates boundaries from the current solution vector from parent. (If
%   parent is updated this could give a mistake as the numbers from which
%   the boundary should be created must be from the same time step as the
%   current one) 

    if G.parent == 0
        g_x = f(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m),0,t)'; % x er en vektor, t er et tall
        g_y = f(0,linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m),t)';
    else
        ratio = G.parent.h/G.h; % Først for ratio = 2
        g_x = zeros(G.m,1);
        g_y = zeros(G.m,1);
        for i = 1:((G.m+1)/ratio) % looper kun til nest siste punkt i hovedgrid som dekker fint grid
            
            %G.u(i,1) = 
            g_x(i*2-1) = G.parent.u(G.location(3)+(i-1),G.location(4));
            g_y(i*2-1) = G.parent.u(G.location(3),G.location(4)+(i-1));
        end
        i=2;
        while i < G.m
            g_x(i) = (g_x(i-1)+g_x(i+1))/2;
            g_y(i) = (g_y(i-1)+g_y(i+1))/2;
            
            if mod(i,ratio) == 0
                i=i+1;
            end
            i=i+1;
        end
    end
    
end

