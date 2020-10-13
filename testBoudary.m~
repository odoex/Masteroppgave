function [g_x,g_y] = testBoudary(G,f,t)
% boundart function that calculates the vectors placed at boundaries of
% a grid. 
%   Creates boundaries from the current solution vector from parent. (If
%   parent is updated this could give a mistake as the numbers from which
%   the boundary should be created must be from the same time step as the
%   current one) 

    x = G.location(3);
    y = G.location(4);

    g_x = f(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m),0,t)'; % Hva gjør jeg med tiden her?
    g_y = f(0,linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m),t)'; % Hva gjør jeg med tiden her?
%     if G.parent ~= 0
%         disp(G.parent.u)
% 
%         disp("Boundary: ")
%         disp(g_x)
%         disp(g_y)
%     end
    
    
    
end