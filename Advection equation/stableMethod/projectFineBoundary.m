function [G] = projectFineBoundary(G_1,G_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    x = G_2.location(3);
    y = G_2.location(4);
    r = G_1.h/G_2.h;

    % Starting at 1 until (m-1)/2 + 1 (removing starting point, dividing by
    % the ratio, and adding the starting point) minus 1 since the loop is
    % starting at 0. This is for mathematical reasons.
    for i = 0:((G_2.m_x - 1)/r)
        for j = 0:((G_2.m_y - 1)/r)
            % Position in coarse grid is given by x,y and increasing each
            % round with i (first round 1). 
            
%             a=1;
%             b=0.5;
%             t=G_1.t+G_1.k;
%             [X,Y] = meshgrid(G_1.location(1):G_1.h:G_1.location(1)+G_1.h*(G_1.m-1));
%             sol=(sin(X-a*t) + sin(Y-b*t))';
            
            if ~(j==0 && i==0)
                G_1.u(x+i,y+j) = G_2.u(r*i+1, r*j+1);
            end
            
            
            
        end
        % For ekstra funksjon g
%         G_1.u(x-1,y+i) = G_2.g(2*i+1,2); % Fiks for større refinement
%         G_1.u(x+i,y-1) = G_2.g(2*i+1,1);
    end
    G = G_1;
end