function [G] = projectFlux(G_1,G_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    x = G_2.location(3);
    y = G_2.location(4);
    r = G_1.h/G_2.h;

    % Starting at 1 until (m-1)/2 + 1 (removing starting point, dividing by
    % the ratio, and adding the starting point) minus 1 since the loop is
    % starting at 0. This is for mathematical reasons.
    for l = 1:length(G_1.u(1,1,:))
        for i = 1:((G_2.m_x - 1)/r -1)
            for j = 1:((G_2.m_y - 1)/r -1)
                % Position in coarse grid is given by x,y and increasing each
                % round with i (first round 0).

                G_1.u(x+i,y+j,l) = G_2.u(r*i+1, r*j+1,l);

            end
        end
    end
    G = G_1;
end

