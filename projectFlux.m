function [G] = projectFlux(G_1,G_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    x = G_2.location(1);
    y = G_2.location(2);

    for i = x:(x + (G_2+1)/2)
        for j = y:(y + (G_2+1)/2)
            G_1.u(i,j) = G_2.u(G_1.h/G_2.h*(i-1)+1,G_1.h/G_2.h*(j-1)+1);
        end
    end
    G = G_1;
end

