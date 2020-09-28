function [G] = projectFlux(G_1,G_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    x = G_2.location(3);
    y = G_2.location(4);
    r = G_1.h/G_2.h;

    for i = x:(G_2.m+1)/2
        for j = y:(G_2.m+1)/2
            l1 =r*(i-x+1)-1;
            l2 = r*(j-y+1)-1;
            G_1.u(i,j) = G_2.u(r*(i-x+1)-1, r*(j-y+1)-1);
        end
    end
    G = G_1;
end

