function [G] = projectFlux2(G_1,G_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    x = G_2.location(3);
    y = G_2.location(4);
    r = G_1.h/G_2.h;
    
    U = zeros(G_2.m_x/r,G_2.m_y/r,4);
    
    for i = 1:r
        for j = 1:r
            U = U + G_2.u(i:r:end-r+i,j:r:end-r+j,:);
        end
    end
    
    G_1.u(x:(x+G_2.m_x/r-1),y:(y+G_2.m_y/r-1),:) = (1/(r^2))*U;

    %G_1.u(x:(x+G_2.m_x/r-1),y) = 
    %G_1.u() = 
    
    G = G_1;
end

