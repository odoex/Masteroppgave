function [u] = explicit_euler(G,t)%u,i,j,h,k,a,b,g_x,g_y,m_x,m_y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    u = G.u;

    u = u + G.k.*rhs(u,t,G);

end

