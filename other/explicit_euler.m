function [u] = explicit_euler(G,t)%u,i,j,h,k,a,b,g_x,g_y,m_x,m_y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    u = G.u;

    if(G.parent ~= 0) 
        u = u + G.k.*rhsFineGridAdv(u,t,G);
    else
        u = u + G.k.*rhsAdv(G.u,t,0.5,1,0,G);
    end
    
end

