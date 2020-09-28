function [U] = RK_4(G,t,a,b,f)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    %RK_4(U,t,h,k,a,b,g_x,g_y,m_x,m_y) før jeg sendte in G istede

    
    U = G.u;
    k=G.k;
    
    k1 = rhs(U,t,a,b,f,G);
    k2 = rhs(U + k/2*k1,t + k/2,a,b,f,G);
    k3 = rhs(U + k/2*k2,t + k/2,a,b,f,G);
    k4 = rhs(U + k*k3,t + k,a,b,f,G);
    
    U = U + (k/6)*(k1 + 2*k2 + 2*k3 + k4);
    
end

