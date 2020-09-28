function [u_ij] = ex_euler(u,i,j,h,k,a,b,g_x,g_y,m_x,m_y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    U = @(x,y) u(x,y);

    u_ij = u(i,j) + k*flux(i,j,U,a,b,g_x,g_y,m_x,m_y,h);

end

