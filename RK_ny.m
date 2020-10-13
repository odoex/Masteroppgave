function [U] = RK_ny(U,t,h,k,a,b,g_x,g_y,m_x,m_y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    %RK_4(U,t,h,k,a,b,g_x,g_y,m_x,m_y) før jeg sendte in G istede

    
    k1 = rhs_2(U,t,h,a,b,g_x,g_y,m_x,m_y);
    k2 = rhs_2(U + k/2*k1,t + k/2,h,a,b,g_x,g_y,m_x,m_y);
    k3 = rhs_2(U + k/2*k2,t + k/2,h,a,b,g_x,g_y,m_x,m_y);
    k4 = rhs_2(U + k*k3,t + k,h,a,b,g_x,g_y,m_x,m_y);
    
    
    U = U + (k/6)*(k1 + 2*k2 + 2*k3 + k4);
    
%       disp(k1)
%      disp(k2)
%     disp(k3)
%     disp(k4)
%     disp(U)
end
