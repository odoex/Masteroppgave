function [U] = timeStep(U,h,k,a,b,g_x,g_y,m_x,m_y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    for i = 1:m_x
        for j = 1:m_y
            U(i,j) = RK(U,i,j,h,k,a,b,g_x(t,x),g_y(t,y),m_x,m_y);
        end 
    end
    
end

