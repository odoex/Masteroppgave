function [f] = flux(i,j,U,a,b,g_x,g_y,m_x,m_y,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    % U(1,j)

    if (i == 1)
        f_x = - a*(U(i+1,j) - g_y(j))/h;
    elseif (i == m_x)
        f_x = - a*(U(i,j) - U(i-1,j))/h;
    else
        f_x = - a*(U(i+1,j) - U(i-1,j))/(2*h);
    end
    
    if (j == 1)
        f_y = - b*(U(i,j+1) - g_x(i))/h;
    elseif (j == m_y)
        f_y = - b*(U(i,j) - U(i,j-1))/h;
    else
        f_y = - b*(U(i,j+1) - U(i,j-1))/(2*h);
    end
    
    f = f_x + f_y;

end

