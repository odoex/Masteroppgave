function [U] = rhs(G,t,a,b,f)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    %rhs(U,t,a,b,f,m_x,m_y,h)

    U = G.u;
    h = G.h;
    m_x = G.m;
    m_y = G.m;
    
    F_x = zeros(m_x);
    F_y = zeros(m_y);
    
    [g_x,g_y] = boundary(G,f,t);

    for i = 1:m_x
        for j = 1:m_y
            
            if (i == 1)
                %F_x(i,j) = - a*(U(i+1,j) - U(i,j))/h;
                F_x(i,j) = - a*(U(i+1,j) - g_y(j))/h;
            elseif (i == m_x)
                F_x(i,j) = - a*(U(i,j) - U(i-1,j))/h;
            else
                F_x(i,j) = - a*(U(i+1,j) - U(i-1,j))/(2*h);
            end

            if (j == 1)
                %F_y(i,j) = - b*(U(i,j+1) - U(i,j))/h;
                F_y(i,j) = - b*(U(i,j+1) - g_x(t,i))/h;
            elseif (j == m_y)
                F_y(i,j) = - b*(U(i,j) - U(i,j-1))/h;
            else
                F_y(i,j) = - b*(U(i,j+1) - U(i,j-1))/(2*h);
            end

            U_n = F_x + F_y;
            
        end
    end
    
    U = U_n;
    
end

