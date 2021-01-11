function [U] = rhsAdv(U,t,a,b,f,G)
% rhs calculates the right side of the equation
%   Collecting the boundary vectors for the grid, the function loops
%   through the grid and approximates the flux. Checking if the point is at
%   the boundary, the function uses the boundary conditions to calculate
%   the flux when at the boundary where boundary condition is needed, and
%   adapts the method at the other boundaries (i = mx or j = my). The
%   solution is inserted to a temporary solution array U_n and is returned
%   in U after the loop is done. 

    h = G.h;
    m_x = G.m_x;
    m_y = G.m_y;
    
    F_x = zeros(m_x);
    F_y = zeros(m_y);
    
    [g_x,g_y] = boundaryAdv(G,t); % Used to be: [g_x,g_y] = boundaryAdv(G,f,t)
    
%     U(:,1)=g_x;
%     U(1,:)=g_y;

    for i = 1:m_x
        for j = 1:m_y
            
            if (i == 1)
                F_x(i,j) = - a*(U(i+1,j) - g_y(j))/h;
            elseif (i == m_x)
                F_x(i,j) = - a*(U(i,j) - U(i-1,j))/h;
            else
                F_x(i,j) = - a*(U(i+1,j) - U(i-1,j))/(2*h);
            end

            if (j == 1)
                F_y(i,j) = - b*(U(i,j+1) - g_x(i))/h;
            elseif (j == m_y)
                F_y(i,j) = - b*(U(i,j) - U(i,j-1))/h;
            else
                F_y(i,j) = - b*(U(i,j+1) - U(i,j-1))/(2*h);
            end
            
        end
    end
    
    U_n = F_x + F_y;
    
    U = U_n;
    
end