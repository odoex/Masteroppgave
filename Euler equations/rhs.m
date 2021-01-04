function [u1] = rhs(u,t,G)
% rhs calculates the right side of the equation
%   Collecting the boundary vectors for the grid, the function loops
%   through the grid and approximates the flux. Checking if the point is at
%   the boundary, the function uses the boundary conditions to calculate
%   the flux when at the boundary where boundary condition is needed, and
%   adapts the method at the other boundaries (i = mx or j = my). The
%   solution is inserted to a temporary solution array U_n and is returned
%   in U after the loop is done. 

% Hvor skal delta bestemmes? 
    h = G.h;
    m_x = G.m_x;
    m_y = G.m_y;
    
    delta = 7; % Avhengig av h? Er delta konstant?
    delta1 = 7;
    
    F_x = zeros(m_x,m_y);
    F_y = zeros(m_x,m_y);
    
    [p_x0,p_xm,p_y0,p_ym] = boundary(G,t);
    
    px = zeros(2,length(p_x0),4);
    py = zeros(2,length(p_y0),4);
    px(1,:,:) = p_x0;
    px(2,:,:) = p_xm;
    py(1,:,:) = p_y0; % Fiks dette rotet
    py(2,:,:) = p_ym; % Fiks dette rotet
    
    f_px = g(px);
    f_py = f(py); % Vektor form i stede for elementvis. Er dette lurt?
    
    U_f = f(u);
    U_g = g(u);
    u1=u;
    
    for l = 1:4
        for i = 1:m_x
            for j = 1:m_y

                if (i == 1)
                    F_x(i,j) = - ((U_f(i+1,j,l) - f_py(1,j,l)) - (delta1/2)*(u(i+1,j,l)-u(i,j,l)))/h;%(-delta1*(u(i+1,j,l)-u(i,j,l)) + delta1*(u(i,j,l))))/h; % Important to have correct sign for boundary value. This is not longer given by this function
                elseif (i == m_x)
                    F_x(i,j) = - ((f_py(2,j,l) - U_f(i-1,j,l)) + (delta1/2)*(u(i,j,l)-u(i-1,j,l)) )/h;%(delta1*(u(i,j,l)) + delta1*(u(i,j,l)-u(i-1,j,l))) )/h;
                else
                    F_x(i,j) = - ((U_f(i+1,j,l) - U_f(i-1,j,l)) - (delta*(u(i+1,j,l)-u(i,j,l)) - delta*(u(i,j,l)-u(i-1,j,l))))/(2*h);
                end

                if (j == 1)
                    F_y(i,j) = - ((U_g(i,j+1,l) - f_px(1,i,l)) - (delta1/2)*(u(i,j+1,l)-u(i,j,l)) )/h;%(delta1*(u(i,j+1,l)-u(i,j,l)) + delta1*u(i,j,l)) )/h;
                elseif (j == m_y)
                    F_y(i,j) = - (f_px(2,i,l) - U_g(i,j-1,l) + (delta1/2)*(u(i,j,l)-u(i,j-1,l)) )/h;%(delta1*u(i,j,l) + delta1*(u(i,j,l)-u(i,j-1,l))) )/h;
                else
                    F_y(i,j) = - ((U_g(i,j+1,l) - U_g(i,j-1,l)) - (delta*(u(i,j+1,l)-u(i,j,l)) - delta*(u(i,j,l)-u(i,j-1,l))))/(2*h);
                end

            end
        end
        u1(:,:,l) = F_x + F_y;
    end

end