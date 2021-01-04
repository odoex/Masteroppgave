function [U] = rhsFineGridAdv2(U,t,a,b,f,G)

    h = G.h;
    m_x = G.m;
    m_y = G.m;
    r = G.parent.h/G.h;
    
    x = G.location(3); % New: might be a better way
    y = G.location(4);
    
    % Change the cases for on the boundary and outside the boundary. Also
    % find the best way to combine these two rhs functions in the method.
    % Remember to change the step size if needed. 
    
    F_x = zeros(m_x); % Number of grid points in fine grid pluss an extra point outside the boundaries
    F_y = zeros(m_y);
    
%     [g_x,g_y] = fineBoundaryAdv(G,f,t);
    

    for i = 1:m_x % 1 + fine grid points
        for j = 1:m_y
            
            if (i == 1 && mod(j,r) == 1) %"Ghost-values" outside the fine grid. 
                
                if (j == 3)
                    F_x(i,j) = -a*((U(i+2,j)*(3/4)+U(i+2,j+1)*(1/4))-G.parent.u(x-1,y+1))/(4*h);
                    F_y(i,j) = -b*(U(i,j+2)-G.parent.u(x,y))/(4*h);
                elseif (j == m_y)
                    F_x(i,j) = -a*((U(i+2,j)*(1/2)+U(i+2,j-1)*(1/2)) - G.parent.u(x-1,y+(j-1)/2))/(4*h); % ??
                    F_y(i,j) = -b*(U(i,j)-U(i,j-2))/(2*h);
                elseif (j > 3)
                    F_x(i,j) = -a*((U(i+2,j)*(1/2)+U(i+2,j+1)*(1/4)+U(i+2,j-1)*(1/4))-G.parent.u(x-1,y+(j-1)/2))/(4*h);
                    F_y(i,j) = -b*(U(i,j+2) - U(i,j-2))/(4*h);
                end
                
            elseif (i == 3 && j >= 3) % Points along the interface in y-direction
                
                if (j == 3)
                    F_x(i,j) = - a*(U(j,i+1) - U(j,i-2))/(3*h);
                elseif (mod(j,r) == 1)
                    F_x(i,j) = - a*(U(i+1,j) - U(i-2,j))/(3*h); % liten h, 2 fordi det går til begge sider.
                else
                    F_x(i,j) = -a*(U(i+1,j) - (U(i-2,j+1)+U(i-2,j-1))*(1/2))/(3*h);
                end
                
            elseif (i == m_x && j >= 3)
                F_x(i,j) = - a*(U(i,j) - U(i-1,j))/h;
            elseif (i > 3 && j >= 3)
                F_x(i,j) = - a*(U(i+1,j) - U(i-1,j))/(2*h); 
            end

            
            if (j == 1 && mod(i,r) == 1) %"Ghost-values" outside the fine grid. 
                
                if (i == 3)
                    F_x(i,j) = -a*((U(i+2,j))-G.parent.u(x,y))/(4*h);
                    F_y(i,j) = -b*((U(i,j+2)*(3/4)+U(i+1,j+2)*(1/4))-G.parent.u(x+1,y-1))/(4*h);
                elseif (i == m_y)
                    F_x(i,j) = -a*(U(i,j)-U(i-2,j))/(2*h);
                    F_y(i,j) = -b*((U(i,j+2)*(1/2)+U(i-1,j+2)*(1/2)) - G.parent.u(x+(i-1)/2,y-1))/(4*h); % ??
                elseif (i > 3)
                    F_x(i,j) = -a*(U(i+2,j) - U(i-2,j))/(4*h);
                    F_y(i,j) = -b*((U(i,j+2)*(1/2)+U(i+1,j+2)*(1/4)+U(i-1,j+2)*(1/4))-G.parent.u(x+(i-1)/2,y-1))/(4*h);
                end
                
            elseif (j == 3 && i >= 3) % points along the interface in x-direction
                
                if (i == 3)
                    F_y(i,j) = - b*(U(j+1,i) - U(j-2,i))/(3*h);
                elseif (mod(i,r) == 1)
                    F_y(i,j) = - b*(U(i,j+1) - U(i,j-2))/(3*h); % liten h, 2 fordi det går til begge sider.
                else
                    F_y(i,j) = - b*(U(i,j+1) - (U(i-1,j-2)*(1/2)+U(i+1,j-2)*(1/2)))/(3*h);
                end
                
            elseif (j == m_y && i >= 3)
                F_y(i,j) = - b*(U(i,j) - U(i,j-1))/h;
            elseif (j > 3 && i >= 3)
                F_y(i,j) = - b*(U(i,j+1) - U(i,j-1))/(2*h);
            end
            
        end
    end
    
    
    U_n = F_x + F_y;
%     if (G.child ~= 0)
%         disp(U_n)
%     end
    U = U_n;
    

end