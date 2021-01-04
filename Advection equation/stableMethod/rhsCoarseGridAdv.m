function [U] = rhsCoarseGridAdv(U,t,a,b,f,G)

    h = G.h;
    m_x = G.m;
    m_y = G.m;
    
    x = G.child.location(3); % New: might be a better way
    y = G.child.location(4);
    
    % Change the cases for on the boundary and outside the boundary. Also
    % find the best way to combine these two rhs functions in the method.
    % Remember to change the step size if needed. 
    
    F_x = zeros(m_x); % Number of grid points in fine grid pluss an extra point outside the boundaries
    F_y = zeros(m_y);
    
    [g_x,g_y] = fineBoundaryAdv(G,f,t);
    
   
    for i = 1:m_x
        for j = 1:m_y
            
            if (i == 1)
                F_x(i,j) = - a*(U(i+1,j) - g_y(j))/h;
            elseif (i == m_x)
                F_x(i,j) = - a*(U(i,j) - U(i-1,j))/h;
            %elseif (i == x-1 && j > y)
                %F_x(i,j) = - a*(G.child.u(1,2*(j-y)+1) - U(i-1,j))/(2*h);
            else
                F_x(i,j) = - a*(U(i+1,j) - U(i-1,j))/(2*h);
            end

            if (j == 1)
                F_y(i,j) = - b*(U(i,j+1) - g_x(i))/h;
            elseif (j == m_y)
                F_y(i,j) = - b*(U(i,j) - U(i,j-1))/h;
            %elseif (j == y-1 && i > x)
                %F_y(i,j) = - a*(G.child.u(2*(i-x)+1,1) - U(i,j-1))/(2*h);
            else
                F_y(i,j) = - b*(U(i,j+1) - U(i,j-1))/(2*h);
            end
            
        end
    end

    
    
%     for i = 1:5
%         for j = 1:5
%             if (i==1) % first point on coarse part
%                 F_x(2*i-1,2*j-1) = - a*(U(2*i+1,2*j-1) - g_y(2*j-1))/(h);
%             elseif (i==2) % outside fine boundary
%                 if (j==1 || j==2)
%                     F_x(2*i-1,2*j-1) = -a*(U(2*i+1,2*j-1)-U(2*i-3,2*j-1))/(2*h);
%                 elseif (j==3)
%                     F_x(2*i-1,2*j-1) = -a*((U(2*i+1,2*j-1)*(3/4)+U(2*i+1,2*j)*(1/4)) - U(2*i-3,2*j-1))/(2*h);
%                 elseif (j==5)
%                     F_x(2*i-1,2*j-1) = -a*((U(2*i+1,2*j-1)*(1/2)+U(2*i+1,2*j-2)*(1/2)) - U(2*i-3,2*j-1))/(2*h);
%                 else
%                     F_x(2*i-1,2*j-1) = -a*((U(2*i+1,2*j-1)*(1/2)+U(2*i+1,2*j)*(1/4)+U(2*i+1,2*j-2)*(1/4)) - U(2*i-3,2*j-1))/(2*h);
%                 end
%             elseif(j==1 || j==2)
%                 if(i==5)
%                     F_x(2*i-1,2*j-1) = -a*(U(2*i-1,2*j-1)-U(2*i-3,2*j-1))/(h);
%                 else
%                     F_x(2*i-1,2*j-1) = -a*(U(2*i+1,2*j-1)-U(2*i-3,2*j-1))/(2*h);
%                 end
%             end
%             
%             if (j==1) % first point on coarse part
%                 F_y(2*i-1,2*j-1) = - b*(U(2*i-1,2*j+1) - g_x(2*i-1))/(h);
%             elseif (j==2) % outside fine boundary
%                 if (i==1 || i==2)
%                     F_x(2*i-1,2*j-1) = -b*(U(2*i-1,2*j+1)-U(2*i-1,2*j-3))/(2*h);
%                 elseif (i==3)
%                     F_x(2*i-1,2*j-1) = -b*((U(2*i-1,2*j+1)*(3/4)+U(2*i,2*j+1)*(1/4)) - U(2*i-1,2*j-3))/(2*h);
%                 elseif (i==5)
%                     F_x(2*i-1,2*j-1) = -b*((U(2*i-1,2*j+1)*(1/2)+U(2*i-2,2*j+1)*(1/2)) - U(2*i-1,2*j-3))/(2*h);
%                 else
%                     F_x(2*i-1,2*j-1) = -b*((U(2*i-1,2*j+1)*(1/2)+U(2*i,2*j+1)*(1/4)+U(2*i-2,2*j+1)*(1/4)) - U(2*i-1,2*j-3))/(2*h);
%                 end
%             elseif(i==1 || i==2)
%                 if(j==5)
%                     F_x(2*i-1,2*j-1) = -b*(U(2*i-1,2*j-1)-U(2*i-1,2*j-3))/(h);
%                 else
%                     F_x(2*i-1,2*j-1) = -b*(U(2*i-1,2*j+1)-U(2*i-1,2*j-3))/(2*h);
%                 end
%             end
%         end
%     end
%                 
%                 
%                 
%             
%             
%             
%                 % Her skal for alle j behandles som om de er tett intill
%                 % den fine boundaryen 
%                 
%     for i = 5:9
%         for j = 5:9
%             % When the grid starts two points away from the boundary
%             if (i == 5) 
%                 F_x(i,j) = -a*((U(1,1)*(3/4)+U(1,2)*(1/4)) - G.parent.u(x-2,y))/(2*2*h);
%             elseif (i == 2)
%                 
%             elseif (i == m_x)
%                 F_x(i,j) = - a*(U(i,j) - U(i-1,j))/h;
%             else
%                 F_x(i,j) = - a*(U(i+1,j) - U(i-1,j))/(2*h); 
%                 
%             end
% 
%             if (j == 1)
%                 F_y(i,j) = - b*(G.parent.u(x-1,y+1)-G.parent.u(x-1,y-1))/(2*2*h); % 2h = h1
%             elseif (j == 2)
%                 
%             elseif (j == m_y)
%                 F_y(i,j) = - b*(U(i,j) - U(i,j-1))/h;
%             else
%                 F_y(i,j) = - b*(U(i,j+1) - U(i,j-1))/(2*h);
%                 if(j==3)
%                 disp(- b*(U(i,j+1) - U(i,j-1))/(2*h))
%                 end
%             end
%             
%         end
%     end
    
    
    U_n = F_x + F_y;
%     if (G.child ~= 0)
%         disp(U_n)
%     end
    U = U_n;
    

end