function [U] = rhsFineGridAdv(U,t,G)
% Only works for a refinement of r=2

    h = G.h;
    m_x = G.m;
    m_y = G.m;
    r = G.parent.h/G.h; % So much to do
    
    x = G.location(3); 
    y = G.location(4);
    xm = G.location(3) + (m_x-1)/r; 
    ym = G.location(4) + (m_y-1)/r;
    
    % Change the cases for on the boundary and outside the boundary. Also
    % find the best way to combine these two rhs functions in the method.
    % Remember to change the step size if needed. 
    
    F_x = zeros(m_x); % Number of grid points in fine grid pluss an extra point outside the boundaries
    F_y = zeros(m_y); 
    
%     [g_x,g_y] = fineBoundaryAdv(G,f,t);
    
    U_f = f(U);
    U_g = g(U);
    
    %...
    
    % Creating two vectors b_x and b_y for boundaries?
    
    % First the points outside the grid: 
    F_x(1,2) = ( ((3/4)*U_f(2,2)+(1/4)*U_f(2,3))-G.parent(x-1,y+1) )/(2*h);
    F_x(1,m_y-1) = ( ((3/4)*U_f(2,m_y-1)+(1/4)*U_f(2,m_y-2))-G.parent(x-1,y+(m_y-1)/2) )/(2*h);
    F_x(m_x,2) = ( G.parent(x+(m_x-1)/2,y+1)-((3/4)*U_f(m_x-1,2)+(1/4)*U_f(m_x-1,3)) )/(2*h);
    F_x(m_x,m_y-1) = ( G.parent(x+(m_x-1)/2+1,y+(m_y-1)/2-1)-((3/4)*U_f(m_x-1,m_y-1)+(1/4)*U_f(m_x-1,m_y-2)) )/(2*h);
    
    F_y(2,1) = ( ((3/4)*U_g(2,2)+(1/4)*U_g(3,2))-G.parent(x+1,y-1) )/(2*h);
    F_y(m_x-1,1) = ( ((3/4)*U_g(m_x-1,2)+(1/4)*U_g(m_x-2,2))-G.parent(x+1,y-1) )/(2*h);
    F_y(2,m_y) = ( G.parent(x+1,y+(m_y-1)/2+1)-((3/4)*U_g(2,m_y-1)+(1/4)*U_g(3,m_y-1)) )/(2*h);
    F_y(m_x-1,m_y) = ( G.parent(x+(m_x-1)/2-1,y+(m_y-1)/2+1)-((3/4)*U_g(m_x-1,m_y-1)+(1/4)*U_g(m_x-2,m_y-1)) )/(2*h);
    
    F_y(1,2) = ( U_g(1,3)-G.parent(x,y) )/(2*h);
    F_y(1,m_y-1) = ( G.parent(x,y+(m_y-1)/2)-U_g(1,m_y-2) )/(2*h);
    F_y(m_x,2) = ( U_g(m_x,3)-G.parent(x+(m_x-1)/2,y) )/(2*h);
    F_y(m_x,m_y-1) = ( G.parent(x+(m_x-1)/2,y+(m_y-1)/2)-U_g(m_x,m_y-2) )/(2*h);
    
    F_x(2,1) = ( U_g(3,1)-G.parent(x,y) )/(2*h);
    F_x(m_x-1,1) = ( G.parent(x+(m_x-1)/2,y)-U_g(m_x-2,1) )/(2*h);
    F_x(2,m_y) = ( U_g(3,m_y)-G.parent(x,y+(m_y-1)/2) )/(2*h);
    F_x(m_x-1,m_y) = ( G.parent(x+(m_x-1)/2,y+(m_y-1)/2)-U_g(m_x-2,m_y) )/(2*h);
    
    
    for i = 4:(m_x-4)/r
    
    for j = 4:(m_y-4)/r
        F_x(1,j) = ( U_f() )
    
    
    %...
    
    
    % First and last points are i = 1 and i = m_x. After that the fine grid
    % starts at i = 

    for i = 1:m_x % 1 + fine grid points
        for j = 1:m_y
            
            if (i == 1 && mod(j,r) == 1) %"Ghost-values" outside the fine grid. 
                
                if (j == 3)
                    F_x(i,j) = - ((U_f(i+2,j)*(3/4)+U_f(i+2,j+1)*(1/4))-f(G.parent.u(x-1,y+1)))/(4*h);
                    F_y(i,j) = - (U_g(i,j+2)-g(G.parent.u(x,y)))/(4*h);
                elseif (j == m_y)
                    F_x(i,j) = - ((U_f(i+2,j)*(1/2)+U_f(i+2,j-1)*(1/2)) - f(G.parent.u(x-1,y+(j-1)/2)))/(4*h); % ??
                    F_y(i,j) = - (U_g(i,j)-U_g(i,j-2))/(2*h);
                elseif (j > 3)
                    F_x(i,j) = - ((U_f(i+2,j)*(1/2)+U_f(i+2,j+1)*(1/4)+U_f(i+2,j-1)*(1/4))-f(G.parent.u(x-1,y+(j-1)/2)))/(4*h);
                    F_y(i,j) = - (U_g(i,j+2) - U_g(i,j-2))/(4*h);
                end
                
            elseif (i == 3 && j >= 3) % Points along the interface in y-direction
                
                if (j == 3)
                    F_x(i,j) = - (U_f(j,i+1) - U_f(j,i-2))/(3*h); % Hvorfor har jeg byttet om indexene her?
                elseif (mod(j,r) == 1)
                    F_x(i,j) = - (U_f(i+1,j) - U_f(i-2,j))/(3*h); % liten h, 2 fordi det går til begge sider.
                else
                    F_x(i,j) = -(U_f(i+1,j) - (U_f(i-2,j+1)+U_f(i-2,j-1))*(1/2))/(3*h);
                end
                
            elseif (i == m_x && j >= 3)
                F_x(i,j) = - (U_f(i,j) - U_f(i-1,j))/h;
            elseif (i > 3 && j >= 3)
                F_x(i,j) = - (U_f(i+1,j) - U_f(i-1,j))/(2*h); 
            end

            
            if (j == 1 && mod(i,r) == 1) % Points outside the fine grid. 
                
                if (i == 3) % First point on this line
                    F_x(i,j) = - ((U_f(i+2,j))-f(G.parent.u(x,y)))/(4*h);
                    F_y(i,j) = - ((U_g(i,j+2)*(3/4)+U_g(i+1,j+2)*(1/4))-g(G.parent.u(x+1,y-1)))/(4*h);
                elseif (i == m_x-2) % Last point on this line
                    F_x(i,j) = - (f(G.parent.u(xm,ym))-U_f(i-2,j))/(4*h); % ok
                    F_y(i,j) = - ((U_g(i,j+2)*(1/2)+U_g(i-1,j+2)*(1/2)) - g(G.parent.u(x+(i-1)/2,y-1)))/(4*h); % ??
                elseif (i > 3 && i < m_x-2)
                    F_x(i,j) = - (U_f(i+2,j) - U_f(i-2,j))/(4*h);
                    F_y(i,j) = - ((U_g(i,j+2)*(1/2)+U_g(i+1,j+2)*(1/4)+U_g(i-1,j+2)*(1/4))-g(G.parent.u(x+(i-1)/2,y-1)))/(4*h);
                end
                
            elseif (j == 3 && i >= 3 && i <= m_x-2) % points along the interface in x-direction
                
                if (i == 3)
                    F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-2))/(3*h);
                elseif (i == m_x-2)
                    F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-2))/(3*h);
                elseif (mod(i,r) == 1)
                    F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-2))/(3*h); % liten h, 2 fordi det går til begge sider.
                else
                    F_y(i,j) = - (U_g(i,j+1) - (U_g(i-1,j-2)*(1/2)+U_g(i+1,j-2)*(1/2)))/(3*h);
                end
                
            elseif (j == m_y && i >= 3)
                F_y(i,j) = - (U-g(i,j) - U_g(i,j-1))/h;
            elseif (j > 3 && j < (m_y - 2) && i >= 3) % within the boundary, on the y-boundary excluding the endpoints (because this is g flux)
                F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-1))/(2*h);
            end
            
        end
    end
    
    
    U_n = F_x + F_y;

    U = U_n;



%     h = G.h;
%     m_x = G.m;
%     m_y = G.m;
%     
%     x = G.location(3); % New: might be a better way
%     y = G.location(4);
%     
%     % Change the cases for on the boundary and outside the boundary. Also
%     % find the best way to combine these two rhs functions in the method.
%     % Remember to change the step size if needed. 
%     
%     F_x = zeros(m_x); % Number of grid points in fine grid p,uss an extra point outside the boundaries
%     F_y = zeros(m_y);
%     
%     [g_x,g_y] = fineBoundaryAdv(G,f,t);
%     
%     c=1;
%     if(G.parent ~= 0)
%         c=3;
%     end
% 
%     for i = 1:m_x
%         for j = 1:m_y
%             % When the grid starts two points away from the boundary
%             if (i == 1)
%                 F_x(i,j) = - a*(U(i+1,j) - g_y(j))/(c*h);
%             elseif (i == m_x)
%                 F_x(i,j) = - a*(U(i,j) - U(i-1,j))/h;
%             else
%                 F_x(i,j) = - a*(U(i+1,j) - U(i-1,j))/(2*h); 
%                 u1=U(i+1,j); % 0.9589
%                 u2=U(i-1,j); % 0.4794
%                 
%             end
% 
%             if (j == 1)
%                 F_y(i,j) = - b*(U(i,j+1) - g_x(i))/(c*h);
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
%     
%     
%     U_n = F_x + F_y;
% %     if (G.child ~= 0)
% %         disp(U_n)
% %     end
%     U = U_n;
    

end