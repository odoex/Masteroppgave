function [U] = rhsFineGridAdv(U,t,G)
% Only works for a refinement of r=2

% Hva er raskest av å sende hele matrisen gjennom funksjonen eller sende
% hvert av elementene (Kun 12 tror jeg) gjennom funksjonen hver gang de
% skal brukes?
    delta = 0;%0*G.h;
    h = G.h;
    m_x = G.m_x;
    m_y = G.m_y;
    r = G.parent.h/G.h; % So much to do
    
    x = G.location(3); 
    y = G.location(4);
    xm = G.location(3) + (m_x-1)/r; 
    ym = G.location(4) + (m_y-1)/r; % BRUK DISSE
    
    % Change the cases for on the boundary and outside the boundary. Also
    % find the best way to combine these two rhs functions in the method.
    % Remember to change the step size if needed. 
    
    F_x = zeros(m_x); % Number of grid points in fine grid pluss an extra point outside the boundaries
    F_y = zeros(m_y); 
    
%     [g_x,g_y] = fineBoundaryAdv(G,f,t);
    
    U_f = f_adv(U);
    U_g = g_adv(U);
    
    Up_f = f_adv(G.parent.u);
    Up_g = g_adv(G.parent.u);
    
    % Dobbeltsjekk h
    
    %...
    
    % Creating two vectors b_x and b_y for boundaries?
    
    % First the points outside the grid: 
    F_x(1,2) = ( ((3/4)*U_f(2,2)+(1/4)*U_f(2,3))-f_adv(G.parent.u(x-1,y+1)) )/(4*h);
    F_x(1,m_y-1) = ( ((3/4)*U_f(2,m_y-1)+(1/4)*U_f(2,m_y-2))-f_adv(G.parent.u(x-1,y+(m_y-3)/2+1)) )/(4*h);
    F_x(m_x,2) = ( f_adv(G.parent.u(x+(m_x-3)/2+3,y+1))-((3/4)*U_f(m_x-1,2)+(1/4)*U_f(m_x-1,3)) )/(4*h);
    F_x(m_x,m_y-1) = ( f_adv(G.parent.u(x+(m_x-3)/r+3,y+(m_y-3)/2+1))-((3/4)*U_f(m_x-1,m_y-1)+(1/4)*U_f(m_x-1,m_y-2)) )/(4*h);
    
    F_y(2,1) = ( ((3/4)*U_g(2,2)+(1/4)*U_g(3,2))-g_adv(G.parent.u(x+1,y-1)) )/(4*h);
    %F_y(m_x-1,1) = ( ((3/4)*U_g(m_x-1,2)+(1/4)*U_g(m_x-2,2))-g_adv(G.parent.u(x+(m_x-3)/r+2,y-1)) )/(4*h);
    F_y(m_x-1,1) = ( ((3/4)*U_g(m_x-1,2)+(1/4)*U_g(m_x-2,2))-g_adv(G.parent.u(x+(m_x-3)/r+1,y-1)) )/(4*h); % Mener denne er riktig
    F_y(2,m_y) = ( g_adv(G.parent.u(x+1,y+(m_y-3)/r+3))-((3/4)*U_g(2,m_y-1)+(1/4)*U_g(3,m_y-1)) )/(4*h);
    %F_y(m_x-1,m_y) = ( g_adv(G.parent.u(x+(m_x-3)/2+3,y+(m_y-3)/2+1))-((3/4)*U_g(m_x-1,m_y-1)+(1/4)*U_g(m_x-2,m_y-1)) )/(4*h);
    F_y(m_x-1,m_y) = ( g_adv(G.parent.u(x+(m_x-3)/2+1,y+(m_y-3)/2+3))-((3/4)*U_g(m_x-1,m_y-1)+(1/4)*U_g(m_x-2,m_y-1)) )/(4*h); % Mener det er denne som er riktig 
    
    F_y(1,2) = ( U_g(1,4)-g_adv(G.parent.u(x,y)) )/(4*h);
    F_y(1,m_y-1) = ( g_adv(G.parent.u(x,y+(m_y-3)/2+2))-U_g(1,m_y-3) )/(4*h);
    F_y(m_x,2) = ( U_g(m_x,4)-g_adv(G.parent.u(x+(m_x-3)/2+2,y)) )/(4*h);
    F_y(m_x,m_y-1) = ( g_adv(G.parent.u(x+(m_x-3)/r+2,y+(m_y-3)/2+2))-U_g(m_x,m_y-3) )/(4*h);
    
    F_x(2,1) = ( U_f(4,1)-f_adv(G.parent.u(x,y)) )/(4*h);
    F_x(m_x-1,1) = ( f_adv(G.parent.u(x+(m_x-3)/r+2,y))-U_f(m_x-3,1) )/(4*h);
    F_x(2,m_y) = ( U_f(4,m_y)-f_adv(G.parent.u(x,y+(m_y-3)/r+2)) )/(4*h);
    F_x(m_x-1,m_y) = ( f_adv(G.parent.u(x+(m_x-3)/2+2,y+(m_y-3)/2+2))-U_f(m_x-3,m_y) )/(4*h);
    
    
    for i = 4:r:(m_x-3) % Change 4 with different ratio
        
        F_x(i,1) = (U_f(i+r,1)-U_f(i-r,1))/(4*h);
        F_x(i,m_y) = (U_f(i+r,m_y)-U_f(i-r,m_y))/(4*h);
        F_y(i,1) = ( ((1/2)*U_g(i,2)+(1/4)*U_g(i+1,2)+(1/4)*U_g(i-1,2)) - g_adv(G.parent.u(x+i/r,y-1)) )/(4*h); % Change with different r
        F_y(i,m_y) = ( -((1/2)*U_g(i,m_y-1)+(1/4)*U_g(i+1,m_y-1)+(1/4)*U_g(i-1,m_y-1)) + g_adv(G.parent.u(x+i/r,ym+2)) )/(4*h);
        
%         F_x(i,1) = ( (1/2)*U_f(r*(i-4)+4,2)+(1/4)*U_f(r*(i-4)+5,2)+(1/4)*U_f(r*(i-4)+3,2) - f_adv(G.parent.u(x+i-1,y-1)) )/(4*h);
%         F_y(i,1) = ( U_g(r*(i-4)+6,1) - U_g(r*(i-4)+2,1) )/(4*h);
%         F_x(i,m_y) = ( -(1/2)*U_f(r*(i-4)+4,m_y-1)+(1/4)*U_f(r*(i-4)+5,m_y-1)+(1/4)*U_f(r*(i-4)+3,m_y-1) + f_adv(G.parent.u(x+i-1,y+(m_y-1)/r+1)) )/(4*h);
%         F_y(i,m_y) = ( U_g(r*(i-4)+6,m_y) - U_g(r*(i-4)+2,m_y) )/(4*h);
    end
    
    for j = 4:r:(m_y-3)
        F_x(1,j) = ( (1/2)*U_f(2,j)+(1/4)*U_f(2,j+1)+(1/4)*U_f(2,j-1) - f_adv(G.parent.u(x-1,y+j/r)) )/(4*h);
        F_x(m_x,j) = ( -((1/2)*U_f(m_x-1,j)+(1/4)*U_f(m_x-1,j+1)+(1/4)*U_f(m_x-1,j-1)) + f_adv(G.parent.u(xm+2,y+j/r)) )/(4*h);
        F_y(1,j) = ( U_g(1,j+r) - U_g(1,j-r) )/(4*h);
        F_y(m_x,j) = ( U_g(m_x,j+r) - U_g(m_x,j-r) )/(4*h);
    end
    
    %...
    
    % Next: boundary. First for refinement 2 only
    
    F_x(2,2) = (U_f(3,2)-U_f(1,2))/(3*h);
    F_y(2,2) = (U_g(2,3)-U_g(2,1))/(3*h);
    
    F_x(2,m_y-1) = (U_f(3,m_y-1)-U_f(1,m_y-1))/(3*h);
    F_y(2,m_y-1) = (U_g(2,m_y)-U_g(2,m_y-2))/(3*h);
    
    F_x(m_x-1,2) = (U_f(m_x,2)-U_f(m_x-2,2))/(3*h);
    F_y(m_x-1,2) = (U_g(m_x-1,3)-U_g(m_x-1,1))/(3*h);
    
    F_x(m_x-1,m_y-1) = (U_f(m_x,m_y-1)-U_f(m_x-2,m_y-1))/(3*h);
    F_y(m_x-1,m_y-1) = (U_g(m_x-1,m_y)-U_g(m_x-1,m_y-2))/(3*h);
    
    for i = 3:(m_x-2)
        if (mod((i-2),r) == 0)
            F_y(i,2) = (U_g(i,3)-U_g(i,1))/(3*h);
            F_y(i,m_y-1) = (U_g(i,m_y)-U_g(i,m_y-2))/(3*h);
        else
            F_y(i,2) = (U_g(i,3)-(1/2)*(U_g(i+1,1)+U_g(i-1,1)))/(3*h);
            F_y(i,m_y-1) = (-U_g(i,m_y-2)+(1/2)*(U_g(i+1,m_y)+U_g(i-1,m_y)))/(3*h);
        end
        
        F_x(i,2) = (U_f(i+1,2)-U_f(i-1,2))/(2*h);
        F_x(i,m_y-1) = (U_f(i+1,m_y-1)-U_f(i-1,m_y-1))/(2*h);
        
    end
    
    for j = 3:(m_x-2)
        if (mod((j-2),r) == 0)
            F_x(2,j) = (U_f(3,j)-U_f(1,j))/(3*h); % Different for r != 2
            F_x(m_x-1,j) = (U_f(m_x,j)-U_f(m_x-2,j))/(3*h);
        else
            F_x(2,j) = (U_f(3,j)-(1/2)*(U_f(1,j+1)+U_f(1,j-1)))/(3*h); % Change these two conditions in some way
            F_x(m_x-1,j) = (-U_f(m_x-2,j)+(1/2)*(U_f(m_x,j+1)+U_f(m_x,j-1)))/(3*h);
        end
        
        F_y(2,j) = (U_g(2,j+1)-U_g(2,j-1))/(2*h);
        F_y(m_x-1,j) = (U_g(m_x-1,j+1)-U_g(m_x-1,j-1))/(2*h);
    end
    
    
    % Inner points on the fine grid: 

    for i = 3:m_x-2 
        for j = 3:m_y-2
            F_x(i,j) = (U_f(i+1,j)-U_f(i-1,j))/(2*h);
            F_y(i,j) = (U_g(i,j+1)-U_g(i,j-1))/(2*h);
        end
    end

    % Trying something else to look for the error, removing some functions
    % where the error may be. Remove when the error is found.
    
%     for i = 3:(m_x-2)
%         for j = 3:(m_y-2)
%             F_x(i,j) = - 0.5*(U(i+1,j)-U(i-1,j))/(2*G.h) - 1*(U(i,j+1)-U(i,j-1))/(2*G.h);
%             F_y(i,j)=0;
%         end
%     end
    
%     for i = 1:m_x 
%         for j = 1:m_y
%             
%             if (i == 1 && mod(j,r) == 1) %"Ghost-values" outside the fine grid. 
%                 
%                 if (j == 3)
%                     F_x(i,j) = - ((U_f(i+2,j)*(3/4)+U_f(i+2,j+1)*(1/4))-f(G.parent.u(x-1,y+1)))/(4*h);
%                     F_y(i,j) = - (U_g(i,j+2)-g(G.parent.u(x,y)))/(4*h);
%                 elseif (j == m_y)
%                     F_x(i,j) = - ((U_f(i+2,j)*(1/2)+U_f(i+2,j-1)*(1/2)) - f(G.parent.u(x-1,y+(j-1)/2)))/(4*h); % ??
%                     F_y(i,j) = - (U_g(i,j)-U_g(i,j-2))/(2*h);
%                 elseif (j > 3)
%                     F_x(i,j) = - ((U_f(i+2,j)*(1/2)+U_f(i+2,j+1)*(1/4)+U_f(i+2,j-1)*(1/4))-f(G.parent.u(x-1,y+(j-1)/2)))/(4*h);
%                     F_y(i,j) = - (U_g(i,j+2) - U_g(i,j-2))/(4*h);
%                 end
%                 
%             elseif (i == 3 && j >= 3) % Points along the interface in y-direction
%                 
%                 if (j == 3)
%                     F_x(i,j) = - (U_f(j,i+1) - U_f(j,i-2))/(3*h); % Hvorfor har jeg byttet om indexene her?
%                 elseif (mod(j,r) == 1)
%                     F_x(i,j) = - (U_f(i+1,j) - U_f(i-2,j))/(3*h); % liten h, 2 fordi det går til begge sider.
%                 else
%                     F_x(i,j) = -(U_f(i+1,j) - (U_f(i-2,j+1)+U_f(i-2,j-1))*(1/2))/(3*h);
%                 end
%                 
%             elseif (i == m_x && j >= 3)
%                 F_x(i,j) = - (U_f(i,j) - U_f(i-1,j))/h;
%             elseif (i > 3 && j >= 3)
%                 F_x(i,j) = - (U_f(i+1,j) - U_f(i-1,j))/(2*h); 
%             end
% 
%             
%             if (j == 1 && mod(i,r) == 1) % Points outside the fine grid. 
%                 
%                 if (i == 3) % First point on this line
%                     F_x(i,j) = - ((U_f(i+2,j))-f(G.parent.u(x,y)))/(4*h);
%                     F_y(i,j) = - ((U_g(i,j+2)*(3/4)+U_g(i+1,j+2)*(1/4))-g(G.parent.u(x+1,y-1)))/(4*h);
%                 elseif (i == m_x-2) % Last point on this line
%                     F_x(i,j) = - (f(G.parent.u(xm,ym))-U_f(i-2,j))/(4*h); % ok
%                     F_y(i,j) = - ((U_g(i,j+2)*(1/2)+U_g(i-1,j+2)*(1/2)) - g(G.parent.u(x+(i-1)/2,y-1)))/(4*h); % ??
%                 elseif (i > 3 && i < m_x-2)
%                     F_x(i,j) = - (U_f(i+2,j) - U_f(i-2,j))/(4*h);
%                     F_y(i,j) = - ((U_g(i,j+2)*(1/2)+U_g(i+1,j+2)*(1/4)+U_g(i-1,j+2)*(1/4))-g(G.parent.u(x+(i-1)/2,y-1)))/(4*h);
%                 end
%                 
%             elseif (j == 3 && i >= 3 && i <= m_x-2) % points along the interface in x-direction
%                 
%                 if (i == 3)
%                     F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-2))/(3*h);
%                 elseif (i == m_x-2)
%                     F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-2))/(3*h);
%                 elseif (mod(i,r) == 1)
%                     F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-2))/(3*h); % liten h, 2 fordi det går til begge sider.
%                 else
%                     F_y(i,j) = - (U_g(i,j+1) - (U_g(i-1,j-2)*(1/2)+U_g(i+1,j-2)*(1/2)))/(3*h);
%                 end
%                 
%             elseif (j == m_y && i >= 3)
%                 F_y(i,j) = - (U-g(i,j) - U_g(i,j-1))/h;
%             elseif (j > 3 && j < (m_y - 2) && i >= 3) % within the boundary, on the y-boundary excluding the endpoints (because this is g flux)
%                 F_y(i,j) = - (U_g(i,j+1) - U_g(i,j-1))/(2*h);
%             end
%             
%         end
%     end
    
    Ux(:,:,1) = U;
    Uy(:,:,1) = U';
    
    U_diff_x = diffusion(delta,Ux,G.parent.u,x,y,r,h);
    U_diff_y = diffusion(delta,Uy,G.parent.u',y,x,r,h);
    
    U_n = - F_x-U_diff_x - F_y-U_diff_y';
    %U_n = - F_y-U_diff_y';

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