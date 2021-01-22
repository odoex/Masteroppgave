function [U] = rhsSubgrid(U,G)
% Only works for a refinement of r=2

    h = G.h;
    m_x = G.m_x;
    m_y = G.m_y;
    r = G.parent.h/G.h; 
    
    delta = 147*(G.parent.h); % Sjekker først med helt lik delta på grovt og fint grid. 
    
    x = G.location(3); 
    y = G.location(4);
    xm = G.location(3) + (m_x-1)/r; 
    ym = G.location(4) + (m_y-1)/r;
    
    variables = length(U(1,1,:));
    F_x = zeros(m_x,m_y,variables); % Number of grid points in fine grid pluss an extra point outside the boundaries
    F_y = zeros(m_x,m_y,variables); 
    
    U_f = f(U);
    U_g = g(U);
    
    Up_f = f(G.parent.u);
    Up_g = g(G.parent.u);
    
    for l = 1:variables
        % First the points outside the grid: 
        F_x(1,2,l) = ( ((3/4)*U_f(2,2,l)+(1/4)*U_f(2,3,l))-(Up_f(x-1,y+1,l)) )/(4*h);
        F_x(1,m_y-1,l) = ( ((3/4)*U_f(2,m_y-1,l)+(1/4)*U_f(2,m_y-2,l))-(Up_f(x-1,y+(m_y-3)/2+1,l)) )/(4*h);
        F_x(m_x,2,l) = ( (Up_f(x+(m_x-3)/2+3,y+1,l))-((3/4)*U_f(m_x-1,2,l)+(1/4)*U_f(m_x-1,3,l)) )/(4*h);
        F_x(m_x,m_y-1,l) = ( (Up_f(x+(m_x-3)/r+3,y+(m_y-3)/2+1,l))-((3/4)*U_f(m_x-1,m_y-1,l)+(1/4)*U_f(m_x-1,m_y-2,l)) )/(4*h);

        F_y(2,1,l) = ( ((3/4)*U_g(2,2,l)+(1/4)*U_g(3,2,l))-(Up_g(x+1,y-1,l)) )/(4*h);
        %F_y(m_x-1,1,l) = ( ((3/4)*U_g(m_x-1,2,l)+(1/4)*U_g(m_x-2,2,l))-(Up_g(x+(m_x-3)/r+2,y-1,l)) )/(4*h);
        F_y(m_x-1,1,l) = ( ((3/4)*U_g(m_x-1,2,l)+(1/4)*U_g(m_x-2,2,l))-(Up_g(x+(m_x-3)/r+1,y-1,l)) )/(4*h);
        F_y(2,m_y,l) = ( (Up_g(x+1,y+(m_y-3)/r+3,l))-((3/4)*U_g(2,m_y-1,l)+(1/4)*U_g(3,m_y-1,l)) )/(4*h);
        %F_y(m_x-1,m_y,l) = ( (Up_g(x+(m_x-3)/2+3,y+(m_y-3)/2+1,l))-((3/4)*U_g(m_x-1,m_y-1,l)+(1/4)*U_g(m_x-2,m_y-1,l)) )/(4*h); %%% Feil? 
        F_y(m_x-1,m_y,l) = ( (Up_g(x+(m_x-3)/r+1,y+(m_y-3)/r+3,l))-((3/4)*U_g(m_x-1,m_y-1,l)+(1/4)*U_g(m_x-2,m_y-1,l)) )/(4*h);

        F_y(1,2,l) = ( U_g(1,4,l)-(Up_g(x,y,l)) )/(4*h);
        F_y(1,m_y-1,l) = ( (Up_g(x,y+(m_y-3)/2+2,l))-U_g(1,m_y-3,l) )/(4*h);
        F_y(m_x,2,l) = ( U_g(m_x,4,l)-(Up_g(x+(m_x-3)/2+2,y,l)) )/(4*h);
        F_y(m_x,m_y-1,l) = ( (Up_g(x+(m_x-3)/r+2,y+(m_y-3)/2+2,l))-U_g(m_x,m_y-3,l) )/(4*h);

        F_x(2,1,l) = ( U_f(4,1,l)-(Up_f(x,y,l)) )/(4*h);
        F_x(m_x-1,1,l) = ( (Up_f(x+(m_x-3)/r+2,y,l))-U_f(m_x-3,1,l) )/(4*h);
        F_x(2,m_y,l) = ( U_f(4,m_y,l)-(Up_f(x,y+(m_y-3)/r+2,l)) )/(4*h);
        F_x(m_x-1,m_y,l) = ( (Up_f(x+(m_x-3)/2+2,y+(m_y-3)/2+2,l))-U_f(m_x-3,m_y,l) )/(4*h);


        for i = 4:r:(m_x-3) % Change 4 with different ratio
            F_x(i,1,l) = (U_f(i+r,1,l)-U_f(i-r,1,l))/(4*h);
            F_x(i,m_y,l) = (U_f(i+r,m_y,l)-U_f(i-r,m_y,l))/(4*h);
            F_y(i,1,l) = ( ((1/2)*U_g(i,2,l)+(1/4)*U_g(i+1,2,l)+(1/4)*U_g(i-1,2,l)) - (Up_g(x+i/r,y-1,l)) )/(4*h); % Change with different r
            F_y(i,m_y,l) = ( -((1/2)*U_g(i,m_y-1,l)+(1/4)*U_g(i+1,m_y-1,l)+(1/4)*U_g(i-1,m_y-1,l)) + (Up_g(x+i/r,ym+2,l)) )/(4*h);
        end

        for j = 4:r:(m_y-3)
            F_x(1,j,l) = ( (1/2)*U_f(2,j,l)+(1/4)*U_f(2,j+1,l)+(1/4)*U_f(2,j-1,l) - (Up_f(x-1,y+j/r,l)) )/(4*h);
            F_x(m_x,j,l) = ( -((1/2)*U_f(m_x-1,j,l)+(1/4)*U_f(m_x-1,j+1,l)+(1/4)*U_f(m_x-1,j-1,l)) + (Up_f(xm+2,y+j/r,l)) )/(4*h);
            F_y(1,j,l) = ( U_g(1,j+r,l) - U_g(1,j-r,l) )/(4*h);
            F_y(m_x,j,l) = ( U_g(m_x,j+r,l) - U_g(m_x,j-r,l) )/(4*h);
        end

        % Next: boundary. First for refinement 2 only

        F_x(2,2,l) = (U_f(3,2,l)-U_f(1,2,l))/(3*h);
        F_y(2,2,l) = (U_g(2,3,l)-U_g(2,1,l))/(3*h);

        F_x(2,m_y-1,l) = (U_f(3,m_y-1,l)-U_f(1,m_y-1,l))/(3*h);
        F_y(2,m_y-1,l) = (U_g(2,m_y,l)-U_g(2,m_y-2,l))/(3*h);

        F_x(m_x-1,2,l) = (U_f(m_x,2,l)-U_f(m_x-2,2,l))/(3*h);
        F_y(m_x-1,2,l) = (U_g(m_x-1,3,l)-U_g(m_x-1,1,l))/(3*h);

        F_x(m_x-1,m_y-1,l) = (U_f(m_x,m_y-1,l)-U_f(m_x-2,m_y-1,l))/(3*h);
        F_y(m_x-1,m_y-1,l) = (U_g(m_x-1,m_y,l)-U_g(m_x-1,m_y-2,l))/(3*h);

        for i = 3:(m_x-2)
            if (mod((i-2),r) == 0)
                F_y(i,2,l) = (U_g(i,3,l)-U_g(i,1,l))/(3*h);
                F_y(i,m_y-1,l) = (U_g(i,m_y,l)-U_g(i,m_y-2,l))/(3*h);
            else
                F_y(i,2,l) = (U_g(i,3,l)-(1/2)*(U_g(i+1,1,l)+U_g(i-1,1,l)))/(3*h);
                F_y(i,m_y-1,l) = (-U_g(i,m_y-2,l)+(1/2)*(U_g(i+1,m_y,l)+U_g(i-1,m_y,l)))/(3*h);
            end

            F_x(i,2,l) = (U_f(i+1,2,l)-U_f(i-1,2,l))/(2*h);
            F_x(i,m_y-1,l) = (U_f(i+1,m_y-1,l)-U_f(i-1,m_y-1,l))/(2*h);

        end

        for j = 3:(m_x-2)
            if (mod((j-2),r) == 0)
                F_x(2,j,l) = (U_f(3,j,l)-U_f(1,j,l))/(3*h); % Different for r != 2
                F_x(m_x-1,j,l) = (U_f(m_x,j,l)-U_f(m_x-2,j,l))/(3*h);
            else
                F_x(2,j,l) = (U_f(3,j,l)-(1/2)*(U_f(1,j+1,l)+U_f(1,j-1,l)))/(3*h); % Change these two conditions in some way
                F_x(m_x-1,j,l) = (-U_f(m_x-2,j,l)+(1/2)*(U_f(m_x,j+1,l)+U_f(m_x,j-1,l)))/(3*h);
            end

            F_y(2,j,l) = (U_g(2,j+1,l)-U_g(2,j-1,l))/(2*h);
            F_y(m_x-1,j,l) = (U_g(m_x-1,j+1,l)-U_g(m_x-1,j-1,l))/(2*h);
        end


        % Inner points on the fine grid: 

        for i = 3:m_x-2 
            for j = 3:m_y-2
                F_x(i,j,l) = (U_f(i+1,j,l)-U_f(i-1,j,l))/(2*h);
                F_y(i,j,l) = (U_g(i,j+1,l)-U_g(i,j-1,l))/(2*h);
            end
        end
    end

    diff_f = diffusion(delta,U,G.parent.u,x,y,r,h);
    diff_g = diffusion(delta,permute(U,[2,1,3]),permute(G.parent.u,[2,1,3]),y,x,r,h);
    
    U_n = - (F_x + diff_f) - (F_y + permute(diff_g,[2,1,3]));
    %U_n =  - (F_y + permute(diff_g,[2,1,3]));

    U = U_n;

end