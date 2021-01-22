function [g_x0,g_xm,g_y0,g_ym] = boundary(G,t)
% boundart function that calculates the vectors placed at boundaries of
% a grid. 
%   Creates boundaries from the current solution vector from parent. (If
%   parent is updated this could give a mistake as the numbers from which
%   the boundary should be created must be from the same time step as the
%   current one) 


% Must fix for mesh refinement due to bc on every side
    x0 = G.location(3);
    y0 = G.location(4); % Should be inserted under else
    
    x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
    y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';


    if G.parent == 0
%         g_x0 = boundaryValuesEuler(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m)',zeros(G.m,1),t); % Hva gjør jeg med tiden her?
%         g_xm = boundaryValuesEuler(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m)',ones(G.m,1),t); % Hva gjør jeg med tiden her?
%         g_y0 = boundaryValuesEuler(zeros(G.m,1),linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m)',t); % Hva gjør jeg med tiden her?
%         g_ym = boundaryValuesEuler(ones(G.m,1),linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m)',t); % Hva gjør jeg med tiden her?
        
        %u = initialConditionsEulerConstant(x,y,t);
        u = exactSolEuler(x,y,t);
        g_x0 = u(:,1,:);
        g_xm = u(:,G.m_x,:);
        g_y0 = u(1,:,:);
        g_ym = u(G.m_y,:,:);
    else
        
        r = G.parent.h/G.h; 
        
        xm = G.location(3)+(G.m_x-1)/r;
        ym = G.location(4)+(G.m_y-1)/r;
        
        g_x0 = zeros(G.m_x,4); % Spesifisert for Euler 
        g_xm = zeros(G.m_x,4);
        g_y0 = zeros(G.m_y,4);
        g_ym = zeros(G.m_y,4);
        
        for j = 1:4
            for i = 0:((G.m_x-1)/r)

                g_x0(i*r+1,j) = G.parent.u(x0+i,y0,j);
                g_xm(i*r+1,j) = G.parent.u(x0+i,ym,j); % 
                
            end
            for i = 0:((G.m_y-1)/r)
                
                g_y0(i*r+1,j) = G.parent.u(x0,y0+i,j);
                g_ym(i*r+1,j) = G.parent.u(xm,y0+i,j); %
                
            end

            i=2;
            gx0=[g_x0(1,j),g_x0(1+r,j)];
            gxm=[g_xm(1,j),g_xm(1+r,j)]; % 

            while i < G.m_x
                s = mod(i-1,r);

                g_x0(i,j) = ((r-s)*gx0(1) + s*gx0(2))/r;
                
                g_xm(i,j) = ((r-s)*gxm(1) + s*gxm(2))/r; % 

                if (s == r-1) && (i+r+1 <= G.m_x)
                    gx0 = [gx0(2),g_x0(i+r+1,j)];
                    
                    gxm = [gxm(2),g_xm(i+r+1,j)]; % 
                    
                    i=i+1;
                end
                i=i+1;
            end
            
            i=2;
            gy0=[g_y0(1,j),g_y0(1+r,j)];
            gym=[g_ym(1,j),g_ym(1+r,j)]; %
            
            while i < G.m_y
                s = mod(i-1,r);
                
                g_y0(i,j) = ((r-s)*gy0(1) + s*gy0(2))/r;
                
                g_ym(i,j) = ((r-s)*gym(1) + s*gym(2))/r; %

                if (s == r-1) && (i+r+1 <= G.m_y)
                    gy0 = [gy0(2),g_y0(i+r+1,j)];
                    
                    gym = [gym(2),g_ym(i+r+1,j)]; %
                    
                    i=i+1;
                end
                i=i+1;
            end
        end
    end
    
%     %
%     u = exactSolEuler(x,y,t);
%         g_x0 = u(:,1,:);
%         g_xm = u(:,G.m_x,:);
%         g_y0 = u(1,:,:);
%         g_ym = u(G.m_y,:,:);
%     %
    
end