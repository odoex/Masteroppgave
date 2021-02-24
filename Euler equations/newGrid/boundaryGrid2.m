function [g_x,g_y] = boundaryGrid2(G,t)
% boundart function that calculates the vectors placed at boundaries of
% a grid. 
%   Creates boundaries from the current solution vector from parent. (If
%   parent is updated this could give a mistake as the numbers from which
%   the boundary should be created must be from the same time step as the
%   current one) 

% The function has an oblique boundary, and the function must find the best
% way to approximate this boundary and impose the most accurate boundary
% conditions
    
    x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
    y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
    
	if G.parent == 0

        g_x(:,1,:) = exactSolEuler(x,0,t);
        g_x(:,2,:) = exactSolEuler(x,1,t);
        g_y(1,:,:) = exactSolEuler(0,y,t);
        g_y(2,:,:) = exactSolEuler(1,y,t);

%         g_x(:,1,:) = initialConditionsEulerConstant(x,0,t);
%         g_x(:,2,:) = initialConditionsEulerConstant(x,1,t);
%         g_y(1,:,:) = initialConditionsEulerConstant(0,y,t);
%         g_y(2,:,:) = initialConditionsEulerConstant(1,y,t);
        
	else
        r = G.parent.h/G.h; 
        
        x0 = G.location(3);
        y0 = G.location(4);
        xm = G.location(3)+(G.m_x-1)/r;
        ym = G.location(4)+(G.m_y-1)/r;
        
        g_x = zeros(G.m_x,2,4); 
        g_y = zeros(2,G.m_y,4);
        
        g_x(1:r:G.m_x,1,:) = G.parent.u(x0:xm,y0,:);
        g_x(1:r:G.m_x,2,:) = G.parent.u(x0:xm,ym,:);
        
        g_y(1,1:r:G.m_y,:) = G.parent.u(x0,y0:ym,:);
        g_y(2,1:r:G.m_y,:) = G.parent.u(xm,y0:ym,:);
        
        i=2;
        gx0=[g_x(1,1,:),g_x(1+r,1,:)];
        gxm=[g_x(1,2,:),g_x(1+r,2,:)];
        
        while i < G.m_x
            s = mod(i-1,r);

            g_x(i,1,:) = ((r-s)*gx0(:,1,:) + s*gx0(:,2,:))/r;
                
            g_x(i,2,:) = ((r-s)*gxm(:,1,:) + s*gxm(:,2,:))/r;

            if (s == r-1) && (i+r+1 <= G.m_x)
                gx0 = [gx0(:,2,:),g_x(i+r+1,1,:)];
                    
                gxm = [gxm(:,2,:),g_x(i+r+1,2,:)];
                    
                i=i+1;
            end
            i=i+1;
        end
            
        i=2;
        gy0=[g_y(1,1,:),g_y(1,1+r,:)];
        gym=[g_y(2,1,:),g_y(2,1+r,:)];
            
        while i < G.m_y
            s = mod(i-1,r);
                
            g_y(1,i,:) = ((r-s)*gy0(:,1,:) + s*gy0(:,2,:))/r;
                
            g_y(2,i,:) = ((r-s)*gym(:,1,:) + s*gym(:,2,:))/r;

            if (s == r-1) && (i+r+1 <= G.m_y)
                gy0 = [gy0(:,2,:),g_y(1,i+r+1,:)];
                    
                gym = [gym(:,2,:),g_y(2,i+r+1,:)];
                    
                i=i+1;
            end
            i=i+1;
        end
        
%         for v = 1:4
%             for i = 0:((G.m_x-1)/r)
% 
%                 g_x(i*r+1,1,v) = G.parent.u(x0+i,y0,v);
%                 g_x(i*r+1,2,v) = G.parent.u(x0+i,ym,v);
%                 
%             end
%             for i = 0:((G.m_y-1)/r)
%                 
%                 g_y(1,i*r+1,v) = G.parent.u(x0,y0+i,v);
%                 g_y(2,i*r+1,v) = G.parent.u(xm,y0+i,v);
%                 
%             end
% 
%             i=2;
%             gx0=[g_x(1,1,v),g_x(1+r,1,v)];
%             gxm=[g_x(1,2,v),g_x(1+r,2,v)];
% 
%             while i < G.m_x
%                 s = mod(i-1,r);
% 
%                 g_x(i,1,v) = ((r-s)*gx0(1) + s*gx0(2))/r;
%                 
%                 g_x(i,2,v) = ((r-s)*gxm(1) + s*gxm(2))/r;
% 
%                 if (s == r-1) && (i+r+1 <= G.m_x)
%                     gx0 = [gx0(2),g_x(i+r+1,1,v)];
%                     
%                     gxm = [gxm(2),g_x(i+r+1,2,v)];
%                     
%                     i=i+1;
%                 end
%                 i=i+1;
%             end
%             
%             i=2;
%             gy0=[g_y(1,1,v),g_y(1,1+r,v)];
%             gym=[g_y(2,2,v),g_y(2,1+r,v)];
%             
%             while i < G.m_y
%                 s = mod(i-1,r);
%                 
%                 g_y(1,i,v) = ((r-s)*gy0(1) + s*gy0(2))/r;
%                 
%                 g_y(2,i,v) = ((r-s)*gym(1) + s*gym(2))/r;
% 
%                 if (s == r-1) && (i+r+1 <= G.m_y)
%                     gy0 = [gy0(2),g_y(1,i+r+1,v)];
%                     
%                     gym = [gym(2),g_y(2,i+r+1,v)];
%                     
%                     i=i+1;
%                 end
%                 i=i+1;
%             end
%         end
	end
    
end