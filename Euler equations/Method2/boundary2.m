function [g_x,g_y] = boundary2(G,t)
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
        g_x = [exactSolEuler(x,0,t),exactSolEuler(x,1,t)];
        g_y = [exactSolEuler(0,y,t);exactSolEuler(1,y,t)];
        
%         g_x = [initialConditionsEulerConstant(x,0,t),initialConditionsEulerConstant(x,1,t)];
%         g_y = [initialConditionsEulerConstant(0,y,t);initialConditionsEulerConstant(1,y,t)];
        
        g_x = g(g_x);
        g_y = f(g_y);
        
%         u = exactSolEuler(x,y,t);
%         g_x0 = u(:,1,:);
%         g_xm = u(:,G.m_x,:);
%         g_y0 = u(1,:,:);
%         g_ym = u(G.m_y,:,:);
    else
        % Linear extrapolation
        
        delta = 0;
        r = G.parent.h/G.h; 
        
        xm = G.location(3)+G.m_x/r-1;
        ym = G.location(4)+G.m_y/r-1;
        
        g_x(:,:,:) = zeros(G.m_x,2,4);
        g_y(:,:,:) = zeros(2,G.m_y,4);
        
%         for i=1:r
%             g_x(i:r:end-r+i,1,:) = (1/2)*(g(G.parent.u(x0:xm,y0-1,:)) - (delta/2)*(G.u(i:r:end-r+i,1,:)-G.parent.u(x0:xm,y0-1,:)) + ( g(G.u(i:r:end-r+i,2,:)) - (delta/2)*(G.u(i:r:end-r+i,2,:)-G.u(i:r:end-r+i,1,:)) ) );
%             g_x(i:r:end-r+i,2,:) = (1/2)*(g(G.parent.u(x0:xm,ym+1,:)) - (delta/2)*(G.parent.u(x0:xm,ym+1,:)-G.u(i:r:end-r+i,end,:)) + ( g(G.u(i:r:end-r+i,end-1,:)) - (delta/2)*(G.u(i:r:end-r+i,end,:)-G.u(i:r:end-r+i,end-1,:)) ) );
%             g_y(1,i:r:end-r+i,:) = (1/2)*(f(G.parent.u(x0-1,y0:ym,:)) - (delta/2)*(G.u(1,i:r:end-r+i,:)-G.parent.u(x0-1,y0:ym,:)) + ( f(G.u(2,i:r:end-r+i,:)) - (delta/2)*(G.u(2,i:r:end-r+i,:)-G.u(1,i:r:end-r+i,:)) ) );
%             g_y(2,i:r:end-r+i,:) = (1/2)*(f(G.parent.u(xm+1,y0:ym,:)) - (delta/2)*(G.parent.u(xm+1,y0:ym,:)-G.u(end,i:r:end-r+i,:)) + ( f(G.u(end-1,i:r:end-r+i,:)) - (delta/2)*(G.u(end,i:r:end-r+i,:)-G.u(end-1,i:r:end-r+i,:)) ) );
%         end
        
%         for i=1:r
%             g_x(i:r:end-r+i,1,:) = g( (1/2)*(G.parent.u(x0:xm,y0-1,:) + G.u(i:r:end-r+i,2,:)) );
%             g_x(i:r:end-r+i,2,:) = g( (1/2)*(G.parent.u(x0:xm,ym+1,:) + G.u(i:r:end-r+i,end-1,:)) );
%             g_y(1,i:r:end-r+i,:) = f( (1/2)*(G.parent.u(x0-1,y0:ym,:) + G.u(2,i:r:end-r+i,:)) );
%             g_y(2,i:r:end-r+i,:) = f( (1/2)*(G.parent.u(xm+1,y0:ym,:) + G.u(end-1,i:r:end-r+i,:)) );
%         end
        
%         for i=1:r
%             g_x(i:r:end-r+i,1,:) = g( G.parent.u(x0:xm,y0,:) );
%             g_x(i:r:end-r+i,2,:) = g( G.parent.u(x0:xm,ym,:) );
%             g_y(1,i:r:end-r+i,:) = f( G.parent.u(x0,y0:ym,:) );
%             g_y(2,i:r:end-r+i,:) = f( G.parent.u(xm,y0:ym,:) );
%         end
        
        h = G.parent.h;
        
        c_y1 = [(G.h/2)/h,(h-G.h/2)/h];
        c_y2 = [(h-G.h/2)/h,(G.h/2)/h];
        
        c_x1 = [((G.h/2)/h),(h-G.h/2)/h];
        c_x2 = [(h-G.h/2)/h,((G.h/2)/h)];
        
        for i=1:r
            g_x(i:r:end-r+i,1,:) = g( c_y1(1).*(c_x1(i).*G.parent.u(x0-2+i:xm-2+i,y0-1,:) + c_x2(i).*G.parent.u(x0-1+i:xm-1+i,y0-1,:)) + c_y2(1).*(c_x1(i).*G.parent.u(x0-2+i:xm-2+i,y0,:) + c_x2(i).*G.parent.u(x0-1+i:xm-1+i,y0,:)) );
            g_x(i:r:end-r+i,2,:) = g( c_y1(2).*(c_x1(i).*G.parent.u(x0-2+i:xm-2+i,ym,:) + c_x2(i).*G.parent.u(x0-1+i:xm-1+i,ym,:)) + c_y2(2).*(c_x1(i).*G.parent.u(x0-2+i:xm-2+i,ym+1,:) + c_x2(i).*G.parent.u(x0-1+i:xm-1+i,ym+1,:)) );
            g_y(1,i:r:end-r+i,:) = f( c_y1(i).*(c_x1(1).*G.parent.u(x0-1,y0-2+i:ym-2+i,:) + c_x2(1).*G.parent.u(x0,y0-2+i:ym-2+i,:)) + c_y2(i).*(c_x1(1).*G.parent.u(x0-1,y0-1+i:ym-1+i,:) + c_x2(1).*G.parent.u(x0,y0-1+i:ym-1+i,:)) );
            g_y(2,i:r:end-r+i,:) = f( c_y1(i).*(c_x1(2).*G.parent.u(xm,y0-2+i:ym-2+i,:) + c_x2(2).*G.parent.u(xm+1,y0-2+i:ym-2+i,:)) + c_y2(i).*(c_x1(2).*G.parent.u(xm,y0-1+i:ym-1+i,:) + c_x2(2).*G.parent.u(xm+1,y0-1+i:ym-1+i,:)) );
        end
        
%         %...
%         
%         g_x0 = zeros(G.m_x,4);
%         g_xm = zeros(G.m_x,4);
%         g_y0 = zeros(G.m_y,4);
%         g_ym = zeros(G.m_y,4);
%         
%         %...
%         
%         for j = 1:4
%             for i = 0:((G.m_x-1)/r)
% 
%                 g_x0(i*r+1,j) = G.parent.u(x0+i,y0,j);
%                 g_xm(i*r+1,j) = G.parent.u(x0+i,ym,j); % 
%                 
%             end
%             for i = 0:((G.m_y-1)/r)
%                 
%                 g_y0(i*r+1,j) = G.parent.u(x0,y0+i,j);
%                 g_ym(i*r+1,j) = G.parent.u(xm,y0+i,j); %
%                 
%             end
% 
%             i=2;
%             gx0=[g_x0(1,j),g_x0(1+r,j)];
%             gxm=[g_xm(1,j),g_xm(1+r,j)]; % 
% 
%             while i < G.m_x
%                 s = mod(i-1,r);
% 
%                 g_x0(i,j) = ((r-s)*gx0(1) + s*gx0(2))/r;
%                 
%                 g_xm(i,j) = ((r-s)*gxm(1) + s*gxm(2))/r; % 
% 
%                 if (s == r-1) && (i+r+1 <= G.m_x)
%                     gx0 = [gx0(2),g_x0(i+r+1,j)];
%                     
%                     gxm = [gxm(2),g_xm(i+r+1,j)]; % 
%                     
%                     i=i+1;
%                 end
%                 i=i+1;
%             end
%             
%             i=2;
%             gy0=[g_y0(1,j),g_y0(1+r,j)];
%             gym=[g_ym(1,j),g_ym(1+r,j)]; %
%             
%             while i < G.m_y
%                 s = mod(i-1,r);
%                 
%                 g_y0(i,j) = ((r-s)*gy0(1) + s*gy0(2))/r;
%                 
%                 g_ym(i,j) = ((r-s)*gym(1) + s*gym(2))/r; %
% 
%                 if (s == r-1) && (i+r+1 <= G.m_y)
%                     gy0 = [gy0(2),g_y0(i+r+1,j)];
%                     
%                     gym = [gym(2),g_ym(i+r+1,j)]; %
%                     
%                     i=i+1;
%                 end
%                 i=i+1;
%             end
%         end



    end
    
end