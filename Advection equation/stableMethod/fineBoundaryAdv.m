function [g_x,g_y] = fineBoundaryAdv(G,f,t) 
    
    x = G.location(3);
    y = G.location(4);

    if G.parent == 0
        g_x = f(linspace(G.location(1),G.location(1)+G.h*(G.m-1),G.m),0,t)'; % Hva gjør jeg med tiden her?
        g_y = f(0,linspace(G.location(2),G.location(2)+G.h*(G.m-1),G.m),t)'; % Hva gjør jeg med tiden her?
    else
        
        g_x = G.g(:,1);
        g_y = G.g(:,2);
        
        
        % r - hvor mange celler er det fine nettet delt inn i 
        r = G.parent.h/G.h; 
        m=(G.m+1)/2;
        p_x = zeros(G.m,1);
        p_y = zeros(G.m,1);
        
        %"Ghost-values" outside the fine grid. 
        %G.g = zeros(G.m,4);

        % Ny plan: lage to vektorer som representerer de to nødvendige
        % punktene utenfor boundarien (kan muligens forbedres senere for å
        % unngå to vektorer) altså det vil si regne de ut i denn metoden. 

        % 1. Dele nettet i to deler og fokusere på det ene interfacet. 

        x = G.location(3); % x coordinate on coarse grid
        y = G.location(4); % y coordinate on coarse grid
        xm = x + (G.m-1)/2;
        ym = y + (G.m-1)/2;
        u = G.u;
        h=G.h*2;
        
        a = 0.5;
        b = 1;

        % First point on the grid interface in y direction: 
        p_y(1) = -a*((G.u(1,1)*(3/4)+G.u(1,2)*(1/4)) - G.parent.u(x-2,y))/(2*h) - b*(G.parent.u(x-1,y+1)-G.parent.u(x-1,y-1))/(2*h);
        u1=G.u(1,1);
        u1b=G.u(1,2);
        u134=(G.u(1,1)*(3/4));
        u1b14=G.u(1,2)*(1/4);
        u1u1b=(G.u(1,1)*(3/4)+G.u(1,2)*(1/4));
        u2=G.parent.u(x-2,y);
        a1=-a*((G.u(1,1)*(3/4)+G.u(1,2)*(1/4)) - G.parent.u(x-2,y))/(2*h);
        b1=- b*(G.parent.u(x-1,y+1)-G.parent.u(x-1,y-1))/(2*h);
        %G.g(1,2) = p_y(1);
        
        % Last point on the grid interface in y direction: -- feil i m?
        % Skal være siste punkt?
        p_y(G.m) = -a*((u(1,G.m)*(1/2)+u(1,G.m-1)*(1/2)) - G.parent.u(x-2,ym))/(2*h) - b*(G.parent.u(x-1,ym)-G.parent.u(x-1,ym-1))/h;
        %G.g(m,2) = p_y(G.m);
        
        % First point on the grid interface in x-direction: 
        p_x(1) = -a*(G.parent.u(x+1,y-1)-G.parent.u(x-1,y-1))/(2*h) - b*((u(1,1)*(3/4)+u(2,1)*(1/4))-G.parent.u(x,y-2))/(2*h);
        %G.g(1,1) = p_x(1);
        
        % Last point on the grid interface in x-direction:
        p_x(G.m) = -a*(G.parent.u(xm,y-1)-G.parent.u(xm-1,y-1))/h - b*((u(G.m-1,1)+u(G.m,1))*(1/2)-G.parent.u(G.m,y-2))/(2*h);
        %G.g(m,1) = p_x(G.m);
        
        
            
        for i = 1:(m-2) % Assuming here that mx=my


            %g_x(i+1,1) = (G.u()) - G.parent.u(x+i,y-2); % x = -1 Ingen
            %interface langs x-aksen enda. 
            p_y(2*i+1) = -a*((u(1,2*i+1)*(1/2)+u(1,2*i+2)*(1/4)+G.u(1,2*i)*(1/4)) - G.parent.u(x-2,y+i))/(2*h) - b*(G.parent.u(x-1,y+i+1)-G.parent.u(x-1,y+i-1))/(2*h); % y = y+2*i/2*i+1, x = -1   G/F
            
            p_x(2*i+1) = -a*(G.parent.u(x+i+1,y-1)-G.parent.u(x+i-1,y-1))/(2*h) - b*((u(2*i+1,1)*(1/2)+u(2*i+2,1)*(1/4)+u(2*i,1)*(1/4))-G.parent.u(x+i,y-2))/(2*h); % y = y+2*i/2*i+1, x = -1   G/F
            
            % Values of even entries in g_x,g_y:
            p_y(2*i) = (1/2)*(p_y(2*i-1)+p_y(2*i+1));
            
            p_y(2*i+2) = (1/2)*(p_y(2*i+1)+p_y(2*i+3));
            
            p_x(2*i) = (1/2)*(p_x(2*i-1)+p_y(2*i+1));
            
            p_x(2*i+2) = (1/2)*(p_x(2*i+1)+p_y(2*i+3));

%             G.g(i+1,1) = p_x(2*i+1);
%             G.g(i+1,2) = p_y(2*i+1);
            
        end
        
        G.g(:,3) = p_x;
        G.g(:,4) = p_y;
        
%         disp(g_x)
%         disp(g_y)
        
    end
end