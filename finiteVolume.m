function [G] = finiteVolume(G,a,b,g_x,g_y,t_0,t_n,f)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % G.t - tidssteget som sist ble regnet ut
    
    t = t_0;  % t = 0, 0
    
    

    while t <= (t_n-G.k) % t_n = k*n, 1, 1/2
        u = RK_4(G,t,a,b,f);
        %u = RK_4(G.u,t,G.h,G.k,a,b,g_x,g_y,G.m,G.m);
        %G.u = RK_4(G.u,t,G.h,G.k,a,b,g_x,g_y,G.m,G.m); % regner ut ved tid t+k, 1, 1/2, 1/4, 1/2
        G.t = t + G.k; % Kan puttes inn i en annen funksjon?
        
        if (G.child ~= 0)
            G = finiteVolume(G.chid,a,b,g_x,g_y,t_0,t+k); % tidssteget der det slapp
            G = projectFlux(G,G.child);
        end
        
        G.u = u;
        t = t + G.k; % 1/4, 1/2
    end
    
end

