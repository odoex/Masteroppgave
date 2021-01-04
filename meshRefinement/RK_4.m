function [U] = RK_4(G,t)
% RK_4 function for approximation with Runge-Kutta 4
%   Using Runge-Kutta 4 to calculate the given time step in the method from
%   time t to time t + k. 
    
    k1 = rhs(G.u,t,G);
    k2 = rhs(G.u + (G.k/2)*k1,t + G.k/2,G);
    k3 = rhs(G.u + (G.k/2)*k2,t + G.k/2,G);
    k4 = rhs(G.u + G.k*k3,t + G.k,G);
    
%     [Y,X] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
%         figure
%         mesh(X,Y,k1(:,:,2))
    
    U = G.u + (G.k/6)*(k1 + 2*k2 + 2*k3 + k4);


    
    % Hvordan setter man inn randbetingelser direkte/uten å regne ut flux
    
    % Husk: randbetingelsene er endret til å opdatere seg i tid med runge
    % kutta. Sjekk hvordan dette påvirker de finere gridene og om koden må
    % endres med tanke på deres rand. Der brukes ikke eksakt løsning slik
    % den gjør i hovedgridet. 
    
end
