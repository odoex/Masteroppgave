function [U] = eRK4(G,t)
% RK_4 function for approximation with Runge-Kutta 4
%   Using Runge-Kutta 4 to calculate the given time step in the method from
%   time t to time t + k. 
    
    k1 = rhs(G.u,t,a,b,f,G);
    k2 = rhs(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
    k3 = rhs(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
    k4 = rhs(G.u + G.k*k3,t + G.k,a,b,f,G);
    

    
    U = G.u + (G.k/6)*(k1 + 2*k2 + 2*k3 + k4);


    
    % Hvordan setter man inn randbetingelser direkte/uten å regne ut flux
    
    % Husk: randbetingelsene er endret til å opdatere seg i tid med runge
    % kutta. Sjekk hvordan dette påvirker de finere gridene og om koden må
    % endres med tanke på deres rand. Der brukes ikke eksakt løsning slik
    % den gjør i hovedgridet. 
    
end
