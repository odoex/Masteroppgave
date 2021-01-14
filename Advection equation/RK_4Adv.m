function [U] = RK_4Adv(G,t,a,b,f)
% RK_4 function for approximation with Runge-Kutta 4
%   Using Runge-Kutta 4 to calculate the given time step in the method from
%   time t to time t + k. 
    
%     if(G.parent ~= 0) 
%         k1 = rhsFineGridAdv(G.u,t,a,b,f,G);
%         k1g1 = G.g(:,3);
%         k1g2 = G.g(:,4);
%         k2 = rhsFineGridAdv(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%         k2g1 = G.g(:,3) + (G.k/2)*k1g1;
%         k2g2 = G.g(:,4) + (G.k/2)*k1g2;
%         k3 = rhsFineGridAdv(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%         k3g1 = G.g(:,3) + (G.k/2)*k2g1;
%         k3g2 = G.g(:,4) + (G.k/2)*k2g2;
%         k4 = rhsFineGridAdv(G.u + G.k*k3,t + G.k,a,b,f,G);
%         k4g1 = G.g(:,3) + (G.k)*k3g1;
%         k4g2 = G.g(:,4) + (G.k)*k3g2;
%         
%         G.g(:,3) = G.g(:,1) + (G.k/6)*(k1g1 + 2*k2g1 + 2*k3g1 + k4g1);
%         G.g(:,4) = G.g(:,2) + (G.k/6)*(k1g2 + 2*k2g2 + 2*k3g2 + k4g2);
%     else
%         k1 = rhsFineGridAdv(G.u,t,a,b,f,G);
%         k2 = rhsFineGridAdv(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%         k3 = rhsFineGridAdv(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%         k4 = rhsFineGridAdv(G.u + G.k*k3,t + G.k,a,b,f,G);
%     end
    
    
%     if(G.parent ~= 0) 
%         k1 = rhsFineGridAdv2(G.u,t,a,b,f,G);
%         k2 = rhsFineGridAdv2(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%         k3 = rhsFineGridAdv2(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%         k4 = rhsFineGridAdv2(G.u + G.k*k3,t + G.k,a,b,f,G);
%     else
%         k1 = rhsAdv(G.u,t,a,b,f,G);
%         k2 = rhsAdv(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%         k3 = rhsAdv(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%         k4 = rhsAdv(G.u + G.k*k3,t + G.k,a,b,f,G);
%     end
        
%     if(G.parent ~= 0) 
%         k1 = rhsFineGridAdv2(G.u,t,a,b,f,G);
%         k2 = rhsFineGridAdv2(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%         k3 = rhsFineGridAdv2(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%         k4 = rhsFineGridAdv2(G.u + G.k*k3,t + G.k,a,b,f,G);
%     else
%         k1 = rhsCoarseGridAdv(G.u,t,a,b,f,G);
%         k2 = rhsCoarseGridAdv(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%         k3 = rhsCoarseGridAdv(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%         k4 = rhsCoarseGridAdv(G.u + G.k*k3,t + G.k,a,b,f,G);
%     end
    
%     k1 = rhsAdv(G.u,t,a,b,f,G);
%     k2 = rhsAdv(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
%     k3 = rhsAdv(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
%     k4 = rhsAdv(G.u + G.k*k3,t + G.k,a,b,f,G);
%      
%     
%     U = G.u + (G.k/6)*(k1 + 2*k2 + 2*k3 + k4);

    if(G.parent ~= 0) 
        k1 = rhsFineGridAdv(G.u,t,G);
        k2 = rhsFineGridAdv(G.u + (G.k/2)*k1,t + G.k/2,G);
        k3 = rhsFineGridAdv(G.u + (G.k/2)*k2,t + G.k/2,G);
        k4 = rhsFineGridAdv(G.u + G.k*k3,t + G.k,G);
    else
        k1 = rhsAdv(G.u,t,a,b,f,G);
        k2 = rhsAdv(G.u + (G.k/2)*k1,t + G.k/2,a,b,f,G);
        k3 = rhsAdv(G.u + (G.k/2)*k2,t + G.k/2,a,b,f,G);
        k4 = rhsAdv(G.u + G.k*k3,t + G.k,a,b,f,G);
    end
    
    U = G.u + (G.k/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % Husk: randbetingelsene er endret til å opdatere seg i tid med runge
    % kutta. Sjekk hvordan dette påvirker de finere gridene og om koden må
    % endres med tanke på deres rand. Der brukes ikke eksakt løsning slik
    % den gjør i hovedgridet. 
    
end

