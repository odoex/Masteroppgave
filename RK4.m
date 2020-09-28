function [u] = RK4(u,k,A,B,P,G)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    k1 = -(A*u + B*u - G - P);
    k2 = -(A*(u + (k/2)*k1) + B*(u + (k/2)*k1) - G - P);
    k3 = -(A*(u + (k/2)*k2) + B*(u + (k/2)*k2) - G - P);
    k4 = -(A*(u + k*k3) + B*(u + k*k3) - G - P);
    
    u = u + (k/6)*(k1 + 2*k2 + 2*k3 + k4);

end

