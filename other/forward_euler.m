function [u] = forward_euler(u,h,A,B,P,G)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes her

    u = u + h*(A*u + B*u - G - P);
    
end

