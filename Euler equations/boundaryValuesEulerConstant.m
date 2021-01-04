function [b] = boundaryValuesEulerConstant(x,y,t)
% boundaryValues returns boundary values for the domain
%   Takes the discrete x and y values, and calculates the boundary
%   at the given time for the euler equations

    if x == 0
        v = y;
    else 
        v = x;
    end
    gamma = 1.4;
    b = zeros(length(v),4);
    
    % Negative since on the left side of the equation
    b(:,1) = 1;
    b(:,2) = b(:,1)*2;
    b(:,3) = b(:,1)*3;
    b(:,4) = 4/(gamma-1)+(1/2)*(b(:,2).^2+b(:,3).^2)./b(:,1);
    
end