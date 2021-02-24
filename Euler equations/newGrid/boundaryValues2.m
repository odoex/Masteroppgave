function [u_n] = boundaryValues2(G,u_n,u,t,delta)
% boundaryValues2 calculates boundary values at a boundary
%   This function takes a grid node, a solution vector u_n that is updated
%   by finite volumes, and a solution vector u which is the solution that
%   is to be updated to the solution of the next step, u_n. It also takes
%   the given time and the diffusion constant delta. The function finds the
%   exact solutions at the boundary written below, and calculates the
%   closest points to the boundary on the grid where these points are
%   imposed. 

    % Have to assume at first that at the boundary y < 1 when x = 0, and it
    % is always in the lower left corner of the domain, that is, the lower
    % values of x and y are outside the boundary. Also that x < 1 when y=0.
    
    % If the absolute value of the slope of the boundary is less than 1,
    % then we must start with the x-values and find an appropriate y-value.
    % Otherwise we must do the oposite and start with the x-values. Let us
    % therefore always assume that the slope is less than 1, and test our
    % case this way. 
    
    % In this function there is two seperate function for each coordinate,
    % x and y, and x_b and y_b. The first two are the indeces for for the
    % points on the grid, and the second two are the exact points on the
    % boundary. 
    
    x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
%     y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
    
    % Boundary 
    a = -0.5;
    b = 0.4;
    eqn_y = @(var_x) a*var_x+b;
        
    y_b = eqn_y(x); % y values at the boundary corresponding to the x values
    
    y_b = y_b(y_b>=0);
    x_b = x(1:length(y_b));
    
	x = round(x_b/G.h)+1; % index for x
    y = round(y_b/G.h)+1; % index for y
    
    u_f = f(u);
    u_g = g(u);

	for i = 1:length(x_b)
        if (x(i) == 1 || i == 1 || u_n(x(i-1),y(i),1) == 0)
            u_n(x(i),y(i),:) = -(u_f(x(i)+1,y(i),:) - f(exactSolEuler(x_b(i),y_b(i),t))- delta*(u(x(i)+1,y(i),:)-u(x(i),y(i),:)))/G.h;
        else 
            u_n(x(i),y(i),:) = -(u_f(x(i)+1,y(i),:) - u_f(x(i)-1,y(i),:)- delta*(u(x(i)+1,y(i),:)-2*u(x(i),y(i),:) + u(x(i)-1,y(i),:)))/(2*G.h);
        end
        
%         if (y(i) == 1 || i == 1 || y(i-1) == 0)
            u_n(x(i),y(i),:) = u_n(x(i),y(i),:) -(u_g(x(i),y(i)+1,:) - g(exactSolEuler(x_b(i),y_b(i),t)) - delta*(u(x(i),y(i)+1,:)-u(x(i),y(i),:)))/G.h;
%         else
%             u_n(x(i),y(i),:) = u_n(x(i),y(i),:) -(u_g(x(i),y(i)+1,:) - u_g(x(i),y(i)-1,:) - delta*(u(x(i),y(i)+1,:)-2*u(x(i),y(i),:)+u(x(i),y(i)-1,:)))/(2*G.h);
%         end
        
        u_n(x(i),1:y(i)-1,:) = zeros(1,length(1:y(i)-1),4);
	end

end