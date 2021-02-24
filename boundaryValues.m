function [u_n] = boundaryValues(G,u_n,u,t,delta)

    % Have to assume at first that at the boundary y < 1 when x = 0, and it
    % is always in the lower left corner of the domain, that is, the lower
    % values of x and y are outside the boundary. Also that x < 1 when y=0.
    
    % If the absolute value of the slope of the boundary is less than 1,
    % then we must start with the x-values and find an appropriate y-value.
    % Otherwise we must do the oposite and start with the x-values. Let us
    % therefore always assume that the slope is less than 1, and test our
    % case this way. 
    
    % !Assuming that there is no case where y = 0 will not be included in
    % the staircase boundary. If it is, the calculations will be wrong as
    % the rest of the boundary, lying at y = 0, will be missing a half
    % stencil at the intersection point with the staircase boundary. 
    
    x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
    y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
    
    % Boundary 
    a = -1;%-0.5;
    b = 0.5;
    eqn_y = @(var_x) a*var_x+b;
	
    y_b = eqn_y(x);
    
    x_b = round(x/G.h)+1; % index for x
    y_b = round(y_b/G.h+1.000000000000000e-05)+1; % index for y
    
    y_b = y_b(y_b>=1);
    x_b = x_b(1:length(y_b));

    u_exactSol = exactSolEuler(x,y,t);
%     u_exactSol = initialConditionsEulerConstant(x,y,t);
    u_exactSol_f = f(u_exactSol);
    u_f = f(u);
    
    u_exactSol_g = g(u_exactSol);
    u_g = g(u);

% Pass på at ikke randen kommer helt bort til x=1 randen
% kan bruke i som x-verdi??
	for i = 1:length(y_b)
        if (x_b(i) == 1 || i == 1 || u_n(x_b(i-1),y_b(i),1) == 0) % Unødvendig mange
            u_n(x_b(i),y_b(i),:) = -(u_f(x_b(i)+1,y_b(i),:) - u_exactSol_f(x_b(i),y_b(i),:)- delta*(u(x_b(i)+1,y_b(i),:)-u(x_b(i),y_b(i),:)))/G.h;
        elseif (x_b(i) == G.m_x)
            u_n(x_b(i),y_b(i),:) = -(u_exactSol_f(x_b(i),y_b(i),:) - u_f(x_b(i)-1,y_b(i),:)+ delta*(u(x_b(i),y_b(i),:) - u(x_b(i)-1,y_b(i),:)))/(G.h);
        else
            u_n(x_b(i),y_b(i),:) = -(u_f(x_b(i)+1,y_b(i),:) - u_f(x_b(i)-1,y_b(i),:)- delta*(u(x_b(i)+1,y_b(i),:)-2*u(x_b(i),y_b(i),:) + u(x_b(i)-1,y_b(i),:)))/(2*G.h);
        end
        
        u_n(x_b(i),y_b(i),:) = u_n(x_b(i),y_b(i),:) -(u_g(x_b(i),y_b(i)+1,:) - u_exactSol_g(x_b(i),y_b(i),:) - delta*(u(x_b(i),y_b(i)+1,:)-u(x_b(i),y_b(i),:)))/G.h;

        u_n(x_b(i),1:y_b(i)-1,:) = zeros(1,length(1:y_b(i)-1),4);
	end

end