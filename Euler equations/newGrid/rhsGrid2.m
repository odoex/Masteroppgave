function [u] = rhsGrid2(u,t,G)
% rhs calculates the right side of the equation
%   Collecting the boundary vectors for the grid, the function loops
%   through the grid and approximates the flux. Checking if the point is at
%   the boundary, the function uses the boundary conditions to calculate
%   the flux when at the boundary where boundary condition is needed, and
%   adapts the method at the other boundaries (i = mx or j = my). The
%   solution is inserted to a temporary solution array U_n and is returned
%   in U after the loop is done. 

    h = G.h;
    m_x = G.m_x;
    m_y = G.m_y;
    
    delta = 0;%6.5;
    
    variables = length(u(1,1,:));
    
    F_x = zeros(m_x,m_y,variables);
    F_y = zeros(m_x,m_y,variables);
    
    [p_x,p_y] = boundaryGrid2(G,t);
    
    f_px = g(p_x);
    f_py = f(p_y);
    
    U_f = f(u);
    U_g = g(u);

	F_x(1,:,:) = - ((U_f(2,:,:) - f_py(1,:,:)) - delta*(u(2,:,:)-u(1,:,:)))/h;
	F_y(:,1,:) = - ((U_g(:,2,:) - f_px(:,1,:)) - delta*(u(:,2,:)-u(:,1,:)))/h;

	F_x(m_x,:,:) = - ((f_py(2,:,:) - U_f(m_x-1,:,:)) + delta*(u(m_x,:,:)-u(m_x-1,:,:)))/h;
    F_y(:,m_y,:) = - (f_px(:,2,:) - U_g(:,m_y-1,:) + delta*(u(:,m_y,:)-u(:,m_y-1,:)))/h;

    F_x(2:end-1,:,:) = - (U_f(3:end,:,:) - U_f(1:end-2,:,:) - delta*(u(3:end,:,:)-u(2:end-1,:,:)) + delta*(u(2:end-1,:,:)-u(1:end-2,:,:)))/(2*h);
    F_y(:,2:end-1,:) = - (U_g(:,3:end,:) - U_g(:,1:end-2,:) - delta*(u(:,3:end,:)-u(:,2:end-1,:)) + delta*(u(:,2:end-1,:)-u(:,1:end-2,:)))/(2*h);

	u_n = F_x + F_y;
    
    u = boundaryValues2(G,u_n,u,t,delta); % Kan putte alle bv inn her?

end

