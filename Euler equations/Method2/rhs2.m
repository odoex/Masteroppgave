function [u] = rhs2(u,t,G)
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
    
    delta = 2;

    F_x = zeros(m_x,m_y);
    F_y = zeros(m_x,m_y);
    
    [p_x,p_y] = boundary2(G,t);
    
%     px = zeros(2,length(p_x0),length(u(1,1,:)));
%     py = zeros(2,length(p_y0),length(u(1,1,:)));
%     px(1,:,:) = p_x0;
%     px(2,:,:) = p_xm;
%     py(1,:,:) = p_y0; % Fiks dette rotet
%     py(2,:,:) = p_ym; % Fiks dette rotet
%     
%     f_px = g(px);
%     f_py = f(py);
    
    U_f = f(u);
    U_g = g(u);

for v = 1:4
    F_x(1,:) = - ((U_f(2,:,v) - p_y(1,:,v)) - delta*(u(2,:,v)-u(1,:,v)))/h;
    F_y(:,1) = - ((U_g(:,2,v) - p_x(:,1,v)) - delta*(u(:,2,v)-u(:,1,v)))/h;
    
    F_x(m_x,:) = - (p_y(2,:,v) - U_f(m_x-1,:,v) + delta*(u(m_x,:,v)-u(m_x-1,:,v)))/h;
    F_y(:,m_y) = - (p_x(:,2,v) - U_g(:,m_y-1,v) + delta*(u(:,m_y,v)-u(:,m_y-1,v)))/h;
    
    F_x(2:end-1,:) = - (U_f(3:end,:,v) - U_f(1:end-2,:,v) - delta*(u(3:end,:,v)-u(2:end-1,:,v)) + delta*(u(2:end-1,:,v)-u(1:end-2,:,v)))/(2*h);
    F_y(:,2:end-1) = - (U_g(:,3:end,v) - U_g(:,1:end-2,v) - delta*(u(:,3:end,v)-u(:,2:end-1,v)) + delta*(u(:,2:end-1,v)-u(:,1:end-2,v)))/(2*h);
    
    u(:,:,v) = F_x + F_y;
end

end