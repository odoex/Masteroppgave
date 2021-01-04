function [u_ij] = RK(u,i,j,h,k,a,b,g_x,g_y,m_x,m_y)
% Runge Kutta 4
%   Runs Runge Kutta for U(i,j) with space step h, time step k, where 
%   g_x and g_y are boundary functions and flux() is the right hand 
%   of the equation.
    
    U = @(x,y) u(x,y);
    
    k1 = @(i1,j1) flux(i1,j1,U,a,b,g_x,g_y,m_x,m_y,h);
    k2 = @(i2,j2) flux(i2,j2, @(e2,f2) U(e2,f2) + k/2*k1(e2,f2),a,b,g_x,g_y,m_x,m_y,h);
    k3 = @(i3,j3) flux(i3,j3, @(e3,f3) U(e3,f3) + k/2*k2(e3,f3),a,b,g_x,g_y,m_x,m_y,h);
    k4 = @(i4,j4) flux(i4,j4, @(e4,f4) U(e4,f4) + k*k3(e4,f4),a,b,g_x,g_y,m_x,m_y,h);
    
    u_ij = u(i,j) + (k/6)*(k1(i,j) + 2*k2(i,j) + 2*k3(i,j) + k4(i,j));
    
end

