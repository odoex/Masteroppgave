% Number of points in space and time
m = 9;
m_x = m;
m_y = m;
n = 500;

% Intervals in space and time
x_0 = 0;
x_m = 1;
y_0 = 0;
y_m = 1;
t_0 = 0;
t_n = 1;

% Space and time steps
h = (x_m-x_0)/(m-1);
k = (t_n-t_0)/(n-1);

% Grid creation (Ikke lagre i skriptet ovenfor)
G = Node(0, [x_0,y_0], h, k, m, n);
x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';

[X,Y] = meshgrid(x,y);

a = 0.5;
b = 1;

% Initial and boundary conditions
F = sin(X) + sin(Y);
ex_sol = @(s) sin(X-a*s) + sin(Y-b*s);

g_x = @(s,i) (sin(x(i)-a*s) + sin(-b*s)); 
g_y = @(s,j) (sin(-a*s) + sin(y(j)-b*s));

% Initializing the solution U (skal lagres i grid)
U = F;

t = linspace(t_0,t_n,n);


for i = 1:(n-1)
    
    U = RK_ny(U,t(i),h,k,a,b,g_x,g_y,m_x,m_y);

     
    %U(1,:) = sin(-a*t(i+1)) + sin(y-b*t(i+1));
    %U(:,1) = sin(x-a*t(i+1)) + sin(-b*t(i+1));
    
%      disp(U)
%      disp(ex_sol(t(i+1))')
    
end

% Exact solution
sol = sin(X - a*t_n) + sin(Y - b*t_n);

% Error estimation
E = abs(U-sol');

error = norm(E);
disp(error)

% error = 0;
% for i=1:m_x
%     for j=1:m_y
%         error = error + h*(U(i,j) - sol(i,j))^2;
%     end
% end
% 
% error = sqrt(error);
% 
% disp(error);


% Plot exact solution and approximated solution
figure 
mesh(sol);

figure
mesh(U)

figure 
mesh(E);