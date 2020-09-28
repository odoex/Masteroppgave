% Number of points in space and time
m = 20;
m_x = m;
m_y = m;
n = 2000;

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

g_x = @(s,z) (sin(z-a*s) + sin(-b*s)); 
g_y = @(s,z) (sin(-a*s) + sin(z-b*s));

% Initializing the solution U (skal lagres i grid)
U = F;
U_n = F;

t = t_0;

while t <= t_n
    for i = 1:m_x
        for j = 1:m_y
            gx = g_x(t,x(i));
            gy = g_y(t,y(j));
            
            U_n(i,j) = RK(U,i,j,h,k,a,b,g_x(t,x),g_y(t,y),m_x,m_y);
        end 
    end
    U = U_n;
    t = t + k;
end

% Exact solution
sol = sin(X - a*t_n) + sin(Y - b*t_n);


%disp(U);
%disp(sol);

E = abs(U-sol');
error = norm(E);
disp(error);

% Plot exact solution and approximated solution
figure 
mesh(sol);

figure
mesh(U)

figure 
mesh(E);