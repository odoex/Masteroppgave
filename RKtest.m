% Number of points in space and time
m = 3;
n = 10;

x_0 = 0;
x_m = 1;
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

round = @(v,d) floor(v*10^d)/10^d;

%% Test 1: initial conditions
g = sin(x(2)-a*0) + sin(y(2)-b*0);
v1=floor(g*1000)/1000;
v2=floor(F(2,2)*1000)/1000;
assert(floor(g*1000)/1000 == floor(F(2,2)*1000)/1000)

%% Test 2: flux function
rhs = flux(1,1,U,a,b,@(i) g_x(k,x(i)),@(j) g_y(k,y(j)),m_x,m_y,h);
w1=round(rhs,3);
w2=round(-(a*U(2,1) - a*U(1,1))/h - (b*U(1,2)-b*U(1,1))/h,3);
assert(round(rhs,3) == round(-(a*U(2,1) - a*U(1,1))/h - (b*U(1,2)-b*U(1,1))/h,3))

