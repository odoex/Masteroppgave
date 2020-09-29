% Number of points in space and time
m = 3;
n = 10;

y_0 = 0;
y_m = 1;
x_0 = 0;
x_m = 1;
t_0 = 0;
t_n = 1;

% Space and time steps
h = (x_m-x_0)/(m-1);
k = (t_n-t_0)/(n-1);

% Grid creation (Ikke lagre i skriptet ovenfor)
% G = Node(0, [x_0,y_0], h, k, m, n);
% x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
% y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';

x = linspace(x_0,x_m,m)';
y = linspace(y_0,y_m,m)';

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

f = @(x,y,s) sin(x - a*s) + sin(y - b*s);

G = Node(0,[0,0,1,1],1/(m-1), k, m, n);

G.u = createSolutionVector(G);

round = @(v,d) floor(v*10^d)/10^d;

%% Test 1: createSolutionVector
val=floor(norm(G.u - F)*1000)/1000;
assert(floor(norm(G.u - F)*1000)/1000 == 0);

%% Test 2: boundary
U = G.u;
[gx,gy] = boundary(G,f,k);
v = g_x(k,x);
assert(round(norm(g_x(k,x)-gx),3)==0)