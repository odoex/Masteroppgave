m = 100;
n = 300;
x_0 = 0;
x_m = 1;
t_0 = 0;
t_n = 1;

x = linspace(x_0,x_m,m)';
y = linspace(0,1,m)';
t = linspace(t_0,t_n,n)';

k = (t_n-t_0)/(n-1);
h = (x_m-x_0)/(m-1);

[X,Y] = meshgrid(x,y);

a = 0.5;
b = 1;

p = @(s,z) (sin(z-a*s) + sin(-b*s));
g = @(s,z) (sin(-a*s) + sin(z-b*s));

% Difference matrix
v = ones((m-1),1);
D_u = diag(v,1);
D_l = diag(-v,-1);
D = D_u + D_l;
D(1,1) = 0;
D(1,2) = 2;
D(m,m) = 2;
D(m,m-1) = -2;

I = eye(m);

A = kron(a/(2*h)*D,I);
B = kron(I,b/(2*h)*D);

e_0 = zeros(m,1);
e_0(1) = 1;

F = sin(X) + sin(Y);

% Initializing u at time 0
u = reshape(F, [m*m,1]);


%
E=zeros(m);
E(1,1)=1;

C = h*eye(m);
C(1,1) = 1/2;
C(m,m) = 1/2;
Cx=kron(inv(C)*E,-eye(m));
Cy=kron(-eye(m),inv(C)*E);
%
%U = reshape(u,[m,m]);
%disp(U)
for i = 2:n
    P = (b/h)*kron(p(t(i-1),x),e_0);
    G = (a/h)*kron(e_0,g(t(i-1),y));
    time=t(i-1);
    g_funk=g(t(i-1),y);
    u = RK4(u,k,A,B,P,G);
    %if (i==3)
    %    break;
    %end
end

U = reshape(u, [m,m]);
sol = sin(X - a*t(n)) + sin(Y - b*t(n));

disp(U)
%%%%%disp(U);
%%%%%disp(sol);

figure 
mesh(sol);

E = abs(sol-U);
%%%%%disp(E);

error = norm(E);
disp(error);

figure
mesh(U)