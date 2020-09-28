m = 10;
n = 50;
x_0 = 0;
x_m = 1;
t_0 = 0;
t_n = 1;

k = (t_n-t_0)/(n-1);
h = (x_m-x_0)/(m-1);

r = (x_m - x_0)/2;
%r1 = (x_m-x_0)/2;
%r2 = sqrt(2*r^2)-r;
%z = ceil(r2*10)/10;

m2 = (x_m+h*ceil(m/2)) - (x_0-h*ceil(m/2)) / h;
x = linspace(x_0-h*ceil(m/2),x_m+h*ceil(m/2),m2)';
y = linspace(x_0-z,x_m+z,m2)';
t = linspace(t_0,t_n,n)';



[X,Y] = meshgrid(x,y);

X = X - (r);
Y = Y - (r);

X2 = X*cos(-30) - Y*sin(-30);
Y2 = Y*cos(-30) + X*sin(-30);
X = X2 + (r);
Y = Y2 + (r);



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

for i = 2:n
    P = (b/h)*kron(p(t(i-1),x),e_0);
    G = (a/h)*kron(e_0,g(t(i-1),y));
    
    u = RK4(u,k,A,B,P,G);
    U = reshape(u, [m,m]);
    
end

 sol = sin(X - a*t(n)) + sin(Y - b*t(n));


disp(U);
disp(sol);

%figure 
%mesh(sol);

E = abs(sol-U);
disp(E);

%figure
%mesh(U)