x = [0,1];
y = [0,1];
m_x = 50;
m_y = 50; 

h = (x(2)-x(1))/(m_x-1);

x = linspace(x(1),x(2),m_x);
y = linspace(y(1),y(2),m_y);

t_int = [0,1];
t = 0; 

n = round(1/(h/2));
k = 1/n;

G = Node(0, [x(1),y(1),1,1], h, k, m_x, m_y, n);
G.t=0;

% Refinement
ratio = 2;
locx = [0.3,0.7];
locy = [0.3,0.7]; 

locx = [(round(locx(1)/h))+1, (round(locx(2)/h))+1];
locy = [(round(locy(1)/h))+1, (round(locy(2)/h))+1];

G.child = Node(G, [(locx(1)-1)*h,(locy(1)-1)*h,locx(1),locy(1)], h/ratio, k, ratio*(locx(2)-locx(1))+1, ratio*(locy(2)-locy(1))+1,n);
G.child.t = 0;

%% Test 1: initial conditions coarse grid
[X,Y] = meshgrid(x,y,G.h);
sol = sin(X) + sin(Y);

G.u = createSolutionVector(G);

assert(max(max(abs(G.u-sol)))<1.0e-12)

%% Test 2: initial conditions fine grid
x_1 = linspace((locx(1)-1)*h,(locx(2)-1)*h,G.child.m_x);
y_1 = linspace((locy(1)-1)*h,(locy(2)-1)*h,G.child.m_x);
[X_1,Y_1] = meshgrid(x_1,y_1,G.child.h);

sol_1 = sin(X_1) + sin(Y_1);

G.child.u = createSolutionVector(G.child);

assert(max(max(abs(G.child.u-sol_1)))<1.0e-12)

%% Test 3: same initial conditions on coarse and fine grid

G.child.u = createSolutionVector(G.child);
G.u = createSolutionVector(G);

A = zeros(G.child.m_x);

for i = 1:((G.child.m_x+1)/2)
    for j = 1:((G.child.m_y+1)/2)
        A(i,j)=G.child.u(2*i-1,2*j-1)-G.u(G.child.location(3)-1+i,G.child.location(4)-1+j);
    end
end

assert(max(max(abs(A))) <1.0e-12)
        