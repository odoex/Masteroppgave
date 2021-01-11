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

G.u = createSolutionVector(G);

%% Test 1: boundary fine grid

gx = zeros(G.child.m_x,1);
gy = zeros(G.child.m_y,1);

gx(1)=G.u(G.child.location(3),G.child.location(4));
gy(1)=G.u(G.child.location(3),G.child.location(4));

for i = 1:(G.child.m_x+1)/2-1
    gx(2*i+1)=G.u(G.child.location(3)+i,G.child.location(4));
    gx(2*i)=(gx(2*i+1)+gx(2*i-1))/2;
end
for j = 1:(G.child.m_y+1)/2-1
    gy(2*j+1)=G.u(G.child.location(3),G.child.location(4)+j);
    gy(2*j)=(gy(2*j+1)+gy(2*j-1))/2;
end

[g_x,g_y] = boundaryAdv(G.child,0);

assert(max(max(abs(g_x-gx)))+max(max(abs(g_y-gy))) < 1.0e-12)