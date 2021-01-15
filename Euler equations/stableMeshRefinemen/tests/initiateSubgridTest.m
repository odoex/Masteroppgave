% Interval in space 
x = [0,1];
y = [0,1];

m = 100; 
m_x = m;
m_y = m;

% Space step
h = (x(2)-x(1))/(m-1);

% Interval in time
t_0 = 0;
t_n = 0.15;  
t = t_0;

% Mesh refinement 
ratio = 2;
locx = [0.3,0.7];
locy = [0.3,0.7]; % Standard: [0.35,0.65];

% Solution vector
x = linspace(x(1),x(2),m_x)';
y = linspace(y(1),y(2),m_y)';

u = exactSolEuler(x,y,0);

% Points in time and time step
c = 331; % speed of sound

maxSpeed = max(max(u(:,:,2)));

k = h/(ratio*(c + maxSpeed)); 

n = floor((t_n-t_0)/k)+1;


% GRID CREATION

% Coarse grid:
G = Node(0, [x(1),x(1),1,1], h, k, m_x, m_y, n);
G.t=t;
G.u = u;

% Fine grid: 

locx = [round(locx(1)/h)+1,round(locx(2)/h)+1];
locy = [round(locy(1)/h)+1,round(locy(2)/h)+1];
location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, (locy(2)-locy(1))*ratio +1, n);
G_1.t = 0;

%% Test 1: initiate subgrid test

r = ratio;
m_x = G_1.m_x;
m_y = G_1.m_y;
x = G_1.location(3); 
y = G_1.location(4);
Mx = x+(m_x-1)/r; % The extra points are not added yet
My = y+(m_y-1)/r;

G.child = initiateSubgrid(G_1,ratio);
    
A = zeros(Mx+1-x,My+1-y);
for k = 1:4
    A(2:end-1,2:end-1,k) = G.child.u(2:r:m_x-1,2:r:m_y-1,k)-G.u(x:Mx,y:My,k);
    
    A(2:end-1,1,k) = G.child.u(2:r:m_x-1,1,k)-G.u(x+1:Mx-1,y,k);
    A(2:end-1,end,k) = G.child.u(2:r:m_x-1,m_y,k)-G.u(x+1:Mx-1,My,k);
    A(1,2:end-1,k) = G.child.u(1,2:r:m_y-1,k)-G.u(x,y+1:My-1,k);
    A(end,2:end-1,k) = G.child.u(m_x,2:r:m_y-1,k)-G.u(Mx,y+1:My-1,k);
end

assert(max(max(max(abs(A)))) < 1.0e-12)
