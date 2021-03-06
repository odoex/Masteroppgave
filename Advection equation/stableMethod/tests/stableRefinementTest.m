% Intervals in time
t_0 = 0;
t_n = 1;
t = t_0;

% Number of points in space and time for main grid
m = 65;
m_x = m;
m_y = m;
n = 2*m*3*t_n; 

% time steps for coarse grid
k = (t_n-t_0)/(n-1);

% PDE coefficients and exact solution
a = 0.5;
b = 1;

f = @(x,y,s) sin(x - a*s) + sin(y - b*s);


% GRID CREATION

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m_x, m_y, n);
G.t=0;
h=G.h;

% Fine grid:
%   Area to be covered by fine grid: Ex locx = [a,b] where a<b
locx = [30,50]; % Correct: this should be decided by location in interval
locy = [30,50];
ratio = 2;

% locx(1) = locx(1)-1; % Adding anothen point to
% locy(1) = locy(1)-1;
% locx(2) = locx(2)+1; 
% locy(2) = locy(2)+1;

location_1 = [(locx(1)-2)*G.h,(locy(1)-2)*G.h,locx(1)-1,locy(1)-1];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +3, (locy(2)-locy(1))*ratio +3, n);
G_1.t = 0;
G.child = G_1;


% SOLUTION VECTORS
G.u = createSolutionVector(G);
if G.child ~= 0
    G.child.u = createSolutionVector(G.child);
    
end

U = G.child.u;

m_x = G.child.m_x;
m_y = G.child.m_y;

r = G.h/G.child.h; % So much to do
    
    x = G.child.location(3); 
    y = G.child.location(4);
    
F_x = zeros(m_x);
F_y = zeros(m_y); 
    
U_f = f_adv(U);
U_g = g_adv(U);

%% Test 1: Points outside grid

h = G.child.h;

Mx = (m_x-3)/r +3;
My = (m_y-3)/r +3;

p_x = zeros(Mx-2,2);
p_y = zeros(2,My-2);

p_xf = [U_f(2:r:m_x-1,1),U_f(2:r:m_x-1,m_y)];
p_xg = [U_g(2:r:m_x-1,1),U_g(2:r:m_x-1,m_y)];
p_yf = [U_f(1,2:r:m_y-1)',U_f(m_y,2:r:m_x-1)']';
p_yg = [U_g(1,2:r:m_y-1)',U_g(m_y,2:r:m_x-1)']';

p_x(1,1) = (p_xf(2,1)-f_adv(G.u(x,y)))/(4*h) + ((3/4)*U_g(2,2)+(1/4)*U_g(3,2) - g_adv(G.u(x+1,y-1)))/(4*h);
p_x(end,1) = (f_adv(G.u(x+Mx-1,y))-p_xf(end-1,1))/(4*h) + ((3/4)*U_g(m_x-1,2)+(1/4)*U_g(m_x-2,2) - g_adv(G.u(x+Mx-1,y-1)))/(4*h); % Feil: skal være - g_adv(G.u(x+Mx-2,y-1)))/(4*h);

p_x(2:end-1,1) = (p_xf(3:end,1) - p_xf(1:end-2,1))/(4*h) + ((1/2)*U_g(4:r:m_x-3,2) + (1/4)*U_g(5:r:m_x-2,2) + (1/4)*U_g(3:r:m_x-4,2) - g_adv(G.u((x+2):(x+Mx-3),y-1)))/(4*h);

p_x(1,2) = ( p_xf(2,2)-f_adv(G.u(x,y+My-1)) )/(4*h) + ( g_adv(G.u(x+1,y+My))-((3/4)*U_g(2,m_y-1)+(1/4)*U_g(3,m_y-1)) )/(4*h);
p_x(end,2) = ( f_adv(G.u(x+Mx-1,y+My-1))-p_xf(end-1,2) )/(4*h) + ( g_adv(G.u(x+Mx-2,y+My))-((3/4)*U_g(m_x-1,m_y-1)+(1/4)*U_g(m_x-2,m_y-1)) )/(4*h);

p_x(2:end-1,2) = (p_xf(3:end,2)-p_xf(1:end-2,2))/(4*h) + (g_adv(G.u(x+2:x+Mx-3,y+My))-((1/2)*U_g(4:r:m_x-3,m_y-1)+(1/4)*U_g(5:r:m_x-2,m_y-1)+(1/4)*U_g(3:r:m_x-4,m_y-1)))/(4*h);

p_y(1,1) = ((3/4)*U_f(2,2)+(1/4)*U_f(2,3) - f_adv(G.u(x-1,y+1)))/(4*h) + (p_yg(1,2)-g_adv(G.u(x,y)))/(4*h);
p_y(1,end) = ((3/4)*U_f(2,m_y-1)+(1/4)*U_f(2,m_y-2) - f_adv(G.u(x-1,y+My-2)))/(4*h) + (g_adv(G.u(x,y+My-1)) - p_yg(1,end-1))/(4*h);

p_y(1,2:end-1) = ((1/2)*U_f(2,4:r:m_y-3)+(1/4)*U_f(2,5:r:m_y-2)+(1/4)*U_f(2,3:r:m_y-4) - f_adv(G.u(x-1,y+2:y+My-3)))/(4*h) + (p_yg(1,3:end)-p_yg(1,1:end-2))/(4*h);

p_y(2,1) = (f_adv(G.u(x+Mx,y+1)) - ((3/4)*U_f(m_x-1,2)+(1/4)*U_f(m_x-1,3)))/(4*h) + (p_yg(2,2)-g_adv(G.u(x+Mx-1,y)))/(4*h);
p_y(2,end) = (f_adv(G.u(x+Mx,y+My-2)) - ((3/4)*U_f(m_x-1,m_y-1)+(1/4)*U_f(m_x-1,m_y-2)))/(4*h) + (g_adv(G.u(x+Mx-1,y+My-1)) - p_yg(2,end-1))/(4*h);

p_y(2,2:end-1) = (f_adv(G.u(x+Mx,y+2:y+My-3)) - ((1/2)*U_f(m_x-1,4:r:m_y-3)+(1/4)*U_f(m_x-1,5:r:m_y-2)+(1/4)*U_f(m_x-1,3:r:m_y-4)))/(4*h) + (p_yg(2,3:end)-p_yg(2,1:end-2))/(4*h);

U_n = rhsFineGridAdv(G.child.u,0,G.child);

assert(max(max(abs(U_n(2:r:end-1,1)-p_x(:,1)))) + max(max(abs(U_n(2:r:end-1,m_y)-p_x(:,2)))) + max(max(abs(U_n(1,2:r:end-1)-p_y(1,:)))) + max(max(abs(U_n(m_y,2:r:end-1)-p_y(2,:)))) < 1.0e-12)

%% Test 2: Boundary points

U = createSolutionVector(G.child);
U2 = U;
v_x = [2,2,m_x-1,m_x-1]';
v_y = [2,m_y-1]';

U2(2,v_y) = (U_f(2+1,v_y)-U_f(2-1,v_y))/(3*h) + (U_g(2,v_y+1)-U_g(2,v_y-1))/(3*h);
U2(m_x-1,v_y) = (U_f(m_x-1+1,v_y)-U_f(m_x-1-1,v_y))/(3*h) + (U_g(m_x-1,v_y+1)-U_g(m_x-1,v_y-1))/(3*h);

% U2(2,2) = (U_f(3,2)-U_f(1,2))/(3*h) + (U_f(2,3)-U_f(2,1))/(3*h);
% U2(2,m_y-1) = (U_f(3,2)-U_f(1,2))/(3*h) + (U_f(2,3)-U_f(2,1))/(3*h);
% U2(m_x-1,2) = 
% U2(m_x-1,m_y-1) = 

U2(3:r:m_x-2,2) = (U_f(4:r:m_x-1,2)-U_f(2:r:m_x-3,2))/(2*h) + (U_g(3:r:m_x-2,3) - ((1/2)*U_g(4:r:m_x-1,1)+(1/2)*U_g(2:r:m_x-3,1)))/(3*h);
U2(4:r:m_x-3,2) = (U_f(5:r:m_x-2,2)-U_f(3:r:m_x-4,2))/(2*h) + (U_g(4:r:m_x-3,3) - (U_g(4:r:m_x-3,1)))/(3*h);

U2(3:r:m_x-2,m_y-1) = (U_f(4:r:m_x-1,m_y-1)-U_f(2:r:m_x-3,m_y-1))/(2*h) + (-U_g(3:r:m_x-2,m_y-2) + ((1/2)*U_g(4:r:m_x-1,m_y)+(1/2)*U_g(2:r:m_x-3,m_y)))/(3*h);
U2(4:r:m_x-3,m_y-1) = (U_f(5:r:m_x-2,m_y-1)-U_f(3:r:m_x-4,m_y-1))/(2*h) + (-U_g(4:r:m_x-3,m_y-2) + (U_g(4:r:m_x-3,m_y)))/(3*h);

U2(2,3:r:m_x-2) = (U_g(2,4:r:m_x-1)-U_g(2,2:r:m_x-3))/(2*h) + (U_f(3,3:r:m_x-2) - ((1/2)*U_f(1,4:r:m_x-1)+(1/2)*U_f(1,2:r:m_x-3)))/(3*h);
U2(2,4:r:m_x-3) = (U_g(2,5:r:m_x-2)-U_g(2,3:r:m_x-4))/(2*h) + (U_f(3,4:r:m_x-3) - (U_f(1,4:r:m_x-3)))/(3*h);

U2(m_x-1,3:r:m_y-2) = (U_g(m_x-1,4:r:m_y-1)-U_g(m_x-1,2:r:m_y-3))/(2*h) + (-U_f(m_x-2,3:r:m_y-2) + ((1/2)*U_f(m_x,4:r:m_y-1)+(1/2)*U_f(m_x,2:r:m_y-3)))/(3*h);
U2(m_x-1,4:r:m_y-3) = (U_g(m_x-1,5:r:m_y-2)-U_g(m_x-1,3:r:m_y-4))/(2*h) + (-U_f(m_x-2,4:r:m_y-3) + (U_f(m_x,4:r:m_y-3)))/(3*h);

U_n = rhsFineGridAdv(G.child.u,0,G.child);

assert(max(max(abs(U_n(2:end-1,2)-U2(2:end-1,2)))) + max(max(abs(U_n(2:end-1,m_y-1)-U2(2:end-1,m_y-1)))) + max(max(abs(U_n(2,2:end-1)-U2(2,2:end-1)))) + max(max(abs(U_n(m_x-1,2:end-1)-U2(m_x-1,2:end-1)))) < 1.0e-12)

%% Test 3: inner points

U = createSolutionVector(G.child);
U3 = U;

U3(3:m_x-2,3:m_y-2) = (U_f(4:m_x-1,3:m_y-2)-U_f(2:m_x-3,3:m_y-2))/(2*h) + (U_g(3:m_x-2,4:m_y-1)-U_g(3:m_x-2,2:m_y-3))/(2*h);

U_n = rhsFineGridAdv(G.child.u,0,G.child);

assert(max(max(abs(U3(3:m_x-2,3:m_y-2)-U_n(3:m_x-2,3:m_y-2)))) < 1.0e-12)
