format long

% Interval in space 
x = [0,1];
y = [0,1];

m = 50; 
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
%   Area to be covered by fine grid: Ex locx = [a,b] where a<b
if (ratio > 1)
    locx = [round(locx(1)/h)+1,round(locx(2)/h)+1];
    locy = [round(locy(1)/h)+1,round(locy(2)/h)+1];
    location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

    G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, (locy(2)-locy(1))*ratio +1, n);
    G_1.t = 0;

    G.child = initiateSubgrid(G_1,ratio);
end

%% Test 1: mesh refinement on points outside grid f flux

h = G.child.h;
U = G.child.u;
r = ratio;

U_f = f(U);
U_g = g(U);
Up_f = f(G.u);

m_x = G.child.m_x;
m_y = G.child.m_y;
x = G.child.location(3); 
y = G.child.location(4);
Mx = (m_x-3)/r +3;
My = (m_y-3)/r +3;

p_row_l = zeros(Mx-2,2,4);
p_row_r = zeros(Mx-2,2,4);
p_col_l = zeros(2,My-2,4);
p_col_r = zeros(2,My-2,4);

p_row = [U_f(2:r:m_x-1,1,:),U_f(2:r:m_x-1,m_y,:)];
p_col = [U_f(1,2:r:m_y-1,:);U_f(m_y,2:r:m_x-1,:)];

p_row_l(1,1,:) = Up_f(x,y,:);
p_row_l(2:end,1,:) = p_row(1:end-1,1,:);

p_row_l(1,2,:) = Up_f(x,y+My-1,:);
p_row_l(2:end,2,:) = p_row(1:end-1,2,:);

p_row_r(end,1,:) = Up_f(x+Mx-1,y,:);
p_row_r(1:end-1,1,:) = p_row(2:end,1,:);

p_row_r(end,2,:) = Up_f(x+Mx-1,y+My-1,:);
p_row_r(1:end-1,2,:) = p_row(2:end,2,:); 

p_col_l(1,:,:) = Up_f(x-1,y+1:y+My-2,:);

p_col_l(2,1,:) = (3/4)*U_f(m_x-1,2,:) + (1/4)*U_f(m_x-1,3,:);
p_col_l(2,end,:) = (3/4)*U_f(m_x-1,m_y-1,:) + (1/4)*U_f(m_x-1,m_y-2,:);
p_col_l(2,2:end-1,:) = (1/2)*U_f(m_x-1,4:r:m_y-3,:) + (1/4)*U_f(m_x-1,5:r:m_y-2,:) + (1/4)*U_f(m_x-1,3:r:m_y-4,:);

p_col_r(1,1,:) = (3/4)*U_f(2,2,:) + (1/4)*U_f(2,3,:);
p_col_r(1,end,:) = (3/4)*U_f(2,m_y-1,:) + (1/4)*U_f(2,m_y-2,:);
p_col_r(1,2:end-1,:) = (1/2)*U_f(2,4:r:m_y-3,:) + (1/4)*U_f(2,5:r:m_y-2,:) + (1/4)*U_f(2,3:r:m_y-4,:);

p_col_r(2,:,:) = Up_f(x+Mx,y+1:y+My-2,:);

p_row = -((p_row_r + p_row) - (p_row + p_row_l))/(4*h); % f_(i+1/2) = f_(i+1) + f_i, f_(i-1/2) = f_i + f_(i-1)
p_col = -((p_col_r + p_col) - (p_col + p_col_l))/(4*h); % f_(i+1/2) - f_(i-1/2)

U_n = rhsFineGrid(U,0,G.child);

%assert(max(max(abs(U_n(2:r:end-1,1,:)-p_row(:,1,:)))) + max(max(abs(U_n(2:r:end-1,m_y,:)-p_row(:,2,:)))) + max(max(abs(U_n(1,2:r:end-1,:)-p_col(1,:,:)))) + max(max(abs(U_n(m_y,2:r:end-1,:)-p_col(2,:,:)))) < 1.0e-8)

%% Test 2: mesh refinement on points outside grid g flux

h = G.child.h;
U = G.child.u;
r = ratio;

U_f = f(U);
U_g = g(U);
Up_g = g(G.u);

m_x = G.child.m_x;
m_y = G.child.m_y;
x = G.child.location(3); 
y = G.child.location(4);
Mx = (m_x-3)/r +3;
My = (m_y-3)/r +3;

p_row_d = zeros(Mx-2,2,4);
p_row_u = zeros(Mx-2,2,4);
p_col_d = zeros(2,My-2,4);
p_col_u = zeros(2,My-2,4);

p_row = [U_g(2:r:m_x-1,1,:),U_g(2:r:m_x-1,m_y,:)];
p_col = [U_g(1,2:r:m_y-1,:);U_g(m_y,2:r:m_x-1,:)];

p_col_d(1,1,:) = Up_g(x,y,:);
p_col_d(1,2:end,:) = p_col(1,1:end-1,:); 

p_col_d(2,1,:) = Up_g(x+Mx-1,y,:);
p_col_d(2,2:end,:) = p_col(2,1:end-1,:);

p_col_u(1,end,:) = Up_g(x,y+My-1,:);
p_col_u(1,1:end-1,:) = p_col(1,2:end,:);

p_col_u(2,end,:) = Up_g(x+Mx-1,y+My-1,:);
p_col_u(2,1:end-1,:) = p_col(2,2:end,:); 

p_row_d(:,1,:) = Up_g(x+1:x+Mx-2,y-1,:);

p_row_d(1,2,:) = (3/4)*U_g(2,m_y-1,:) + (1/4)*U_g(3,m_y-1,:);
p_row_d(end,2,:) = (3/4)*U_g(m_x-1,m_y-1,:) + (1/4)*U_g(m_x-2,m_y-1,:);
p_row_d(2:end-1,2,:) = (1/2)*U_g(4:r:m_x-3,m_y-1,:) + (1/4)*U_g(5:r:m_x-2,m_y-1,:) + (1/4)*U_g(3:r:m_x-4,m_y-1,:);

p_row_u(1,1,:) = (3/4)*U_g(2,2,:) + (1/4)*U_g(3,2,:);
p_row_u(end,1,:) = (3/4)*U_g(m_x-1,2,:) + (1/4)*U_g(m_x-2,2,:);
p_row_u(2:end-1,1,:) = (1/2)*U_g(4:r:m_x-3,2,:) + (1/4)*U_g(5:r:m_x-2,2,:) + (1/4)*U_g(3:r:m_x-4,2,:);

p_row_u(:,2,:) = Up_g(x+1:x+Mx-2,y+My,:);

p_row = -((p_row_u + p_row) - (p_row + p_row_d))/(4*h); % f_(i+1/2) = f_(i+1) + f_i, f_(i-1/2) = f_i + f_(i-1)
p_col = -((p_col_u + p_col) - (p_col + p_col_d))/(4*h); % f_(i+1/2) - f_(i-1/2)

U_n = rhsFineGrid(U,0,G.child);

assert(max(max(abs(U_n(2:r:end-1,1,:)-p_row(:,1,:)))) + max(max(abs(U_n(2:r:end-1,m_y,:)-p_row(:,2,:)))) + max(max(abs(U_n(1,2:r:end-1,:)-p_col(1,:,:)))) + max(max(abs(U_n(m_y,2:r:end-1,:)-p_col(2,:,:)))) < 1.0e-12)

%% Test 3: Points outside grid diffusion

h = G.child.h;
U = G.child.u;
r = ratio;

m_x = G.child.m_x;
m_y = G.child.m_y;
x = G.child.location(3); 
y = G.child.location(4);
Mx = (m_x-3)/r +3;
My = (m_y-3)/r +3;

p_row_l = zeros(Mx-2,2,4);
p_row_r = zeros(Mx-2,2,4);

p_col_l = zeros(2,My-2,4);
p_col_r = zeros(2,My-2,4);


p_row = [U(2:r:m_x-1,1,:),U(2:r:m_x-1,m_y,:)];
p_col = [U(1,2:r:m_y-1,:);U(m_y,2:r:m_x-1,:)];

p_row_l(1,1,:) = G.u(x,y,:);
p_row_l(2:end,1,:) = p_row(1:end-1,1,:);

p_row_l(1,2,:) = G.u(x,y+My-1,:);
p_row_l(2:end,2,:) = p_row(1:end-1,2,:);

p_row_r(end,1,:) = G.u(x+Mx-1,y,:);
p_row_r(1:end-1,1,:) = p_row(2:end,1,:);

p_row_r(end,2,:) = G.u(x+Mx-1,y+My-1,:);
p_row_r(1:end-1,2,:) = p_row(2:end,2,:); 

p_col_l(1,:,:) = G.u(x-1,y+1:y+My-2,:);

p_col_l(2,1,:) = (3/4)*U(m_x-1,2,:) + (1/4)*U(m_x-1,3,:);
p_col_l(2,end,:) = (3/4)*U(m_x-1,m_y-1,:) + (1/4)*U(m_x-1,m_y-2,:);
p_col_l(2,2:end-1,:) = (1/2)*U(m_x-1,4:r:m_y-3,:) + (1/4)*U(m_x-1,5:r:m_y-2,:) + (1/4)*U(m_x-1,3:r:m_y-4,:);

p_col_r(1,1,:) = (3/4)*U(2,2,:) + (1/4)*U(2,3,:);
p_col_r(1,end,:) = (3/4)*U(2,m_y-1,:) + (1/4)*U(2,m_y-2,:);
p_col_r(1,2:end-1,:) = (1/2)*U(2,4:r:m_y-3,:) + (1/4)*U(2,5:r:m_y-2,:) + (1/4)*U(2,3:r:m_y-4,:);

p_col_r(2,:,:) = G.u(x+Mx,y+1:y+My-2,:);

p_row = -((p_row_r - p_row) - (p_row - p_row_l));
p_col = -((p_col_r - p_col) - (p_col - p_col_l));

U_diff = diffusion(1,U,G.u,x,y,r,h);

assert(max(max(abs(U_diff(2:r:end-1,1,:)-p_row(:,1,:)))) + max(max(abs(U_diff(2:r:end-1,m_y,:)-p_row(:,2,:)))) + max(max(abs(U_diff(1,2:r:end-1,:)-p_col(1,:,:)))) + max(max(abs(U_diff(m_y,2:r:end-1,:)-p_col(2,:,:)))) < 1.0e-12)

%% Test 4: Diffusion on boundary and inner points

h = G.child.h;
U = G.child.u;
r = ratio;

m_x = G.child.m_x;
m_y = G.child.m_y;
x = G.child.location(3); 
y = G.child.location(4);

U_pluss = zeros(m_x,m_y,4);

U_pluss(m_x-1,2:r:m_y-1,:) = U(m_x,2:r:m_y-1,:);
U_pluss(m_x-1,3:r:m_y-2,:) = (1/2)*U(m_x,2:r:m_y-3,:) + (1/2)*U(m_x,4:r:m_y-1,:);

U_pluss(2:m_x-2,2:m_y-1,:) = U(3:m_x-1,2:m_y-1,:);


U_minus = zeros(m_x,m_y,4);

U_minus(2,2:r:m_y-1,:) = U(1,2:r:m_y-1,:);
U_minus(2,3:r:m_y-2,:) = (1/2)*U(1,2:r:m_y-3,:) + (1/2)*U(1,4:r:m_y-1,:);

U_minus(3:m_x-1,2:m_y-1,:) = U(2:m_x-2,2:m_y-1,:);

U_n = -(U_pluss(2:m_x-1,2:m_y-1,:) - U(2:m_x-1,2:m_y-1,:)) - (-(U(2:m_x-1,2:m_y-1,:) - U_minus(2:m_x-1,2:m_y-1,:)));

U_diff = diffusion(1,U,G.u,x,y,r,h);

assert(max(max(max(abs(U_n-U_diff(2:m_x-1,2:m_y-1,:))))) < 1.0e-12)

%% Test 5: diffusion on g flux

h = G.child.h;
U = G.child.u;
r = ratio;

m_x = G.child.m_x;
m_y = G.child.m_y;
x = G.child.location(3); 
y = G.child.location(4);
Mx = (m_x-3)/r +3;
My = (m_y-3)/r +3;

U_pluss = zeros(m_x,m_y,4);

U_pluss(2:r:m_x-1,m_y-1,:) = U(2:r:m_x-1,m_y,:);
U_pluss(3:r:m_x-2,m_y-1,:) = (1/2)*U(2:r:m_x-3,m_y,:) + (1/2)*U(4:r:m_x-1,m_y,:);

U_pluss(2:m_x-1,2:m_y-2,:) = U(2:m_x-1,3:m_y-1,:);


U_minus = zeros(m_x,m_y,4);

U_minus(2:r:m_x-1,2,:) = U(2:r:m_x-1,1,:);
U_minus(3:r:m_x-2,2,:) = (1/2)*U(2:r:m_x-3,1,:) + (1/2)*U(4:r:m_x-1,1,:);

U_minus(2:m_x-1,3:m_y-1,:) = U(2:m_x-1,2:m_y-2,:);

U_n = -(U_pluss(2:m_x-1,2:m_y-1,:) - U(2:m_x-1,2:m_y-1,:)) - (-(U(2:m_x-1,2:m_y-1,:) - U_minus(2:m_x-1,2:m_y-1,:)));

U_diff = diffusion(1,permute(U,[2,1,3]),permute(G.u,[2,1,3]),y,x,r,h);
U_diff = permute(U_diff,[2,1,3]);

assert(max(max(max(abs(U_n-U_diff(2:m_x-1,2:m_y-1,:))))) < 1.0e-12)