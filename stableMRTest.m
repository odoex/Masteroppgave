% Creates two grids and inserts the fine grid into the coarse grid in the
% way of the stable mesh refinement. Can be used to test all functions with
% this purpose

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

%u = exactSolEuler(x,y,0); % Change equation here 
u = exactSolAdv(x,y,0);

% Points in time and time step

if length(u(1,1,:)) > 1
    maxSpeed = max(max(u(:,:,2)));
    c = 331; % speed of sound
else
    maxSpeed = 1;
    c = 0; % speed of sound
end

k = h/(ratio*(c + maxSpeed)); 

n = floor((t_n-t_0)/k)+1;

% GRID CREATION

% Coarse grid:
G = Node(0, [x(1),x(1),1,1], h, k, m_x, m_y, n);
G.t=t;
G.u = u;

locx = [round(locx(1)/h)+1,round(locx(2)/h)+1];
locy = [round(locy(1)/h)+1,round(locy(2)/h)+1];
location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, (locy(2)-locy(1))*ratio +1, n);
G_1.t = 0;

variables = length(u(1,1,:));
G_1.u = randi(10,G_1.m_x,G_1.m_y,variables);%initiateSubgrid(G_1,ratio);

G.child = G_1;

% Properties of the subgrid
h = G.child.h;
U = G.child.u;
r = ratio;

%U_n = rhsSubgrid(U,G.child); % Change function here
U_n = rhsFineGridAdv(U,0,G.child);

U_f = f_adv(U); % Matrix of f-flux % Change flux function here
U_g = g_adv(U); % Matrix of g-flux % Change flux function here
Up_f = f_adv(G.u);
Up_g = g_adv(G.u);

m_x = G.child.m_x;
m_y = G.child.m_y;
x = G.child.location(3);
y = G.child.location(4);
Mx = (m_x-3)/r +3;
My = (m_y-3)/r +3;

%% Test 1: Points outside grid f-flux

p_row_l = zeros(Mx-2,2,variables);
p_row_r = zeros(Mx-2,2,variables);
p_col_l = zeros(2,My-2,variables);
p_col_r = zeros(2,My-2,variables);

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

assert(max(max(abs(U_n(2:r:end-1,1,:)-p_row(:,1,:)))) + max(max(abs(U_n(2:r:end-1,m_y,:)-p_row(:,2,:)))) + max(max(abs(U_n(1,2:r:end-1,:)-p_col(1,:,:)))) + max(max(abs(U_n(m_y,2:r:end-1,:)-p_col(2,:,:)))) < 1.0e-8)

%% Test 2: Boundary and inner points

U_l = zeros(m_x-2,m_y-2,variables); % Smaller than U because they don't include points from outside the boundary
U_r = zeros(m_x-2,m_y-2,variables);

% finite volumes 
U_l(1,2:r:end-1,:) = (1/2)*U_f(1,2:r:m_y-3,:) + (1/2)*U_f(1,4:r:m_y-1,:); % Points inbetween the coarse gridpoints
U_l(1,1:r:end) = U_f(1,2:r:m_y-1,:); % Points overlapping the coarse gridpoints

U_l(2:end,:,:) = U_f(2:m_x-2,2:m_y-1,:); % The rest

U_r(end,2:r:end-1,:) = (1/2)*U_f(m_x,2:r:m_y-3,:) + (1/2)*U_f(m_x,4:r:m_y-1,:); % Points inbetween the coarse gridpoints
U_r(end,1:r:end) = U_f(m_x,2:r:m_y-1,:); % Points overlapping the coarse gridpoints

U_r(1:end-1,:,:) = U_f(3:m_x-1,2:m_y-1,:); % The rest

H = ones(length(U_r(:,1,1)),length(U_r(1,:,1)),length(U_r(1,1,:)));
H(2:end-1,:,:) = (1/(2*h)).*H(2:end-1,:,:);
H(1,:,:) = (1/(3*h)).*H(1,:,:);
H(end,:,:) = (1/(3*h)).*H(end,:,:);

U_t = -((U_r + U_f(2:end-1,2:end-1,:)) - (U_l + U_f(2:end-1,2:end-1,:))).*H;

assert(max(max(max(abs(U_t-U_n(2:end-1,2:end-1,:))))) < 1.0e-8)

%% Test 3: Points outside grid g-flux

p_row_d = zeros(Mx-2,2,variables);
p_row_u = zeros(Mx-2,2,variables);
p_col_d = zeros(2,My-2,variables);
p_col_u = zeros(2,My-2,variables);

p_row = [U_g(2:r:m_x-1,1,:),U_g(2:r:m_x-1,m_y,:)];
p_col = [U_g(1,2:r:m_y-1,:);U_g(m_y,2:r:m_x-1,:)];

p_col_d(1,1,:) = Up_g(x,y,:); %% 
p_col_d(1,2:end,:) = p_col(1,1:end-1,:); %%

p_col_d(2,1,:) = Up_g(x+Mx-1,y,:); %%
p_col_d(2,2:end,:) = p_col(2,1:end-1,:); %%

p_col_u(1,end,:) = Up_g(x,y+My-1,:); %%
p_col_u(1,1:end-1,:) = p_col(1,2:end,:); %%

p_col_u(2,end,:) = Up_g(y+My-1,x+Mx-1,:); %%
p_col_u(2,1:end-1,:) = p_col(2,2:end,:); %%

p_row_d(:,1,:) = Up_g(x+1:x+Mx-2,y-1,:); %%%%%

p_row_d(1,2,:) = (3/4)*U_g(2,m_y-1,:) + (1/4)*U_g(3,m_y-1,:); %%
p_row_d(end,2,:) = (3/4)*U_g(m_x-1,m_y-1,:) + (1/4)*U_g(m_x-2,m_y-1,:); %%
p_row_d(2:end-1,2,:) = (1/2)*U_g(4:r:m_x-3,m_y-1,:) + (1/4)*U_g(5:r:m_x-2,m_y-1,:) + (1/4)*U_g(3:r:m_x-4,m_y-1,:); %% 

p_row_u(1,1,:) = (3/4)*U_g(2,2,:) + (1/4)*U_g(3,2,:); %%
p_row_u(end,1,:) = (3/4)*U_g(m_x-1,2,:) + (1/4)*U_g(m_x-2,2,:); %%
p_row_u(2:end-1,1,:) = (1/2)*U_g(4:r:m_x-3,2,:) + (1/4)*U_g(5:r:m_x-2,2,:) + (1/4)*U_g(3:r:m_x-4,2,:); %%

p_row_u(:,2,:) = Up_g(x+1:x+Mx-2,y+My,:); %%

p_row = -((p_row_u + p_row) - (p_row + p_row_d))/(4*h); % f_(i+1/2) = f_(i+1) + f_i, f_(i-1/2) = f_i + f_(i-1)
p_col = -((p_col_u + p_col) - (p_col + p_col_d))/(4*h); % f_(i+1/2) - f_(i-1/2)

assert(max(max(abs(U_n(2:r:end-1,1,:)-p_row(:,1,:)))) + max(max(abs(U_n(2:r:end-1,m_y,:)-p_row(:,2,:)))) + max(max(abs(U_n(1,2:r:end-1,:)-p_col(1,:,:)))) + max(max(abs(U_n(m_y,2:r:end-1,:)-p_col(2,:,:)))) < 1.0e-8)