% Number of points in space and time for main grid
m = 5;

m_x = m;
m_y = m;
n = 10;

% Intervals in time
t_0 = 0;
t_n = 1;

% time steps for coarse grid
k = (t_n-t_0)/(n-1);

% PDE coefficients and exact solution
a = 0.5;
b = 1;

f = @(x,y,s) sin(x - a*s) + sin(y - b*s);


% Grid creation 

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);

% Fine grid: 
%   Area to be covered by fine grid: Ex locx = [a,b] where a<b
locx = [2,4];
locy = [2,4];
ratio = 3;

location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k/ratio, (locx(2)-locx(1))*ratio +1, n);
G.child = G_1;


% Creating solution vectors
G.u = createSolutionVector(G);
G.child.u = createSolutionVector(G.child);

t = t_0;
   
[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1)); 
mesh(X,Y,G.u) 
hold on

% Running the scheme: 
U = finiteVolume(G,a,b,t_0,t_n,f);

% Boundary conditions i det fine gridet holder tilbake fordi de er ved
% forrige tidssteg i forrige coarse grid mens det fine gridet beveger seg
% litt etter litt frem i tid. Hvordan fikser jeg dette?

error = calculateError(G,a,b);
%error1 = calculateError(G.child);

% Exact solution
%sol = sin(X - a*t_n) + sin(Y - b*t_n);


%disp(U);
%disp(sol);

%E = abs(U-sol');
% error = norm(E);
disp(error);


% Plot exact solution and approximated solution
% figure 
% mesh(sol);

% [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
% %[X1,Y1] = meshgrid(G.child.location(1):G.child.h:G.child.location(1)+G.child.h*(G.child.m-1));
% 
% mesh(X,Y,G.u) 
% hold on
% %mesh(X1,Y1,G.child.u)


% figure 
% mesh(E);