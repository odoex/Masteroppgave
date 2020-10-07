% Number of points in space and time for main grid
m = 17;
m_x = m;
m_y = m;
n = 50;

% Opdatere ved forrige eller neste tidssteg?

% Intervals in time
t_0 = 0;
t_n = 1;
t = t_0;

% time steps for coarse grid
k = (t_n-t_0)/(n-1);

% PDE coefficients and exact solution
a = 0.5;
b = 1;

f = @(x,y,s) sin(x - a*s) + sin(y - b*s);


% Grid creation 

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);
G.t=0;
h=G.h;

% Fine grid: 
%   Area to be covered by fine grid: Ex locx = [a,b] where a<b
locx = [2,7];
locy = [2,7];
ratio = 2;

location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, n);
G_1.t = 0;
%G.child = G_1;

% Creating solution vectors
G.u = createSolutionVector(G);
if G.child ~= 0
    G.child.u = createSolutionVector(G.child);
end

figure
   
% [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1)); 
% mesh(X,Y,G.u) 
% hold on

% Running the scheme: 
G = finiteVolume(G,a,b,t_0,t_n,f);


error = calculateError(G,a,b);
%error1 = calculateError(G.child);

% Exact solution
%sol = sin(X - a*t_n) + sin(Y - b*t_n);


%disp(U);
%disp(sol);

%E = abs(U-sol');
% error = norm(E);
disp(error);
disp(G.h);

figure

[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
interval = linspace(t_0,t_n,n);
% for i = 1:length(interval)
%     mesh((sin(X-a*interval(i)) + sin(Y-b*interval(i)))');
%     hold on
%     pause
% end
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