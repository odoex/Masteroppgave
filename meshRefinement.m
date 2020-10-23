% Intervals in time
t_0 = 0;
t_n = 1;
t = t_0;

% Number of points in space and time for main grid
m = 9;
m_x = m;
m_y = m;
n = 2*m*3*t_n; 

% time steps for coarse grid
k = (t_n-t_0)/(n-1);



%%%%%%%%%%%%


% PDE coefficients and exact solution
a = 0.5;
b = 1;

f = @(x,y,s) sin(x - a*s) + sin(y - b*s);

%%%%%%%%%%%



% GRID CREATION

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);
G.t=0;
h=G.h;

% Fine grid: 
%   Area to be covered by fine grid: Ex locx = [a,b] where a<b
locx = [3,7]; % Correct: this should be decided by location in interval
locy = [3,7];
ratio = 2;

location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, n);
G_1.t = 0;
G.child = G_1;


% SOLUTION VECTORS
G.u = initialConditionsEuler(G); % Change with equation
if G.child ~= 0
    G.child.u = initialConditionsEuler(G.child);
end


% Plotting the function .. this requires four plots 

% figure
% [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1)); 
% mesh(X,Y,G.u)
% hold on
% [X,Y] = meshgrid(G.child.location(1):G.child.h:G.child.location(1)+G.child.h*(G.child.m-1));
% mesh(X,Y,G.child.u)
% hold on


% Running the scheme: 
G = finiteVolume(G,t_0,t_n);

disp(G.u)
disp(G.t)

%Calculating the error
error = calculateError(G,a,b);

disp(error);

% figure
[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
interval = linspace(t_0,t_n,n);
