
% Rydd opp i koden
% Laget med [30,50] og m=65
% Er i hvert fall HELT sikker på de indre punktene på det fine gridet for
% de har jeg testet tre ganger nå (Selv om jeg klarte å gjøre den samme
% feilen to av gangene).

% Intervals in time
t_0 = 0;
t_n = 1;
t = t_0;

% Number of points in space and time for main grid
m = 200;
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
locx = [80,120]; % Correct: this should be decided by location in interval
locy = [80,120]; 
ratio = 2;

% locx(1) = locx(1)-1; % Adding anothen point to
% locy(1) = locy(1)-1;
% locx(2) = locx(2)+1; 
% locy(2) = locy(2)+1;

location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, (locy(2)-locy(1))*ratio +1, n);
G_1.t = 0;
G.child = G_1;


% SOLUTION VECTORS
G.u = createSolutionVector(G);
if G.child ~= 0
    
    G.child.u = createSolutionVector(G.child);
    u = zeros(G.child.m_x+2,G.child.m_y+2);
    u(2:G.child.m_x+1,2:G.child.m_y+1) = G.child.u;
    
    u(2:G.child.m_x+1,1)=exactSolAdv(linspace((locx(1)-1)*G.h,(locx(2)-1)*G.h,G.child.m_x)',(locy(1)-2)*G.h,0);
    u(2:G.child.m_x+1,G.child.m_y+2)=exactSolAdv(linspace((locx(1)-1)*G.h,(locx(2)-1)*G.h,G.child.m_x)',locy(2)*G.h,0);
    u(1,2:G.child.m_y+1)=exactSolAdv((locx(1)-2)*G.h,linspace((locy(1)-1)*G.h,(locy(2)-1)*G.h,G.child.m_y)',0);
    u(G.child.m_y+2,2:G.child.m_y+1)=exactSolAdv(locx(2)*G.h,linspace((locy(1)-1)*G.h,(locy(2)-1)*G.h,G.child.m_y)',0);
    
    G.child.location = [(locx(1)-2)*G.h,(locy(1)-2)*G.h,locx(1)-1,locy(1)-1];
    G.child.m_x = G.child.m_x+2;
    G.child.m_y = G.child.m_y+2;
    G.child.u = u;
    
end

% Running the scheme: 
G = finiteVolumeAdv(G,a,b,t_0,t_n,f);

% Plotting the solution
figure
[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m_x-1)); 
mesh(X,Y,G.u)
hold on

%Calculating the error
error = calculateError(G,a,b);

disp(error);
