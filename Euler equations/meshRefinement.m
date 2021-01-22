format long

% Change three variables for advection: u (u for fine grid) and k (and
% remove maxSpeed) and maybe t_n and error estimation

% Interval in space 
x = [0,1];
y = [0,1];

m = 50; % > 50 Sjekk med 80 og 160 etterpå med mr
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
locx = [0.35,0.65];
locy = [0.35,0.65]; % Standard: [0.35,0.65];

% Solution vector
x = linspace(x(1),x(2),m_x)';
y = linspace(y(1),y(2),m_y)';

u = exactSolEuler(x,y,0);
%u = initialConditionsEulerConstant(x,y,0);
%u = exactSolAdv(x,y,0); 

% Points in time and time step
c = 331; % speed of sound

maxSpeed = max(max(u(:,:,2)));

k = h/(ratio*(c + maxSpeed)); 
%k = h/(ratio*2); % For advection

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
    
    x_1 = linspace(G_1.location(1),G_1.location(1) + (G_1.m_x-1)*G_1.h,G_1.m_x)';
    y_1 = linspace(G_1.location(2),G_1.location(2) + (G_1.m_y-1)*G_1.h,G_1.m_y)';
    G_1.u = exactSolEuler(x_1,y_1,0);
    %G_1.u = exactSolAdv(x_1,y_1,0);
    G.child = G_1;
end

% Plot of initial conditions
[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m_x-1));
X = X';
Y = Y';
figure
quiver(X,Y,G.u(:,:,2),G.u(:,:,3));
figure
mesh(X,Y,G.u(:,:,1))

% Running the scheme: 
G = finiteVolume(G,t_0,t_n);

% Plot of final result
figure
quiver(X,Y,G.u(:,:,2),G.u(:,:,3));
figure
mesh(X,Y,G.u(:,:,1))

% Error estimation
sol=exactSolEuler(linspace(0,1,G.m_x)',linspace(0,1,G.m_x)',G.t);
%sol=exactSolAdv(linspace(0,1,G.m_x)',linspace(0,1,G.m_x)',G.t);

E = abs(G.u-sol);
    
error = zeros(4,1);

l_n = length(u(1,1,:));
for l = 1:l_n
    for i = 1:G.m_x
        for j = 1:G.m_x
            error(l) = error(l) + E(i,j,l)^2*G.h^2;
        end
    end
    error2 = error;
    error(l) = sqrt(error(l));
end
    
disp(error)
if l_n > 1
    disp(sqrt(error2(1)+error2(2)+error2(3)+error2(4)))
end


% For m = 100 CFL <= 1.5