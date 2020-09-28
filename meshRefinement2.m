% Number of points in space and time
% For main grid m = 9
m = 9;
m_x = m;
m_y = m;
n = 1000;

% Intervals in time
t_0 = 0;
t_n = 1;

% time steps
k = (t_n-t_0)/(n-1);

a = 0.5;
b = 1;

% Exact solution
f = @(x,y,s) sin(x - a*s) + sin(y - b*s);

% Grid creation (Ikke lagre i skriptet ovenfor)
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);

% Fine grid: 
location_1 = [3*G.h,3*G.h,4,4];
ratio = 2;
m_1 = 5;

G_1 = Node(G, location_1, G.h/ratio, k, m_1, n); 
G.child = G_1;
% PROBLEM TIL I MORGEN: vil ikke lage node ut av G.child, setter den til null av en rar grunn
% prøv å lage noden først, se hva den blir og så sette G.child til å være
% den skapte node.
% Resultat: Det går an å lage ny node og sette G.child til å være denne.
% Dette er veldig rart. Aner ikke hvorfor matlab er sånn


% Creating solution vectors
G.u = createSolutionVector(G);
G.child.u = createSolutionVector(G.child);

% x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
% y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';
% 
% x1 = linspace(G_1.location(1),G_1.location(1) + (G_1.m-1)*G_1.h,G_1.m)';
% y1 = linspace(G_1.location(2),G_1.location(2) + (G_1.m-1)*G_1.h,G_1.m)';
% 
% [X,Y] = meshgrid(x,y);
% [X1,Y1] = meshgrid(x1,y1);

% Initial and boundary conditions
%F = sin(X) + sin(Y);

%g_x = @(s,z) (sin(z-a*s) + sin(-b*s)); 
%g_y = @(s,z) (sin(-a*s) + sin(z-b*s));

% Initializing the solution U (skal lagres i grid)
%U_n = F;

t = t_0;%linspace(t_0,t_n,n);

%U = finiteVolume(U,h,k,a,b,@(s) g_x(s,x),@(s) g_y(s,y),m_x,m_y,t_0,t_n);




    
    U = finiteVolume(G,a,b,t_0,t_n,f);
    
%     for i = 1:m_x
%         for j = 1:m_y
%             gx = g_x(t,x(i));
%             gy = g_y(t,y(j));
%             
%             U_n(i,j) = RK(U,i,j,h,k,a,b,g_x(t,x),g_y(t,y),m_x,m_y);
%         end 
%     end
%    U = U_n;



% Exact solution
%sol = sin(X - a*t_n) + sin(Y - b*t_n);


%disp(U);
%disp(sol);

%E = abs(U-sol');
% error = norm(E);
% disp(error);

% Plot exact solution and approximated solution
% figure 
% mesh(sol);

figure
mesh(G.u)

% figure 
% mesh(E);