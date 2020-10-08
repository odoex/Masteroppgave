% Opdatere ved forrige eller neste tidssteg?
% Spørsmål: Har virkelig randen (på det grove gridet) så mye å si at
% metoden ikke får riktig konvergens dersom 

% 1. sjekk at konvergensen er riktig nå, før du går videre til å
% finne ut av rand og mesh refinement. 

% 2. Sjekk konvergens på det fine gridet med bcs fra det grove gridet. 
% - tror først

% 3. Sjekk om boundary conditions er riktige (sjekk i det hele tatt ut litt
% om boundary conditions). Denne settes på pause ettersom riktig konvergens
% er oppnådd med de eksisterende bcs. (Finn ut hva det vil si at de ikke
% opdaterer seg. Lær også litt mer om hvordan stabilitenen påvirkes rundt
% dette.

% 4. Endre metoden til å ikke ta med videre boundary fra fint grid. (Gir
% mening å ikke ta med denne, men skjønner fortsatt ikke helt hvorfor den
% ble helt shit når jeg gjorde det).

% Intervals in time
t_0 = 0;
t_n = 1;
t = t_0;

% Number of points in space and time for main grid
m = 5;
m_x = m;
m_y = m;
n = 50*t_n;

% time steps for coarse grid
k = (t_n-t_0)/(n-1);

% PDE coefficients and exact solution
a = 0.5;
b = 1;

f = @(x,y,s) sin(x - a*s) + sin(y - b*s);


% GRID CREATION

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);
G.t=0;
h=G.h;

% Fine grid: 
%   Area to be covered by fine grid: Ex locx = [a,b] where a<b
locx = [2,4];
locy = [2,4];
ratio = 16;

location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, n);
G_1.t = 0;
G.child = G_1;

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


error = calculateError(G.child,a,b);
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
%     zlim([-2 2])
%     %hold on
%     pause
% end

% Plot exact solution and approximated solution
% figure 
% mesh(sol);


% [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
% %[X1,Y1] = meshgrid(G.child.location(1):G.child.h:G.child.location(1)+G.child.h*(G.child.m-1));
% 
 %mesh(X,Y,G.u)
 [X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
 mesh(X,Y,G.u)
 hold on
% %mesh(X1,Y1,G.child.u)

disp(G.child.h)
% figure 
% mesh(E);