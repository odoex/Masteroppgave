% Opdatere ved forrige eller neste tidssteg?
% Spørsmål: Har virkelig randen (på det grove gridet) så mye å si at
% metoden ikke får riktig konvergens dersom 

% 1. Konvergens for grovt grid gitt i convergence2.mat er 2

% 2. Konvergens på fint grid er ikke kontstant. Først stor og minker
% deretter. 
% - først veldig nøyaktig sammenlignet med ingen forfining, deretter har
% forfiningen mindre å si ettersom det ikke kommer ny informasjon. Verdiene
% som ligger rundt er bare nøyaktig til en viss grad. 

% 3. Boundary conditions stemmer. Selv med eksakte bcs har det fine gridet
% dårlig konvergens på det indre av det store gridet. Dersom det strekkes
% ut i kantene har det konvergens 2.  (Finn ut hva det vil si at de ikke
% opdaterer seg. Sjekk også litt mer om hvordan stabilitenen påvirkes rundt
% dette.)


% Idea: Since Euler euations are four, send it into the function four times
% instead of modifying the function fit the euler equations.

% 1. Med fint grid vel inni det grove og kun halverte ruter sjekk om
% metoden blir mer presis. Rettelse: helver hele gridet i fint og grovt.
% Sykt stress med fire sider. Sjekk med to sider først.

% 2. Endre fluxen i rhs funksjonen til å være en funksjon f og g??? Dette
% er kanskje ikke nødvendig siden dette er et test program som alltid
% kommer til å brukes på advection.

% 3. Hvis metoden er bedre: endre til helt generelle tilfeller. 

% 4. Undersøk tilfellet av interface som to boundaries. 

% Problem: Må ta vare på laget utenfor det lille gridet og overføre det til
% det grove. 

% Intervals in time
t_0 = 0;
t_n = 1;
t = t_0;

% Number of points in space and time for main grid
m = 25;
m_x = m;
m_y = m;
%n = 2*m*3*t_n; 

% time steps for coarse grid
k = (1/(m-1))/(2); %(t_n-t_0)/(n-1); % CHANGE BACK %k = (1/(m-1))/(2);
n = (t_n-t_0)/k;

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
locx = [10,20]; % Correct: this should be decided by location in interval
locy = [10,20];
ratio = 2;

location_1 = [(locx(1)-1)*G.h,(locy(1)-1)*G.h,locx(1),locy(1)];

G_1 = Node(G, location_1, G.h/ratio, k, (locx(2)-locx(1))*ratio +1, (locy(2)-locy(1))*ratio +1, n);
G_1.t = 0;
G.child = G_1;


% SOLUTION VECTORS
G.u = createSolutionVector(G);
if G.child ~= 0
    G.child.u = createSolutionVector(G.child);
end


% % "Ghost-values" for fine grid
% G.child.g = zeros(G.child.m,4);
% G.child.g(1,1) = G.u(G.child.location(3),G.child.location(4)-1);
% G.child.g(1,2) = G.u(G.child.location(3)-1,G.child.location(4));
% 
% for i = 1:(G.child.m-1)/2
%     
%     G.child.g(2*i+1,1) = G.u(G.child.location(3)+i,G.child.location(4)-1);
%     G.child.g(2*i+1,2) = G.u(G.child.location(3)-1,G.child.location(4)+i);
% 
%     G.child.g(2*i,1) = (1/2)*G.child.g(2*i-1,1) + (1/2)*G.child.g(2*i+1,1);
%     G.child.g(2*i,2) = (1/2)*G.child.g(2*i-1,2) + (1/2)*G.child.g(2*i+1,2);
%     
% end

% Plotting the function

figure
[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m_x-1));  % Change: different size in x and y direction
mesh(X,Y,G.u)
hold on
if G.child ~= 0
    [X,Y] = meshgrid(G.child.location(1):G.child.h:G.child.location(1)+G.child.h*(G.child.m_x-1)); % Change: different size in x and y direction
    mesh(X,Y,G.child.u)
    hold on
end


% Running the scheme: 
G = finiteVolumeAdv(G,a,b,t_0,t_n,f);

%disp(G.u)
disp(G.t)

%Calculating the error
error = calculateError(G,a,b);

disp(error);

figure
[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m_x-1)); % Change: different size in x and y direction
interval = linspace(t_0,t_n,n);
mesh(X,Y,G.u)
