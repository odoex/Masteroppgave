t_0 = 0;
t_n = 1;
x_0 = 0;
x_m = 1;

m = 50;
n = 2*m*t_n; 

t = t_0;

% time steps for coarse grid
k = (t_n-t_0)/(n-1);


% GRID CREATION

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);
G.t=0;
h=G.h;

% Solution vector
G.u = initialConditionsEuler(G);

%% Test 1: check Runge Kutta k1

U = G.u;

p = (U(:,:,4) - (1/2)*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1)).*(1.4-1);

F(:,:,1) = U(:,:,2); % rho*v1
F(:,:,2) = (U(:,:,2).^2)./U(:,:,1) + p;
F(:,:,3) = U(:,:,2).*U(:,:,3)./U(:,:,1);
F(:,:,4) = (U(:,:,4) + p).*U(:,:,2)./U(:,:,1);

Gg(:,:,1) = U(:,:,3); % rho*v2
Gg(:,:,2) = U(:,:,2).*U(:,:,3)./U(:,:,1);
Gg(:,:,3) = (U(:,:,3).^2)./U(:,:,1) + p;
Gg(:,:,4) = (U(:,:,4) + p).*U(:,:,3)./U(:,:,1);

k1 = - (F(3:50,2:49,:)-F(1:48,2:49,:))./(2*h) - (Gg(2:49,3:50,:)-Gg(2:49,1:48,:))./(2*h); 
k1_m = rhs(G.u,t,G);
F_m = f(G.u);

N = (F_m(2:49,2:49,:) - F(2:49,2:49,:));

%assert(max(max(max(F_m(2:49,2:49,:) - F(2:49,2:49,:)))) == 0)

assert(max(max(max(k1_m(2:49,2:49,:) - k1))) == 0) % Must be checked in another way, substract

