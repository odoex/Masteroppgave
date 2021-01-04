% GRID CREATION
m = 9;
n = 20;
k = 1/(n-1);

% Coarse grid: 
G = Node(0, [0,0,1,1], 1/(m-1), k, m, n);
G.t = 0;
h = G.h;


% Center of vortex
    x_c = 0.5;
    y_c = 0.5;

    % Constants
    Ma = 0.1; % Mach number 
    rho_inf = 1; % Freestream state
    T_inf = 273.15; % Freestream state
    R_g = 287.15; % Ideal gas constant
    gamma = 1.4; % Adiabatic exponent
    R_v = 0.1; % Radius of vortex 
    beta = 1; % Vortex strength
     v_inf = Ma*331;


u = initialConditionsEuler(G);

u_x = u(:,:,2)./u(:,:,1);
u2=u(1,1,2);
u1=u(1,1,1);
ur=u(1,1,2)/u(1,1,1);
u_y = u(:,:,3)./u(:,:,1);

f = @(x,y) ((x - x_c).^2 + (y - y_c).^2)./R_v;

dv_1 = @(x,y) -(v_inf*beta).*((y-y_c)./R_v).*exp(-f(x,y)./2);
dv_2 = @(x,y) (v_inf*beta).*((x-x_c)./R_v).*exp(-f(x,y)./2);

[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));
f1=-f(X,Y)./2;
e1=exp(-f(X,Y)./2);
y1=((Y-y_c)./R_v);
x1=((X-x_c)./R_v);
s1=-(v_inf*beta);
s2=(v_inf*beta);
a=s1*y1;
b=s2*x1;

figure;
contourf(X,Y,f(X,Y));

% k=1;
% v1=-k*((X - x_c).^2 + (Y - y_c).^2)./R_v;
% v2=k*((X - x_c).^2 + (Y - y_c).^2)./R_v;
% figure
% quiver(X,Y,v1,v2);




% figure; 
% quiver(X,Y,dv_1(X,Y),dv_2(X,Y));
% figure; 
% quiver(X,Y,dv_2(X,Y));

% v_1 = (2/R_v).*(Y-y_c);
% v_2 = -(2/R_v).*(X-x_c);

v_1 = v_inf + dv_1(X,Y);
v_2 = dv_2(X,Y);
figure;
% quiver(X,Y,v_1,v_2);
for i = 1:30
%     v_1 = (v_1 + dv_1(X,Y));
%     v_2 = v_2 + dv_2(X,Y);
    quiver(X,Y,v_1,v_2);
    pause
    % quiver(X,Y,f(X,Y));
    v_1 = (v_1 + dv_1(X,Y));
    v_2 = v_2 + dv_2(X,Y);
end


%[X,Y] = meshgrid(G.location(1):G.h:G.location(1)+G.h*(G.m-1));


% mesh(X,Y,u(:,:,1))
% hold on
% mesh(X,Y,u(:,:,2))
% hold on
% mesh(X,Y,u(:,:,3))
% hold on
% mesh(X,Y,u(:,:,1))
% hold on