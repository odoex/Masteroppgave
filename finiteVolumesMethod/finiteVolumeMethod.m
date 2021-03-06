% Number of points in space and time
m = 9; 
m_x = m;
m_y = m;
n = 54;

% Intervals in space and time
x_0 = 0;
x_m = 1;
y_0 = 0;
y_m = 1;
t_0 = 0;
t_n = 1;

% Space and time steps
h = (x_m-x_0)/(m-1);
k = (t_n-t_0)/(n-1);

% Grid creation (Ikke lagre i skriptet ovenfor)
%G = Node(0, [0.25,0.25], 0.125/2, k, m, n);
G = Node(0, [x_0,y_0], h, k, m, n);
x = linspace(G.location(1),G.location(1) + (G.m-1)*G.h,G.m)';
y = linspace(G.location(2),G.location(2) + (G.m-1)*G.h,G.m)';

%
h=G.h;
%

[X,Y] = meshgrid(x,y);

a = 0.5;
b = 1;

% Initial and boundary conditions
F = sin(X) + sin(Y);
ex_sol = @(s) sin(X-a*s) + sin(Y-b*s);

g_x = @(s,i) (sin(x(i)-a*s) + sin(-b*s)); 
g_y = @(s,j) (sin(-a*s) + sin(y(j)-b*s));

% Initializing the solution U (skal lagres i grid)
U = F;
%disp(U)
t = linspace(t_0,t_n,n);

disp(U)
for i = 1:(n-1)
    
    U = RK_ny(U,t(i),h,k,a,b,g_x,g_y,m_x,m_y);
%     disp(sin(X(4:10,4:10) - a*t(i+1)) + sin(Y(4:10,4:10) - b*t(1+i)))
%     disp(U(4:10,4:10))
%     disp(abs(U(4:10,4:10)- sin(X(4:10,4:10) - a*t(1+i)) - sin(Y(4:10,4:10) - b*t(i+1))))
%     disp(U)
%     if i == 1
%         disp(U)
%     end 
    
    %U(1,:) = sin(-a*t(i+1)) + sin(y-b*t(i+1));
    %U(:,1) = sin(x-a*t(i+1)) + sin(-b*t(i+1));
%     if i <= 1
%         disp(U)
%     end
%      disp(ex_sol(t(i+1))')
%    break
end

disp(U)
disp(t(i))

% Exact solution
sol = sin(X - a*t_n) + sin(Y - b*t_n);

% Error estimation
E = abs(U-sol');

%error = norm(E);
error = 0;
for i = 1:G.m
    for j = 1:G.m
        error = error + E(i,j)^2*h^2;
    end
end
error = sqrt(error);

disp(error)
%disp(E)

% error = 0;
% for i=1:m_x
%     for j=1:m_y
%         error = error + h*(U(i,j) - sol(i,j))^2;
%     end
% end
% 
% error = sqrt(error);
% 
% disp(error);


% Plot exact solution and approximated solution
figure 
mesh(sol);

%figure
mesh(U)
%disp(U)
figure 
mesh(E);