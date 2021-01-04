
m_x = 50; 
m_y = 50; 

x = [0,1];
y = [0,1];
t = [0,0.15];

c = 331;

h = (x(2)-x(1))/(m_x-1); % When grid is structured 
k = h/(c);

n = floor(t(2)-t(1))/k+1;

G = Node(0,[x(1),y(1),1,1], h, k, m_x, m_y, n);

x = linspace(x(1),x(2),m_x)';
y = linspace(y(1),y(2),m_y)';

G.u = exactSolEuler(x,y,0);

%% Test 1: boundary values

[g_x0, g_xm, g_y0, g_ym] = boundary(G,0);
bc(1:m_x,:) = g_x0;
bc(m_x+1:2*m_x,:) = g_xm;
bc(2*m_x+1:3*m_y,:) = g_y0;
bc(3*m_x+1:4*m_y,:) = g_ym;

ex(1:m_x,:) = [G.u(:,1,1),G.u(:,1,2),G.u(:,1,3),G.u(:,1,4)];
ex(m_x+1:2*m_x,:) = [G.u(:,m_y,1),G.u(:,m_y,2),G.u(:,m_y,3),G.u(:,m_y,4)];
ex(2*m_x+1:3*m_y,:) = [G.u(1,:,1)',G.u(1,:,2)',G.u(1,:,3)',G.u(1,:,4)'];
ex(3*m_x+1:4*m_y,:) = [G.u(m_x,:,1)',G.u(m_x,:,2)',G.u(m_x,:,3)',G.u(m_x,:,4)'];

assert(norm(bc-ex)<1.0e-16)

%% Test 2: subgrid boundary values

G.u = exactSolEuler(x,y,0);

ratio = 2;

locx = [0.30,0.60];
locy = [0.30,0.60];
i = round(locx(1)/h)+1;
j = round(locy(1)/h)+1;

G1 = Node(G,[h*(i-1),h*(i-1),i,j], h/ratio, k, (round(locx(2)/h)-(i-1))*ratio+1, (round(locy(2)/h)-(j-1))*ratio+1, n); 
G1.t = 0;
G.child = G1;

x1 = linspace(G1.location(1),G1.location(1)+G1.h*(G1.m_x-1),G1.m_x);
y1 = linspace(G1.location(2),G1.location(2)+G1.h*(G1.m_y-1),G1.m_y);

[g_x0, g_xm, g_y0, g_ym] = boundary(G1,0);

for i = G1.location(3):(G1.location(3)+(G1.m_x-1)/2)
    
    A((i-G1.location(3)+1)*2-1,:) = g_x0((i-G1.location(3)+1)*2 -1,:) - permute(G.u(i,G1.location(3),:),[1,3,2]);
    A(G1.m_x+(i-G1.location(3)+1)*2-1,:) = g_xm((i-G1.location(3)+1)*2 -1,:) - permute(G.u(i,(G1.location(3)+(G1.m_x-1)/2),:),[1,3,2]);
        
	if (i > G1.location(3) && i < (G1.location(3)+(G1.m_x-1)/2))
        A((i-G1.location(3)+1)*2,:) = g_x0((i-G1.location(3)+1)*2,:) - (g_x0((i-G1.location(3)+1)*2+1,:) + g_x0((i-G1.location(3)+1)*2-1,:))/2;
        A(G1.m_x+(i-G1.location(3)+1)*2,:) = g_xm((i-G1.location(3)+1)*2,:) - (g_xm((i-G1.location(3)+1)*2+1,:) + g_xm((i-G1.location(3)+1)*2-1,:))/2;
	end
end
    
for j = G1.location(4):(G1.location(4)+(G1.m_y-1)/2)
	A(2*G1.m_x+(j-G1.location(4)+1)*2-1,:) = g_y0((j-G1.location(4)+1)*2 -1,:) - permute(G.u(G1.location(4),j,:),[1,3,2]);
	A(2*G1.m_x+G1.m_y+(j-G1.location(4)+1)*2-1,:) = g_ym((j-G1.location(4)+1)*2 -1,:) - permute(G.u((G1.location(4)+(G1.m_y-1)/2),j,:),[1,3,2]);
        
    if (j > G1.location(4) && j < (G1.location(4)+(G1.m_x-1)/2))
        A(2*G1.m_x+(j-G1.location(4)+1)*2,:) = g_y0((j-G1.location(4)+1)*2,:) - (g_y0((j-G1.location(4)+1)*2+1,:) + g_y0((j-G1.location(4)+1)*2-1,:))/2;
        A(2*G1.m_x+G1.m_y+(j-G1.location(4)+1)*2,:) = g_ym((j-G1.location(4)+1)*2,:) - (g_ym((j-G1.location(4)+1)*2+1,:) + g_ym((j-G1.location(4)+1)*2-1,:))/2;
    end
end


assert(norm(norm(A))<1.0e-16)
