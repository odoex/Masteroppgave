
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

% G.u = exactSolEuler(x,y,0);
G.u = zeros(G.m_x,G.m_y,4);

%% Test 1: projectFlux

locx = [0.35,0.65];
locy = [0.35,0.65];
r = 2;

ix = [round(locx(1)/(h))+1, round(locx(2)/h)+1];
jy = [round(locy(1)/(h))+1, round(locy(2)/(h))+1];

locx = [(ix(1)-1)*(h), (ix(2)-1)*(h)];
locy = [(jy(1)-1)*(h), (jy(2)-1)*(h)];

G.child = Node(G,[locx(1), locy(1), ix(1), jy(1)], h/2, k, (ix(2)-ix(1))*r+1, (jy(2)-jy(1))*r+1, n);

G.child.u = randi([1,10],[G.child.m_x, G.child.m_y,4]);

G = projectFlux(G,G.child);

ic = 1;
jc = 1; 
for i = G.child.location(3)+1:ix(2)-1
    for j = G.child.location(4)+1:jy(2)-1
        for l = 1:4
            A(1+i-G.child.location(3),1+j-G.child.location(4),l) = G.u(i,j,l) -G.child.u((i-G.child.location(3)+1)*2-1,(j-G.child.location(4)+1)*2-1,l);
        end
    end
end

assert(max(max(max(A)))<1.0e-16)

