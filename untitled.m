Upluss(2:end-1,:,:) = ( U_f(i+1)+U_f(i) )/2
% U_i+(1/2)


U_pluss(1,2,:) = ( (3/4)*U_f(2,2,:)+(1/4)*U_f(2,3,:) ) + U_f(1,2,:);
U_pluss(1,m_y-1,:) = ( (3/4)*U_f(2,m_y-1,:)+(1/4)*U_f(2,m_y-2,:) ) + U_f(1,m_y-1,:);

% for i = [2,m_x-1;1,1]
%     U_pluss(i(1),i(2),:) = ( (3/4)*U_g(i(1),i(2)+1,:)+(1/4)*U_g(i(1)+1,i(2)+1,:) ) + U_g(i(1),i(2),:);
% end

% for j = [2,m_y-1]
%     U_pluss(m_x,j,:) = ( Up_f(x+Mx-1,j,:) + U_f(m_x,j,:) ); %%%
% end
% 
% for i = [2,m_x-1;m_y,m_y]
%     U_pluss(i(1),i(2),:) = ( Up_g(i(1),i(2)+1,:) + U_g(i(1),i(2),:) );
% end

for j = [y,y+My-1]
    U_pluss(m_x-1,j,:) = ( Up_f(x+Mx-1,j,:)+U_f(m_x-1,j,:) );
end

U_pluss(1,4:r:m_y-3,:) = ( (1/2)*U_f(2,4:r:m_y-3,:)+(1/4)*U_f(5:r:m_y-2)+(1/4)*U_f(3:r:m_y-4) ) + U_f(1,4:r:m_y-3,:);

U_pluss(m_x,2:r:m_y-1,:) = ( Up_f(x+Mx,y+1:y+My-2,:) + U_f(m_x,2:r:m_y-1,:) );

U_pluss(m_x-1,2:r:m_y-1,:) = ( U_f(m_x,2:r:m_y-1,:) + U_f(m_x-1,2:r:m_y-1,:) );

%For both 
% for j = [m_x,m_x,m_x-1,m_x-1;2,m_y-1,1,m_y]
%     U_pluss(j(1),j(2),:) = ( Up_f()+U_f(j(1),j(2),:) );
% end

% 
% U_pluss(1,2,:)     = ( (3/4)*U_f(2,2,:)+(1/4)*U_f(2,3,:) ) - U_f(1,2,:);
% U_pluss(1,m_y-1,:) = ( (3/4)*U_f(2,m_y-1,:)+(1/4)*U_f(2,m_y-2,:) ) - U_f(1,m_y-1,:);

%U_pluss(index(i,1),index(i,2),:) = ( (3/4)*U_f(index(i,1)+1,index(i,2),:)+(1/4)*U_f(index(i,1)+1,index(i,2)+1,:) ) + U_f(index(i,1),index(i,2),:);


U_f(3:end-2,2:end-1,:) = f(U(2:m_x-1,1:my,:));
U_f(1,:,:) = Up_f(x-1,y-1:y+My,:);
U_f(end,:,:) = Up_f(x+Mx,y-1:y+My,:);
U_f(:,1,:) = Up_f(x-1:x+Mx,y-1,:);
U_f(:,end,:) = Up_f(x-1:x+Mx,y+My,:);