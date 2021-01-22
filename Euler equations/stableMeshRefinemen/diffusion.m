function [U_diff] = diffusion(delta,U,Up,x,y,r,h)
%DIFFUSION Summary of this function goes here
%   Detailed explanation goes here

    m_x = length(U(:,1,1));
    m_y = length(U(1,:,1));
    variables = length(U(1,1,:));
    
    Mx = (m_x-3)/r+3;
    My = (m_y-3)/r+3;

    U_diff = zeros(m_x,m_y,variables);

    U_diff(1,2,:) = ( (3/4)*U(2,2,:)+(1/4)*U(2,3,:) ) - U(1,2,:);
    U_diff(1,m_y-1,:) = ( (3/4)*U(2,m_y-1,:)+(1/4)*U(2,m_y-2,:) ) - U(1,m_y-1,:); % top and bottom of left column at x = 1

%     for j = [y,y+My-1]
%         U_diff(m_x-1,j,:) = ( Up(x+Mx-1,j,:)+U(m_x-1,j,:) );
%     end

    U_diff(m_x-1,1,:) = ( Up(x+Mx-1,y,:) - U(m_x-1,1,:) ); % End of bottom row 
    U_diff(m_x-1,m_y,:) = ( Up(x+Mx-1,y+My-1,:) - U(m_x-1,m_y,:) ); % End of top row

    U_diff(1,4:r:m_y-3,:) = ( (1/2)*U(2,4:r:m_y-3,:)+(1/4)*U(2,5:r:m_y-2,:)+(1/4)*U(2,3:r:m_y-4,:) ) - U(1,4:r:m_y-3,:); % Left column at x = 1

    U_diff(m_x,2:r:m_y-1,:) = ( Up(x+Mx,y+1:y+My-2,:) - U(m_x,2:r:m_y-1,:) ); % Rigth column at x = m_x
    
    U_diff(m_x-1,2:r:m_y-1,:) = ( U(m_x,2:r:m_y-1,:) - U(m_x-1,2:r:m_y-1,:) );
    U_diff(m_x-1,3:r:m_y-2,:) = ( (1/2)*U(m_x,4:r:m_y-1,:)+(1/2)*U(m_x,2:r:m_y-3,:) ) - U(m_x-1,3:r:m_y-2,:); % right boundary
    
    for i=[1,m_y]
        U_diff(2:r:m_x-3,i,:) = U(4:r:m_x-1,i,:) - U(2:r:m_x-3,i,:); % Top and bottom row
    end
    
    U_diff(2:m_x-2,2:m_y-1,:) = U(3:m_x-1,2:m_y-1,:) - U(2:m_x-2,2:m_y-1,:);
    
    U_diff = -delta.*U_diff;
    
    % negative direction: For the flux f_i-1/2
    
    U_minus = zeros(m_x,m_y,variables);
    
    U_minus(m_x,2,:) = U(m_x,2,:) - ( (3/4)*U(m_x-1,2,:)+(1/4)*U(m_x-1,3,:) );
    U_minus(m_x,m_y-1,:) = U(m_x,m_y-1,:) - ( (3/4)*U(m_x-1,m_y-1,:)+(1/4)*U(m_x-1,m_y-2,:) ); % top and bottom of right column at x = m_x

    U_minus(2,1,:) = ( U(2,1,:) - Up(x,y,:) ); % Left end of bottom row 
    U_minus(2,m_y,:) = ( U(2,m_y,:) - Up(x,y+My-1,:) ); % Left end of top row

    U_minus(1,2:r:m_y-1,:) = U(1,2:r:m_y-1,:) - Up(x-1,y+1:y+My-2,:); % Left column at x = 1

    U_minus(m_x,4:r:m_y-3,:) = U(m_x,4:r:m_y-3,:) - ( (1/2)*U(m_x-1,4:r:m_y-3,:)+(1/4)*U(m_x-1,5:r:m_y-2,:)+(1/4)*U(m_x-1,3:r:m_y-4,:) ); % Rigth column at x = m_x
    
    U_minus(2,2:r:m_y-1,:) = ( U(2,2:r:m_y-1,:) - U(1,2:r:m_y-1,:) );
    U_minus(2,3:r:m_y-2,:) = U(2,3:r:m_y-2,:) - ( (1/2)*U(1,4:r:m_y-1,:)+(1/2)*U(1,2:r:m_y-3,:) ); % left boundary
    
    for i=[1,m_y]
        U_minus(4:r:m_x-1,i,:) = U(4:r:m_x-1,i,:) - U(2:r:m_x-3,i,:); % Top and bottom row
    end
    
    U_minus(3:m_x-1,2:m_y-1,:) = U(3:m_x-1,2:m_y-1,:) - U(2:m_x-2,2:m_y-1,:);
    
    U_minus = -delta.*U_minus;
    
    U_diff = U_diff - U_minus;
    
    H = ones(m_x,m_y,variables);
    H = H.*1/(2*h);
    H(2:m_x-1,1,:) = 1/(4*h);
    H(2:m_x-1,m_y,:) = 1/(4*h);
    H(1,2:m_y-1,:) = 1/(4*h);
    H(m_x,2:m_y-1,:) = 1/(4*h);
    H(2,2:m_y-1,:) = 1/(3*h);
    H(m_x-1,2:m_y-1,:) = 1/(3*h);
    
    U_diff = U_diff.*H; % multiplying wiht 1/2 inside H
    
    % Ide: lage kun det som skal trekkes fra/ legges til, og så trekke det
    % fra/legge det til det gjeldende punktet (hele matrisen) etterpå.
    
end

