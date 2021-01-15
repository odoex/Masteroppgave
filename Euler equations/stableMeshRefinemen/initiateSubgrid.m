function [G] = initiateSubgrid(G,r)
%INITIATESUBGRID Initiates the subgrids of coarser grids 
%   The fine grid must include a layer of gridpoints from the coarse grid
%   outside the boundary of the refinement in order to get a stable
%   transition. The grid is first initiated in the inner part by the
%   function exactSolutionEuler(). Then four vectors are initiated with the
%   solution of the equation around the fine subgrid. These are stored
%   outside the boundary in the fine grid vector u, but the cornerpoints
%   are not included as they are not needed. 

    h = G.h*r; % stepsize of parentgrid
    m_x = G.m_x;
    m_y = G.m_y;
    locx = [G.location(3),round(G.location(3)+(m_x-1))];
    locy = [G.location(4),round(G.location(3)+(m_y-1))];

    x = linspace(G.location(1),G.location(1) + (G.m_x-1)*G.h,G.m_x)';
    y = linspace(G.location(2),G.location(2) + (G.m_y-1)*G.h,G.m_y)';
    G.u = exactSolEuler(x,y,0);

    % Modifying the grid to include the frame of points around the boundary

    u(2:G.m_x+1,2:G.m_y+1,:) = G.u;
    
    u(2:m_x+1,1,:)=exactSolEuler(x,G.location(2),0);
    u(2:m_x+1,m_y+2,:)=exactSolEuler(x,G.location(2)+(G.m_y-1)*G.h,0);
    u(1,2:m_y+1,:)=exactSolEuler(G.location(1),y,0);
    u(m_y+2,2:m_y+1,:)=exactSolEuler(G.location(1) + (G.m_x-1)*G.h,y,0);
    
    G.location = [(locx(1)-2)*h,(locy(1)-2)*h,locx(1)-1,locy(1)-1];
    G.m_x = G.m_x+2;
    G.m_y = G.m_y+2;
    G.u = u;

end

