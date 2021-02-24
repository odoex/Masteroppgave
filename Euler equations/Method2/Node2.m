classdef Node2 < handle
    % Node
    %   Node in a tree representing a hierarchy of grids. Each node is a
    %   grid with the following properties: 
    %       parent: the coarser grid above the current grid in the tree
    %       child: a finer grid below the current grid in the tree
    %       sibling: a grid of the same level with the same parent
    %       location: the bottom left coordinates of the grid
    %       h: space step in the grid
    %       k: time step in the grid
    %       t: the point in time where the grid currently is 
    %       n: number of time steps 
    %       m_x: number of space steps in x-direction
    %       m_y: number of space steps in y-direction
    %       u: solution vector for the grid
    %       dF: vector for coarse side of interface on coarse grid
    
    properties
        parent
        child
        sibling
        location
        h
        k
        t
        n
        m_x
        m_y
        u
        dF
    end
    
    methods
        function obj = Node(parent,location, h, k, m_x, m_y, n)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.parent = parent;
            obj.location = location;
            obj.h = h;
            obj.k = k;
            obj.m_x = m_x;
            obj.m_y = m_y;
            obj.n = n;
        end
        
        function set.child(obj, child)
            obj.child = child;
        end 
        
        function set.sibling(obj, sibling)
            obj.sibling = sibling;
        end
        
        function parent = get.parent(obj)
            parent = obj.parent;
        end
        %%%
%         function child = get.child(obj)
%             child = obj.parent;
%         end
        
    end
end

