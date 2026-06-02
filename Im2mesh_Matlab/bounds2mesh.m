function [vert,tria,tnum,vert2,tria2,etri] = bounds2mesh( bounds, hmax, grad_limit, opt )
% bounds2mesh: generate meshes of parts defined by polygonal boundary.
%        	   Mesh generator: MESH2D (https://github.com/dengwirda/mesh2d)
%              See demo17 for the usage example of function bounds2mesh.
%
% usage:
%   [vert,tria,tnum] = bounds2mesh( bounds, hmax, grad_limit );
%   [vert,tria,tnum,vert2,tria2,etri] = bounds2mesh( bounds, hmax, grad_limit );
%   [vert,tria,tnum,vert2,tria2,etri] = bounds2mesh( bounds, hmax, grad_limit, opt );
% 
% input:
%   bounds - a nested cell array of 2d polygonal boundaries.
%            Polygons in bounds{i} belong to the i-th part or phase.
%            bounds{i}{j} is one of the polygons in the i-th part.
%            bounds{i}{j} is a n-by-2 array for x y coordinates of vertices
%            in a polygon. You can use 
%            plot( bounds{i}{j}(:,1), bounds{i}{j}(:,2) ) to view the
%            polygon. Use plotBounds( bounds ) to view all polygons.
%   
%   hmax - Maximum mesh edge lengths. This is an approximate upper bound
%		   on the mesh edge lengths. Range: > 0
%   
%   grad_limit - A limit on the gradient of mesh-size function.
%                Range: > 0. Typical value: 0.15 - 0.5
%                http://persson.berkeley.edu/pub/persson04gradlim.pdf
%
%   opt - a structure array. It is the options for bounds2mesh.
%         It stores extra parameter settings for bounds2mesh.
%
%   opt.mesh_kind - Method used to create mesh-size functions based on 
%                   an estimate of the "local-feature-size".
%                   Value: 'delaunay' or 'delfront' 
%                   Delaunay-refinement or Frontal-Delaunay
%                   Default value: 'delaunay'
%
%   opt.tf_smooth - Boolean. Value: 0 or 1. Whether improve triangulation 
%                   quality by adjusting the vertex positions and mesh 
%                   topology (hill-climbing type optimisation). 
%                   Default value: 1
%
%   opt.num_split - Number of splitting for refining mesh.
%                   Each triangle is split into four new sub-triangles.
%                   Default value: 0
%
%   opt.outerbound_size - Element size at the outermost polygonal boundary.
%                         This is used to refine mesh near the outermost 
%                         boundary.
%                         If you don t need to refine mesh near outermost 
%                         boundary, you can set opt.outerbound_size to 0.
%                         Default value: 0
%
%   opt.bound_size - Element size at constraint edges (i.e., polygonal 
%                    boundary). This is used to refine mesh near all 
%                    polygonal boundary.
%                    If you don't need to refine mesh near boundary, you 
%                    can set opt.bound_size to 0.
%                    Default value: 0
%                    In each edge, the seeds are inserted according to the 
%                    following equation. 'len' is the length of an edge.
%                       numSegment = round( len / targetSize );
%                       actualSize = len / numSegment;
%
%   opt.local_max - n-by-2 array, used to specify max mesh size in a part.
%                   '[2, 0.5; 3, 0.15]' means that max mesh size in part 2
%                   is 0.5; max mesh size in part 3 is 0.15.
%                   When set as [], this parameter will be ignored.
%                   Default value: []
%
%   opt.pnt_size - p-by-3 array, used to specify mesh size at a point.
%                  Each row is a point and its corresponding mesh size.
%                  '[2, 3, 0.4; 5, 1, 0.15]' means that mesh size at point
%                  (2, 3) is 0.4; mesh size at point (5, 1) is 0.15.
%                  Default value: []
%
%   opt.interior_poly - c-by-1 cell array, used to specify interior edge
%                       constraints. See Im2mesh package demo17.
%                       opt.interior_poly{i} is a n-by-2 array, for x and y
%                       coordinates of a 2d polyline.
%                       Experimental feature. May fail in some cases.
%                       Default value: {}
%
%   opt.disp - Verbosity. Set as 'inf' to mute verbosity.
%              Default value: 10
%
%   opt.hinitial - initial mesh size when creating local feature size
%                  function. It's used to avoid atrraction issue when 
%                  specifying mesh size at multiple points.
%                  Default value: []
%
%   opt.tf_preserve - Boolean. Value: 0 or 1. Whether preserve the
%                     input edges with minimal refinement. When this option
%                     is set as 1, MESH2D attempts to retain the initial 
%                     edges without further subdivision, and edges are 
%                     split only to satisfy basic geomertical conformance.
%                     Default value: 0
%
% output:
%   vert, tria define linear elements. vert2, tria2 define 2nd order elements.
%
%   vert: Mesh nodes (for linear element). It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   tria: Mesh elements (for linear element). For triangular elements, 
%         it s a Ne-by-3 matrix, where Ne is the number of elements in 
%         the mesh. Each row in tria contains the indices of the nodes 
%         for that mesh element.
%     
%   tnum: Label of phase. Ne-by-1 array, where Ne is the number of elements
%       tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%     
%   vert2: Mesh nodes (for quadratic element). It’s a Nn-by-2 matrix.
%     
%   tria2: Mesh elements (for quadratic element). For triangular 
%          elements, it s a Ne-by-6 matrix.
%
%   etri: C-by-2 array of constraining edges, where each row defines an edge
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % ---------------------------------------------------------------------
    % check the number of inputs
    if nargin == 3
        opt = [];
    elseif nargin == 4
        % use as input
    else
        error("Check the number of inputs");
    end
    
    % ---------------------------------------------------------------------
    % verify field names and set values for opt
    opt = setOption( opt );

    % ---------------------------------------------------------------------
    % Add uniform seeds to the outer boundary according to opt.outerbound_size

    if opt.outerbound_size > 0
        bounds = refineOuterBoundary( bounds, opt.outerbound_size );
    end

    % ---------------------------------------------------------------------
    % Add uniform seeds to all boundaries according to opt.bound_size
    
    if opt.bound_size > 0
        space = opt.bound_size;  % space between seeds
        n_digit = 2;
        
        for i = 1: length(bounds)
	        for j = 1: length(bounds{i})
		        bounds{i}{j} = insertEleSizeSeed( bounds{i}{j}, space );
                
		        bounds{i}{j} = round( bounds{i}{j}, n_digit );
	        end
        end
    end
    
    % ---------------------------------------------------------------------
    % create geometry (planar straight-line graph)

    % get nodes and edges (cell array) of polygonal boundary
    [ poly_node, poly_edge ] = getPolyNodeEdge( bounds );
    
    % create planar straight-line graph
    [ node, edge, part ] = regroup( poly_node, poly_edge );
    
    % node, edge - array. Nodes and edges of all polygonal boundary
    % node, edge doesn't record phase info. Phase info is recorded by part.
    % node - V-by-2 array. x,y coordinates of vertices. 
    %        Each row is one vertex.
    % edge - E-by-2 array. Node numbering of two connecting vertices of
    %        edges. Each row is one edge.
    % part - cell array. Used to record phase info.
    %        part{i} is edge indexes of the i-th phase, indicating which 
    %        edges make up the boundary of the i-th phase.
    
    % ---------------------------------------------------------------------
    % add extra nodes according to opt.pnt_size

    if ~isempty( opt.pnt_size )
        node = addNode( node, opt.pnt_size );
    end

    % ---------------------------------------------------------------------
    % add extra nodes and interior edges according to opt.interior_poly
    
    if ~isempty( opt.interior_poly )
        [ node, edge ] = addInteriorEdge( node, edge, opt.interior_poly );
    end

    % ---------------------------------------------------------------------
    % Create mesh size function
    % LFSHFN2 routine is used to create mesh-size functions based on an 
    % estimate of the local-feature-size associated with a polygonal domain 
    
    optLfs.kind = opt.mesh_kind;    % Method for mesh-size functions
                                    % Value: 'delaunay' or 'delfront' 
                                    % default value is 'delaunay'
    
    optLfs.dhdx = grad_limit;   % dhdx is scalar gradient-limit
                                % default +0.2500
    
    optLfs.disp = opt.disp;
    
    if opt.hinitial <= 0,  opt.hinitial = [];  end
    hinitial = opt.hinitial;
    
    [vlfs,tlfs,hlfs] = lfshfn2( node, edge, part, optLfs, hinitial );
    
    % ---------------------------------------------------------------------
    % modify mesh size field hlfs according to opt.local_max
    % opt.local_max is mesh size in a part

    if ~isempty(opt.local_max)
        hlfs = setLocalMax(node,edge,part,hmax,vlfs, hlfs, opt.local_max );
        
        hlfs = limhfn2(vlfs,tlfs,hlfs,grad_limit);  % push gradient limits
    end
    
    % ---------------------------------------------------------------------
    % modify mesh size field hlfs according to opt.pnt_size
    % opt.pnt_size is mesh size at a point
    
    if ~isempty( opt.pnt_size )
        hlfs = setPntSize( vlfs, hlfs, opt.pnt_size );
        
        hlfs = limhfn2(vlfs,tlfs,hlfs,grad_limit);  % push gradient limits
    end

    % ---------------------------------------------------------------------
    hlfs = min(hmax,hlfs);

    [slfs] = idxtri2(vlfs,tlfs);
    
    % ---------------------------------------------------------------------
    % Mesh genenration using refine2
    
    hfun = @trihfn2;
    
    % set opt for refine2
    optRef.disp = opt.disp;
    if opt.tf_preserve == 1
        optRef.ref1 = 'preserve';
    end

    [vert,etri,tria,tnum] = refine2(node,edge,part,optRef,hfun, ...
                                    vlfs,tlfs,slfs,hlfs);
    
    % ---------------------------------------------------------------------
    % Smooth mesh
    % SMOOTH2 routine provides iterative mesh "smoothing" capabilities, 
    % seeking to improve triangulation quality by adjusting the vertex 
    % positions and mesh topology.
    
    if opt.tf_smooth ~= 0
        optSmo.disp = opt.disp;
        [vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum,optSmo);
    end
    
    % ---------------------------------------------------------------------
    % Refine by splitting (quadtree-type mesh refinement)
    % TRIDIV2 routine is used to refine existing trangulations. 
    % Each triangle is split into four new sub-triangles, such that 
    % element shape is preserved.
    
    if opt.num_split > 0
        vnew = vert;    
        enew = etri;
        tnew = tria;
        
        for i = 1: opt.num_split
            [vnew,enew,tnew,tnum] = tridiv2(vnew,enew,tnew,tnum);
        end
        
        if opt.tf_smooth ~= 0
            optSmo.disp = opt.disp;
            [vnew,enew,tnew,tnum] = smooth2(vnew,enew,tnew,tnum,optSmo);
        end
        
        vert = vnew;
        etri = enew;
        tria = tnew;
    end
    
    % ---------------------------------------------------------------------
    % Get 2nd order elements
    [ vert2, tria2 ] = insertNode( vert, tria );

    % ---------------------------------------------------------------------
end


function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    % initialize new_opt with default field names & value 
    new_opt.mesh_kind = 'delaunay';
    new_opt.tf_smooth = true;
    new_opt.num_split = 0;
    new_opt.outerbound_size = 0;
    new_opt.bound_size = 0;
    new_opt.local_max = [];
    new_opt.pnt_size = [];
    new_opt.interior_poly = {};
    new_opt.disp = 10;
    new_opt.hinitial = [];
    new_opt.tf_preserve = false;

    if isempty(opt)
        return
    end

    if ~isstruct(opt)
        error("opt is not a structure array. Not valid input.")
    end

    % get the field names of opt
    nameC = fieldnames(opt);

    % verify field names in opt and set values in new_opt
    % compare the field name of opt with new_opt using for loop
    % if a field name of opt exist in new_opt, assign the that field value 
    % in opt to new_opt
    % if a field name of opt not exist in new_opt, show error

    for i = 1: length(nameC)
        if isfield( new_opt, nameC{i} )
            value = getfield( opt, nameC{i} );
            new_opt = setfield( new_opt, nameC{i}, value );
        else
            error("Field name %s in opt is not correct.", nameC{i});
        end
    end

end


function [ node, edge ] = addInteriorEdge( node, edge, polyCell )
% addInteriorEdge: add extra nodes and interior edges according to polyCell
% polyCell is a cell array of polyline

    % convert polylines (polyCell) to node, edge (PSLG)
    node_ex = [];
    edge_ex = [];
    
    for i = 1:length(polyCell)
        node_t = polyCell{i};      % t means temp
        nn = length(node_t);       % nn is the number of nodes
        edge_t = [(1:nn-1)', (2:nn)'];
        [ node_ex, edge_ex ] = joinNodeEdge( node_ex,edge_ex, node_t,edge_t );
    end
    
    % add to global
    [ node, edge ] = joinNodeEdge( node,edge, node_ex,edge_ex );
end


function hlfs = setLocalMax( node, edge, part, hmax, vlfs, hlfs, local_max )
% setLocalMax: modify mesh size field hlfs according to local_max
% local_max is mesh size in a part

    % create a vector for local max mesh size
    size_vec = hmax * ones( size(hlfs,1), 1 );
    
    numPart2Refine = size( local_max, 1 );
    
    % set size_vec
    for i = 1: numPart2Refine
        idx = local_max(i,1);    % part index
        lmax = local_max(i,2);   % local max mesh size
        
        tf_in = inpoly2( vlfs, node, edge(part{idx},:) );
        size_vec(tf_in) = lmax;
    end
    
    hlfs = min( size_vec, hlfs );
    
end

function hlfs = setPntSize( vlfs, hlfs, pnt_size )
% setPntSize: modify mesh size field hlfs according to pnt_size
% pnt_size is mesh size at a point
    
    xys = pnt_size( :, 1:2 );   % point x y
    lsize = pnt_size( :, 3 );   % local mesh size at a point

    [tf_vec, loc] = isvertex( xys, vlfs );
    if ~all(tf_vec)
        error('Wierd case. Probably the point is outside of polygon.');  
    end
    
    % modify mesh size field hlfs
    hlfs(loc) = lsize;
end

function node = addNode( node, pnt_size )
% addNode: add extra nodes according to pnt_size

    xys = pnt_size( :, 1:2 );
    tf_vec = isvertex( xys, node );
    
    xys = xys( ~tf_vec, : );    % find 'xys' not existing in 'node'
    node = [node; xys];         % append
end








