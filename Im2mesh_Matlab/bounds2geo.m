function bounds2geo( bounds, path_to_geo, opt )
% bounds2geo: generate geo file for Gmsh based on polygonal boundary
%
% geo file is a text file that contains Gmsh’s own scripting language. 
% It s used to define geometry and mesh generation instructions for Gmsh.
%
% usage1:
%   path_to_geo = 'test.geo';
%   bounds2geo( bounds, path_to_geo );
%
% usage2:
%   path_to_geo = 'test.geo';
%   opt = [];   % reset
%   opt.sizeMax = 5;
%   bounds2geo( bounds, path_to_geo, opt );
%
% usage3:
%   path_to_geo = 'test.geo';
%   opt = [];   % reset
%   opt.sizeMin = 0.1;
%   opt.sizeMax = 500;
%   opt.algthm = 6;
%   opt.recombAll = 0;
%   opt.recombAlgthm = 3;
%   opt.eleOrder = 1;
%   opt.scalingFactor = 1;
%   opt.num_split = 0;
%   bounds2geo( bounds, path_to_geo, opt );
%
% usage4:
%   path_to_geo = 'test.geo';
%   opt = [];   % reset
%   opt.sizeMin = 0.1;
%   opt.sizeMax = 500;
%   opt.algthm = 6;
%   opt.recombAll = 0;
%   opt.recombAlgthm = 3;
%   opt.eleOrder = 1;
%   opt.scalingFactor = 1;
%   opt.num_split = 0;
%   opt.grad_mode = 1;
%   opt.sizeAtBound = 2;
%   opt.sizeSlope = 0.25;
%   bounds2geo( bounds, path_to_geo, opt );
%
% input:
%  bounds - a nested cell array of 2d polygonal boundaries.
%            Polygons in bounds{i} belong to the i-th part or phase.
%            bounds{i}{j} is one of the polygons in the i-th part.
%            bounds{i}{j} is a n-by-2 array for x y coordinates of vertices
%            in a polygon. You can use 
%            plot( bounds{i}{j}(:,1), bounds{i}{j}(:,2) ) to view the
%            polygon. Use plotBounds( bounds ) to view all polygons.
%
%  path_to_geo - file name of geo file
%
%  opt - a structure array. It is the options for bounds2geo.
%         It stores parameter settings for bounds2geo.
% 
%  opt.sizeMin - Minimum mesh element size. Default value: 0.1
%				 Gmsh command: Mesh.MeshSizeMin
% 
%  opt.sizeMax - Maximum mesh element size. Default value: 500
%				 Gmsh command: Mesh.MeshSizeMax
% 
%  opt.algthm - 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 
%       3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 
%       8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 
%       11: Quasi-structured Quad).
%       Please refer to section 1.2.1 in Gmsh manual about mesh algorithm. 
%       Default value: 6
%		Gmsh command: Mesh.Algorithm
% 
%  opt.recombAll - Boolean. Value: 0 or 1. Whether to apply recombination 
%                  algorithm to all surfaces, ignoring per-surface spec. 
%                  Default value: 0
%				   Gmsh command: Mesh.RecombineAll
% 
%  opt.recombAlgthm - Mesh recombination algorithm (0: simple, 1: blossom, 
%                     2: simple full-quad, 3: blossom full-quad). 
%                     Default value: 3
%					  Gmsh command: Mesh.RecombinationAlgorithm
% 
%  opt.eleOrder - 1 or 2. Default value: 1
%				  Gmsh command: Mesh.ElementOrder
%
%  opt.scalingFactor - Global scaling factor applied to the saved mesh. 
%                      Default value: 1.
%					   Gmsh command: Mesh.ScalingFactor
%
%  opt.num_split - Number of splitting for refining mesh.
%                  Each triangle is split into four new sub-triangles.
%                  Default value: 0
%				   Gmsh command: RefineMesh;
% 
%  opt.grad_mode - Gradient mode. Integer. Value: 0, 1, 2.
%				   A parameter introduced by me.
%                  Default value: 0. 
%                  When grad_mode is 0, mesh size field is not specified. 
%                  When grad_mode is 1, the element size linearly increases
%                  as the distance from constraint edges (i.e., polygonal 
%                  boundary) increases. See the equation below.
%                  Element size = SizeAtBound + SizeSlope × distance from edges.
%                  When grad_mode is 2, the mesh size field is created 
%                  based on the local-feature-size function of MESH2D.
% 
%  opt.sizeAtBound - Element size at constraint edges (i.e., polygonal 
%                    boundary). 
%                    Default value: 500
% 
%  opt.sizeSlope - Slope in mesh size field.
%                  Default value: 0.2
%
%  Argument opt.sizeAtBound & opt.sizeSlope only work when opt.grad_mode
%  is set as 1 or 2.
%
%  Argument opt.local_max, opt.pnt_size, opt.interior_poly, and 
%  opt.hinitial only work when opt.grad_mode is set as 2. Their meaning and
%  usage is the same as the ones in function bounds2mesh (see demo17). They
%  are used to specify local feature size. Their default value are empty.
%      opt.local_max = [];
%      opt.pnt_size = [];
%      opt.interior_poly = {};
%      opt.hinitial = [];
% 
%  opt.tf_includeMesh - Boolean. Value: 0 or 1. Whether to include meshing
%                       parameters in the generated geo file.
%                       Default value: 1
%
% Notes:
%  1. To generate triangular mesh, you can set opt.algthm to 1-7 and set 
%     opt.recombAll to 0.
%  2. To generate unstructured quadrilateral mesh, you can set opt.algthm 
%     to 8 and set opt.recombAll to 1. You can try different 
%     opt.recombAlgthm.
%  3. To generate quasi-structured quadrilateral mesh, you can set 
%     opt.algthm to 11 and set opt.recombAll to 0.
%  4. When opt.grad_mode is 0, you can adjust the value of opt.sizeMax to 
%     achieve better mesh quality. Typical value of opt.sizeMax is 4 to 20.
%  5. When opt.grad_mode is 1, you can adjust opt.sizeAtBound, 
%     opt.sizeSlope, and opt.sizeMax to achieve better mesh quality.
%  6. When opt.grad_mode is 2, you can set opt.SizeAtBound to a larger 
%     value (e.g., 10-500) and adjust opt.SizeSlope.
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % ---------------------------------------------------------------------
    % check the number of inputs
    if nargin == 2
        opt = [];
    elseif nargin == 3
        % use as input
    else
        error("Check the number of inputs");
    end

    %---------------------------------------------------------------------
    % verify field names and set values for opt
    opt = setOption( opt );

    %---------------------------------------------------------------------
    % check opt.grad_mode

    if opt.grad_mode == 2       % create mesh size field based on MESH2D
        hmax = opt.sizeMax;
        grad_limit = opt.sizeSlope;
        optLfs.bound_size = opt.sizeAtBound;
        optLfs.local_max = opt.local_max;
        optLfs.pnt_size = opt.pnt_size;
        optLfs.interior_poly = opt.interior_poly;
        optLfs.hinitial = opt.hinitial;
        path_to_pos = 'sizes.pos';
        
        [vlfs,hlfs] = getLfsMESH2D( bounds, hmax, grad_limit, optLfs );
        % save the mesh size field to .pos file
        printSizePos( vlfs, hlfs, path_to_pos );
    end

    %---------------------------------------------------------------------
    % convert a cell array of polygonal boundaries to a nesting cell array 
    % for storing multiple loops (Gmsh).
    [ C, point, line ] = bound2SurfaceLoop( bounds );
    
    %---------------------------------------------------------------------
    % Open a file for writing
    fid = fopen(path_to_geo, 'wW');
    
    fprintf(fid, '// Gmsh input file generated by Im2mesh package. \n');
    fprintf(fid, '\n');

    %---------------------------------------------------------------------
    % print point, line, line loop, plane surface, physical surface
    printGeometry( fid, C, point, line );

    %---------------------------------------------------------------------
    % print meshing parameters
    fprintf(fid, '// Mesh \n');
    fprintf(fid, 'Mesh.ScalingFactor = %.8f;\n',    opt.scalingFactor );
    fprintf(fid, '\n');

    if opt.tf_includeMesh == true
        fprintf(fid, 'Mesh.MeshSizeMin = %.3f;\n',            opt.sizeMin );
        fprintf(fid, 'Mesh.MeshSizeMax = %.3f;\n',            opt.sizeMax );
        fprintf(fid, 'Mesh.Algorithm = %d;\n',                opt.algthm );
        fprintf(fid, 'Mesh.RecombineAll = %d;\n',             opt.recombAll );
        fprintf(fid, 'Mesh.RecombinationAlgorithm = %d;\n',   opt.recombAlgthm );
        fprintf(fid, 'Mesh.ElementOrder = %d;\n',             opt.eleOrder );
        fprintf(fid, '\n');

        if opt.grad_mode == 1   
            % element size linearly increases as the distance from edges.
            fprintf(fid, 'Field[1] = Distance;\n');
            fprintf(fid, 'Field[1].EdgesList = {1:%d};\n', size(line,1) );
            fprintf(fid, 'Field[2] = MathEval;\n');
            fprintf(fid, 'Field[2].F = "F1 * %.4f + %.4f";\n', opt.sizeSlope, opt.sizeAtBound );
            fprintf(fid, 'Background Field = 2;\n');
            fprintf(fid, '\n');
            fprintf(fid, 'Mesh.MeshSizeFromPoints = 0;\n');
            fprintf(fid, 'Mesh.MeshSizeFromCurvature = 0;\n');
            fprintf(fid, 'Mesh.MeshSizeExtendFromBoundary = 0;\n');

        elseif opt.grad_mode == 2   
            % use mesh size from "sizes.pos" which is based on MESH2D.
            fprintf(fid, 'Mesh.MeshSizeFromPoints = 0;\n');
            fprintf(fid, 'Mesh.MeshSizeFromCurvature = 0;\n');
            fprintf(fid, 'Mesh.MeshSizeExtendFromBoundary = 0;\n');
            fprintf(fid, 'Merge "sizes.pos";\n');
            fprintf(fid, 'Field[1] = PostView;\n');
            fprintf(fid, 'Field[1].ViewTag = 1;\n');
            fprintf(fid, 'Background Field = 1;\n');
        end
        
        fprintf(fid, '\n');
        fprintf(fid, 'Mesh 2;\n');
        fprintf(fid, '\n');
        
        if opt.num_split > 0
            for i = 1: opt.num_split
                fprintf(fid, 'RefineMesh;\n');
            end
            fprintf(fid, '\n');
        end
    end

    %---------------------------------------------------------------------
    fprintf(fid, '// set format\n');
    fprintf(fid, 'Mesh.Format=10;\n');
    fprintf(fid, '\n');
    
    %---------------------------------------------------------------------
    fprintf(fid, '// End of geo file\n');
    
    fclose(fid);

    disp('bounds2geo Done! Check the geo file!');
    %---------------------------------------------------------------------
end

function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    %---------------------------------------------------------------------
    % initialize new_opt with default field names & value 
    new_opt.sizeMin = 0.1;
    new_opt.sizeMax = 500;
    new_opt.algthm = 6;
    new_opt.recombAll = 0;
    new_opt.recombAlgthm = 3;
    new_opt.eleOrder = 1;
    new_opt.scalingFactor = 1;
    new_opt.num_split = 0;

    new_opt.grad_mode = 0;
    new_opt.sizeAtBound = 500;
    new_opt.sizeSlope = 0.2;
    
    new_opt.local_max = [];
    new_opt.pnt_size = [];
    new_opt.interior_poly = {};
    new_opt.hinitial = [];

    new_opt.tf_includeMesh = true;

    %---------------------------------------------------------------------
    if isempty(opt)
        return
    end

    if ~isstruct(opt)
        error("opt is not a structure array. Not valid input.")
    end

    %---------------------------------------------------------------------
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
    %---------------------------------------------------------------------
end

function [vlfs,hlfs] = getLfsMESH2D( bounds, hmax, grad_limit, opt )
% getLfsMESH2D
%
%     opt.bound_size
%     opt.local_max
%     opt.pnt_size
%     opt.interior_poly
%     opt.hinitial

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
    
    optLfs.kind = 'delaunay';   % Method for mesh-size functions
    
    optLfs.dhdx = grad_limit;   % dhdx is scalar gradient-limit
                                % default +0.2500
    
    optLfs.disp = +inf;

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

    % ---------------------------------------------------------------------
end

function printSizePos(vlfs, hlfs, filename)
% printSizePos: Export a "size cloud" to a Gmsh .pos file
%
%   printSizePos(vlfs, hlfs)                writes to 'size.pos'
%   printSizePos(vlfs, hlfs, filename)      writes to the file you supply
%
%   vlfs : N-by-2 matrix of (x,y) coordinates
%   hlfs : N-by-1 vector of target element sizes at those points
%
%   The function produces a file whose contents look like
%       View "h_target" {
%         SP(x, y, z){h};
%         ...
%       };
%

    % ---------------------------------------------------------------------
    % argument checking
    narginchk(2,3);
    
    if size(vlfs,2) ~= 2
        error('vlfs must be N-by-2 (x,y).');
    end

    if numel(hlfs) ~= size(vlfs,1)
        error('hlfs must have the same number of rows as vlfs.');
    end
    
    if nargin < 3 || isempty(filename)
        filename = 'sizes.pos';
    end
    
    % ---------------------------------------------------------------------
    vlfs(:, 3) = zeros(size(vlfs,1),1);

    % ---------------------------------------------------------------------
    % file writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open %s for writing.', filename);
    end

    fprintf(fid, 'View "h_target" {\n');
    fmt = '  SP(%.15g, %.15g, %.15g){%.15g};\n';
    fprintf(fid, fmt, [vlfs hlfs].');   % write all points at once
    fprintf(fid, '};\n');
    fclose(fid);
    % ---------------------------------------------------------------------
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

function printGeometry(fid, C, point, line )
% printGeometry

    %---------------------------------------------------------------------
    % ---------------------------
    % 1. Print out the Points
    % ---------------------------
    for i = 1:size(point, 1)
        x = point(i, 1);
        y = point(i, 2);
        % Gmsh syntax: Point(ID) = {x, y, z, lc};
        fprintf(fid, 'Point(%d) = {%.6f, %.6f, 0};\n', ...
                i, x, y );
    end
    fprintf(fid, '\n');
    
    % ---------------------------
    % 2. Print out the Lines
    % ---------------------------
    for i = 1:size(line, 1)
        p1 = line(i, 1);
        p2 = line(i, 2);
        % Gmsh syntax: Line(ID) = {startPt, endPt};
        fprintf(fid, 'Line(%d) = {%d, %d};\n', i, p1, p2);
    end
    fprintf(fid, '\n');
    
    fprintf(fid, '//----------------------------------------\n');
    fprintf(fid, '\n');

    %---------------------------------------------------------------------
    % Initialize counters for Gmsh IDs
    phyID     = 1;  % Physical surface ID
    surfaceID = 1; % Plane surface ID
    loopID    = 1;  % Line loop ID
    
    % Loop over each physical surface
    for i = 1:numel(C)
        % We will collect the plane surfaces (their IDs) that go into this physical surface
        surfacesInThisPhysical = [];
        
        % Each C{i} is a cell array of plane surfaces
        for j = 1:numel(C{i})
            
            % We will collect the loop IDs in this plane surface
            loopsInThisSurface = [];
            
            % Each C{i}{j} is a cell array of loops
            for k = 1:numel(C{i}{j})
                lineIndices = C{i}{j}{k};  % Nx1 array of line indices for the k-th loop
                
                % -------------------------------
                % 1. Print the Line Loop
                % -------------------------------
                % Gmsh syntax: Line Loop(loopID) = {line1, line2, ...};
                fprintf(fid, 'Line Loop(%d) = {', loopID);
                % Print each line index, separated by commas
                for n = 1:numel(lineIndices)
                    if n == numel(lineIndices)
                        fprintf(fid, '%d};\n', lineIndices(n));
                    else
                        fprintf(fid, '%d, ', lineIndices(n));
                    end
                end
                
                % Store this loop ID so we can reference it in the plane surface
                loopsInThisSurface = [loopsInThisSurface, loopID];
                loopID = loopID + 1;
                
            end % k
            
            % -------------------------------
            % 2. Print the Plane Surface
            % -------------------------------
            % A plane surface can have one or more loops: the first is the outer boundary,
            % and subsequent loops can be holes (with negative sign in advanced usage).
            % But here we assume each plane surface uses all the loops in loopsInThisSurface.
            
            % Gmsh syntax: Plane Surface(surfaceID) = {loop1, loop2, ...};
            fprintf(fid, 'Plane Surface(%d) = {', surfaceID);
            for u = 1:numel(loopsInThisSurface)
                if u == numel(loopsInThisSurface)
                    fprintf(fid, '%d};\n', loopsInThisSurface(u));
                else
                    fprintf(fid, '%d, ', loopsInThisSurface(u));
                end
            end
            fprintf(fid, '\n');
            
            % Store this surface ID so we can reference it in the physical surface
            surfacesInThisPhysical = [surfacesInThisPhysical, surfaceID];
            surfaceID = surfaceID + 1;
            
        end % j
        
        % -------------------------------
        % 3. Print the Physical Surface
        % -------------------------------
        % Gmsh syntax: Physical Surface(phyID) = {surface1, surface2, ...};
        
        fprintf(fid, 'Physical Surface(%d) = {', phyID);
        for v = 1:numel(surfacesInThisPhysical)
            if v == numel(surfacesInThisPhysical)
                fprintf(fid, '%d};\n', surfacesInThisPhysical(v));
            else
                fprintf(fid, '%d, ', surfacesInThisPhysical(v));
            end
        end
        fprintf(fid, '\n');
        fprintf(fid, '//----------------------------------------\n');
        fprintf(fid, '\n');

        phyID = phyID + 1;
        
    end % i
    %---------------------------------------------------------------------
end

