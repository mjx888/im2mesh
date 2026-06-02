function [vert,ele,tnum,vert2,ele2] = pixelMesh( im, opt )
% pixelMesh: Convert 2d multi-phase image to pixel-based finite element 
% mesh (4-node quadrilateral element). See demo10 for usage example.
%
% usage:
%   [vert,ele,tnum,vert2,ele2] = pixelMesh( im );
%   % OR
%   opt.select_phase = [1 3];
%   [vert,ele,tnum,vert2,ele2] = pixelMesh( im, opt );
%
% input: 
%   im - is grayscale uint8 2d image.
%
%   opt - a structure array. It is the options for pixelMesh.
%         It stores parameter settings for pixelMesh.
%
%   opt.select_phase - Select phase for meshing
%                      Parameter type: vector
%                      If 'select_phase' is [], all the phases will be
%                      chosen to perform meshing
%                      'select_phase' is an index vector for sorted 
%                      grayscales (ascending order) in an image.
%                      For example, an image with grayscales of 40, 90,
%                      200, 240, 255. If u're interested in 40, 200, and
%                      240, then set 'select_phase' as [1 3 4]. Those 
%                      phases corresponding to grayscales of 40, 200, 
%                      and 240 will be chosen to perform meshing.
%                      Default value: []
%                      See demo08 for usage example.
%
% output:
%   vert, ele define linear elements. vert2, ele2 define 2nd order elements.
%
%     vert: Mesh nodes (for linear element). It’s a Nn-by-2 matrix, where 
%           Nn is the number of nodes in the mesh. Each row of vert 
%           contains the x, y coordinates for that mesh node.
%     
%     ele: Mesh elements (for linear element). For quadrilateral elements, 
%           it s a Ne-by-4 matrix, where Ne is the number of elements in 
%           the mesh. Each row in ele contains the indices of the nodes 
%           for that mesh element.
%     
%     tnum: Label of phase. Ne-by-1 array, where Ne is the number of 
%           elements
%       tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%     
%     vert2: Mesh nodes (for quadratic element). It’s a Nn-by-2 matrix.
%     
%     ele2: Mesh elements (for quadratic element). For quadrilateral 
%           elements, it s a Ne-by-8 matrix.
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    %----------------------------------------------------------------------
    % check the number of inputs
    if nargin == 1
        opt = [];
    elseif nargin == 2
        % normal case
    else
        error("check the number of inputs");
    end

    % verify field names and set values for opt
    opt = setOption( opt );

    %----------------------------------------------------------------------
    % pre-process
    im = flip(im,1);	% in FEM software using right-hand coordinate, 
                        % to coincide with that, must flip in row direction
                        % so the origin of coordinates is at bottom-left

    num_row = size( im, 1 );
    num_col = size( im, 2 );
    integer_type = getIntType( num_col, num_row );
    
    %----------------------------------------------------------------------
    % get unique intensities from image
    intensity = unique( im );     % column vector

    % select phase
    if isempty(opt.select_phase)
        % = do nothing = all phases will be chosen
    elseif ~isvector(opt.select_phase)
        error("select_phase is not a vector")
    elseif length(opt.select_phase) > length(intensity)
        error("length of select_phase is larger than the number of phases")
    else
        % update intensity vector
        intensity = intensity( opt.select_phase );
    end

    %----------------------------------------------------------------------
    % total number of elements
    num_phase = length( intensity );
    num_ele = 0;
    for i = 1: num_phase
        num_ele = num_ele + sum(sum( im==intensity(i) ));
    end
    
    % initialize ele, tnum
    ele = zeros( num_ele, 4, integer_type );
    tnum = zeros( num_ele, 1, 'uint8' );

    %----------------------------------------------------------------------
    % get 4-node numbering of each element
    k = 1; % counter
    for i = 1: num_row
        for j = 1: num_col
            row = i;
            col = j;
            
            % check intensity
            if ~ismember( im(row,col), intensity )
                continue
            end
            
            % node numbering of 4 corner 
            Lind_4corner = [ 
                             (col-1)*(num_row+1) + row, ...
                             col*(num_row+1) + row, ...
                             col*(num_row+1) + row + 1, ...
                             (col-1)*(num_row+1) + row + 1
                             ];
            
            ele(k,:) = Lind_4corner;
            tnum(k,:) = find( im(row,col) == intensity );
            
            k = k + 1;
        end
    end

    %----------------------------------------------------------------------
    % get all node numbering 
    unique_node_ind_v = unique(ele);

    % get list of node coordinates, corresponding to unique_node_ind_v
    % nodecoor_list(i,:) = [ node_numbering, x, y ]
    nodecoor_list = getNodelist( unique_node_ind_v, num_col, num_row );
    
    %----------------------------------------------------------------------
    % update node numbering in ele by mapping: nodecoor_list(i,1) -> i
    % so we can safely discard the 1st column of nodecoor_list in next step
    
    % new_ele = ele;
    % 
    % for i = 1: size(nodecoor_list,1)
    %     old_ind = nodecoor_list(i,1);
    %     new_ind = i;
    %     new_ele( ele == old_ind ) = new_ind;
    % end
    % 
    % ele = new_ele;

    %----------------------------------------------------------------------
    % speed up the above process
    ind_vec = nodecoor_list(:,1);
    
    nInd = max(ind_vec);
    mapping = zeros(nInd,1);
    for i = 1: length(ind_vec)
        mapping( ind_vec(i) ) = i;
    end
    
    numE = size(ele,1);
    for i = 1: numE
        for j = 1: 4
            ele(i,j) = mapping( ele(i,j) );
        end
    end

    %----------------------------------------------------------------------
    % x y coordinates of vertices
    vert = nodecoor_list(:,2:3);

    %----------------------------------------------------------------------
    % convert linear to quadratic element
    [vert2, ele2] = insertNode(vert, ele);
    
    %----------------------------------------------------------------------
end


function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    % initialize new_opt with default field names & value 
    new_opt.select_phase = [];

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


function integer_type = getIntType( num_col, num_row )
% get the suitable integer type for storing node number

    total_num_node = (num_row+1)*(num_col+1);
    if total_num_node >0 && total_num_node < 2^64
        
        if total_num_node < 2^8
           integer_type = 'uint8';
        elseif total_num_node < 2^16
           integer_type = 'uint16';
        elseif total_num_node < 2^32
           integer_type = 'uint32';
        else
           integer_type = 'uint64';
        end
    else
        error('unexpected number of nodes');
    end
end


function nodecoor_list = getNodelist( unique_node_ind_v, num_col, num_row )
% getNodelist: get list of all nodes
%
% nodecoor_list(i,:) = [ node_number, x, y ]
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Nov 2019
%

    % generate x y coordinate of all nodes
    % can be accessed by X( row, col, sli ), Y( row, col, sli ), 
    xs = 0.5: num_col+0.5;
    ys = 0.5: num_row+0.5;
    [ X, Y ] = meshgrid( xs, ys );
    
    % reshape into vector
    % can be accessed by X(i), Y(i)
    X = X(:);
    Y = Y(:);
    
    % extract certain nodes
    X = X( unique_node_ind_v );
    Y = Y( unique_node_ind_v );
    
    num_node = length( unique_node_ind_v );
    % temporary list
    temp_list = zeros( num_node, 2 );
    
    for i = 1: num_node
        temp_list( i, : ) = [ X(i), Y(i) ];
    end
    
    % create point list, storing x y coordinate of all nodes
    % nodecoor_list(i,:) = [ node_number, x, y ]
    nodecoor_list = zeros( num_node, 3 );
    nodecoor_list( :, 1 ) = unique_node_ind_v;
    nodecoor_list( :, 2:3 ) = temp_list;

end







