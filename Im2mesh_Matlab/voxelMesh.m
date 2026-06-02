function [vert,ele,tnum,vert2,ele2] = voxelMesh( im, opt )
% voxelMesh: Convert 3d multi-phase image to voxel-based finite element 
% mesh (8-node brick element). See demo19 for usage example.
%
% usage:
%   [vert,ele,tnum,vert2,ele2] = voxelMesh( im );
%   % OR
%   opt.select_phase = [1 3];
%   [vert,ele,tnum,vert2,ele2] = voxelMesh( im, opt );
%
% input: 
%   im - is grayscale uint8 3d image.
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
%     ele: Mesh elements (for linear element). For brick elements, 
%           it s a Ne-by-8 matrix, where Ne is the number of elements in 
%           the mesh. Each row in ele contains the indices of the nodes 
%           for that mesh element.
%     
%     tnum: Label of phase. Ne-by-1 array, where Ne is the number of 
%           elements
%       tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%     
%     vert2: Mesh nodes (for quadratic element). It’s a Nn-by-2 matrix.
%     
%     ele2: Mesh elements (for quadratic element). For brick 
%           elements, it s a Ne-by-20 matrix.
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    %----------------------------------------------------------------------
    % check the number of inputs
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    im = flip( im, 1 ); % FEM software use right-hand coordinate
    im = rot90(im, -1);
    
    num_row = size( im, 1 );
    num_col = size( im, 2 );
    num_sli = size( im, 3 );    % slice
    integer_type = getIntType( num_col, num_row, num_sli );
    
    %----------------------------------------------------------------------
    % get unique intensities from image
    %----------------------------------------------------------------------
    intensity = unique( im(:) );     % column vector
    
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
    % get 8-node numbering of each element
    % This section is improved by Gemini based on my own code.
    % My own code is easier to understand but much slower. :D
    %----------------------------------------------------------------------
    % Find indices of all voxels belonging to the target phases
    % 'ismember' returns a mask of the same size as 'im'
    mask = ismember( im, intensity );
    valid_lin_ind = find( mask );
    num_ele = length( valid_lin_ind );
    
    % Calculate Phase IDs (tnum)
    % Use the 2nd output of ismember to map image values to phase indices
    [~, loc] = ismember( im(valid_lin_ind), intensity );
    tnum = cast(loc, 'uint8'); 
    
    % Calculate Nodal Connectivity (ele)
    if num_ele > 0
        % Convert linear indices to subscripts [row, col, sli]
        [r, c, s] = ind2sub([num_row, num_col, num_sli], valid_lin_ind);
        
        % Define strides (increments) for the node grid (num_row+1 x num_col+1)
        stride_row = 1;
        stride_col = num_row + 1;
        stride_sli = (num_row + 1) * (num_col + 1);
        
        % Calculate the linear index of the first corner (Node 1) for EVERY
        % element. Formula: (c-1)*stride_col + r + (s-1)*stride_sli
        n1 = (c-1)*stride_col + r + (s-1)*stride_sli;
        
        % Define offsets for the 8 corners relative to Node 1
        node_offsets = [0, ...                                           % Node 1
                        stride_col, ...                                  % Node 2
                        stride_col + stride_row, ...                     % Node 3
                        stride_row, ...                                  % Node 4
                        stride_sli, ...                                  % Node 5
                        stride_sli + stride_col, ...                     % Node 6
                        stride_sli + stride_col + stride_row, ...        % Node 7
                        stride_sli + stride_row];                        % Node 8
        
        % Use implicit expansion (broadcasting) to create the element 
        % matrix. n1 is (NumEle x 1), offsets is (1 x 8) -> Result is 
        % (NumEle x 8)
        ele = cast(n1 + node_offsets, integer_type);
    else
        error('Num of element < 1')
    end
    
    %----------------------------------------------------------------------
    % get all node numbering
    %----------------------------------------------------------------------
    unique_node_ind_v = unique(ele);

    % get list of node coordinates, corresponding to unique_node_ind_v
    % nodecoor_list(i,:) = [ node_numbering, x, y ]
    nodecoor_list = getNodelist( unique_node_ind_v, num_col, num_row, num_sli );
    
    %----------------------------------------------------------------------
    % update node numbering in ele by mapping: nodecoor_list(i,1) -> i
    % so we can safely discard the 1st column of nodecoor_list in next step
    %----------------------------------------------------------------------
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
        for j = 1: 8
            ele(i,j) = mapping( ele(i,j) );
        end
    end
    
    %----------------------------------------------------------------------
    % x y z coordinates of vertices
    %----------------------------------------------------------------------
    vert = nodecoor_list(:,2:4);

    %----------------------------------------------------------------------
    % convert linear to quadratic element
    %----------------------------------------------------------------------
    [vert2, ele2] = insertNode3d(vert, ele);
    
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

function integer_type = getIntType( num_col, num_row, num_sli )
% get the suitable integer type for storing node number

    total_num_node = (num_row+1)*(num_col+1)*(num_sli+1);
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

function nodecoor_list = getNodelist( unique_node_ind_v, num_col, num_row, num_sli )
% getNodelist: Vectorized calculation of node coordinates
% This function is improved by Gemini based on my own code.
% My own code is easier to understand but slower. :D
%
% Inputs:
%   unique_node_ind_v: Vector of unique linear node indices
%   num_col, num_row, num_sli: Dimensions of the voxel grid (elements)
%
    
    % Ensure input is a column vector
    unique_node_ind_v = unique_node_ind_v(:);
    
    % Define the grid dimensions for nodes.
    % Note: A grid of RxCxS elements has (R+1)x(C+1)x(S+1) nodes.
    dim_nodes = [num_row + 1, num_col + 1, num_sli + 1];
    
    % Convert linear indices to Grid Subscripts (r, c, s)
    [r, c, s] = ind2sub(dim_nodes, unique_node_ind_v);
    
    % Map subscripts to physical coordinates. Formula: Coord = Index - 0.5
    % This means Index 1 -> 0.5, Index 2 -> 1.5, etc.
    % Meshgrid maps Row indices (dim 1) to Y and Col indices (dim 2) to X.
    y = r - 0.5; 
    x = c - 0.5; 
    z = s - 0.5; 
    
    % Combine into final list: [node_ID, x, y, z]
    % Concatenate the vectors directly.
    nodecoor_list = [double(unique_node_ind_v), x, y, z];

end