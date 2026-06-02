function plotMeshes( vert, ele, tnum, color_code, opt )
% plotMeshes: plot triangular or quadrilateral finite element mesh
% Works for linear and quadratic elements. 
% Works for triangular and quadrilateral mesh. 
%
% Nodes must be counter-clockwise ordering in an linear element.
%
% See demo07 for usage examples.
%
% usage 1:
%   plotMeshes( vert, ele );        % one phase
%   plotMeshes( vert, ele, tnum );  % multiple phases
%
% usage 2:
%   plotMeshes( vert, ele, [], color_code );    % one phase
%   plotMeshes( vert, ele, tnum, color_code );  % multiple phases
%
% usage 3:
%   color_code = 2;
%   opt = [];   % reset
%   opt.mode = 2;
%   plotMeshes( vert,ele,tnum, color_code, opt )
%
% usage 4:
%   color_code = 2;
%   opt = [];   % reset
%   opt.mode = 1;
%   opt.wid = 0.5;
%   opt.alpha = 0.2;
%   opt.beta = 0.5;
%   opt.tf_gs = 1;
%   opt.q_thresh = 0.6;
%
%   plotMeshes( vert,ele,tnum, color_code, opt )
%
% input:
%   Argument tnum, color_code, and opt are optional.
%
%   vert: Mesh nodes. It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   ele: Mesh elements. 
%        For linear triangular elements, it s a Ne-by-3 matrix. 
%        For linear quadrilateral elements, it s a Ne-by-4 matrix
%
%        Ne is the number of elements in the mesh. Each row in ele 
%        contains the indices of the nodes for that mesh element.
%   
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%         When omitted, assign one phase.
%
%   color_code: Color code for selecting colormap.
%               Interger. Value: 0-10. Default value: 0
%			    0: grayscale, 1: lines, 2: parula, 3: turbo, 4: jet, 5: hot
%			    6: cool, 7: summer, 8: winter, 9: bone, 10: pink.
%
%   opt: a structure array. It is the extra options for plotMeshes.
%        It stores extra parameter settings for plotMeshes.
%
%   opt.mode: display mode. Value: 1 or 2.
% 			  When opt.mode=1, plot faces and edges. Slower.
% 			  When opt.mode=2, plot edges only. Faster.
%             Default value: 1
%
%	opt.wid: line width of the plotted edges. Positive value.
%            Default value: 0.5
%
%   opt.alpha: edge line transparency. It's a scalar value in range [0,1].
%              When the mesh data is huge, I usually set alpha to 0.2-0.3.
%              Default value: 0.5
%
%   opt.beta: brightness adjustment of colormap. Scalar value in range 
%             [-1, 1]. The colors brighten when beta >0. The colors darken 
%             when beta <0. The magnitude of the color change is 
%             proportional to the magnitude of beta.
%             Default value: 0
%
%   opt.tf_gs: Boolean value. Whether to use graphics smoothing when
%              plotting. For some computer hardwares, set opt.tf_gs to 0
%              may plot mesh faster.
%              Default value: 1
%
%   opt.q_thresh: threshold for mesh quality. Those elements with mesh 
%                 quality lower than opt.q_thresh will be highlighted. You
%                 can set color_code to 0 to faciliate the visualization.
%                 Default value: 0
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    %--------------------------------------------------------------------   
    % Check the number of inputs. If missing, set as empty. 
    if nargin < 2
        error("Not enough input arguments.");
    end
    
    if nargin < 3
        tnum = [];
    end
    
    if nargin < 4
        color_code = [];
    end
    
    if nargin < 5
        opt = [];
    end

    % ---------------------------------------------------------------------
    % verify field names and set values for opt
    opt = setOption( opt );

    % ---------------------------------------------------------------------
    % If input is empty, assign defaualt value to input
    if isempty(tnum)
        tnum = ones( size(ele,1), 1 );
    end

    if isempty(color_code)
        color_code = 0;
    end

    %--------------------------------------------------------------------
    % check element type
    ele_wid = size(ele,2);

    if ele_wid == 3         % linear triangle
        range_vec = 1:3;
    elseif ele_wid == 6     % quadratic triangle
        range_vec = [1 4 2 5 3 6];
    elseif ele_wid == 4     % linear quadrilateral
        range_vec = 1:4;
    elseif ele_wid == 8     % quadratic quadrilateral
        range_vec = [1 5 2 6 3 7 4 8];
    else
        error("ele - wrong size.")
    end

    %--------------------------------------------------------------------
    tvalue = unique( tnum );
    num_phase = length( tvalue );
    
    %--------------------------------------------------------------------
    % set color
    % Create variable 'colors' - num_phase-by-3 array.
    % Each row in 'colors' is one rgb color.

    switch color_code
        case 0
            % grayscale
            if num_phase == 1
                col = 0.98;
            elseif num_phase > 1
                col = 0.3: 0.68/(num_phase-1): 0.98;
                col = col(:);
            end
            colors = [col, col, col];

        case 1
            colors = lines( num_phase );
        case 2
            colors = parula( num_phase );
        case 3
            colors = turbo( num_phase );
        case 4
            colors = jet( num_phase );
        case 5
            colors = hot( num_phase );
        case 6
            colors = cool( num_phase );
        case 7
            colors = summer( num_phase );
        case 8
            colors = winter( num_phase );
        case 9
            colors = bone( num_phase );
        case 10
            colors = pink( num_phase );
        otherwise
            error('Input argument color_code is out of range.')
    end

    colors = brighten( colors, opt.beta );

    %--------------------------------------------------------------------
    % plot mesh
    figure;
    hold on;
    axis image off;

    if opt.tf_gs == 0
        set(gcf,'GraphicsSmoothing','off');
    end
    
    if opt.mode == 1
        % use function patch to plot faces and edges
        for i = 1: num_phase
            current_phase = tvalue(i);
            patch( ...
                'faces', ele( tnum==current_phase, range_vec ), ...
                'vertices', vert, ...
                'facecolor', colors(i,:), ...
                'edgecolor', [.1,.1,.1], ...
                'linewidth', opt.wid, ...
                'edgealpha', opt.alpha ...
                );
        end
    elseif opt.mode == 2
        % use function patch to plot edges
        for i = 1: num_phase
            current_phase = tvalue(i);
            patch( ...
                'faces', ele( tnum==current_phase, range_vec ), ...
                'vertices', vert, ...
                'facecolor', 'none', ...
                'edgecolor', 0.8*colors(i,:), ...
                'linewidth', opt.wid, ...
                'edgealpha', opt.alpha ...
                );
        end
    else
        error('opt.mode has to be 1 or 2.')
    end
    
    if opt.q_thresh > 0 && opt.q_thresh < 1
        % evaluate mesh quality
        if ele_wid == 3 || ele_wid == 6
            q = triscr2( vert, ele );
        elseif ele_wid == 4 || ele_wid == 8
            q = MeshQualityQuads( ele(:,1:4), vert(:,1), vert(:,2) );
        end
    
        % highlight element with quality lower than opt.q_thresh
        tf_bad_vec = q < opt.q_thresh;
        if any( tf_bad_vec )
		    if opt.mode == 1
			    patch( ...
				    'faces', ele( tf_bad_vec, range_vec ), ...
				    'vertices', vert, ...
				    'facecolor', 'r', ...
				    'edgecolor', 'g', ...
				    'linewidth', 3, ...
				    'edgealpha', 1 ...
				    );
		    elseif opt.mode == 2
			    patch( ...
				    'faces', ele( tf_bad_vec, range_vec ), ...
				    'vertices', vert, ...
				    'facecolor', 'none', ...
				    'edgecolor', [0.8500 0.3250 0.0980], ...
				    'linewidth', 3, ...
				    'edgealpha', 1 ...
				    );
		    end
        end
    end

    hold off
    %--------------------------------------------------------------------
    
end


function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    % initialize new_opt with default field names & value 
    new_opt.mode = 1;
    new_opt.wid = 0.5;
    new_opt.alpha = 0.5;
    new_opt.beta = 0;
    new_opt.tf_gs = 1;
    new_opt.q_thresh = 0;
    
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



