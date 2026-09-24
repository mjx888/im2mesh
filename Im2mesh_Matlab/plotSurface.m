function plotSurface( V, faces, color_code, opt )
% plotSurface: plot surfaces
% The second argument has to be an N-by-1 cell array.
%
% usage 1:
%   plotSurface( V, faces );
%
% usage 2:
%   color_code = 2;
%   plotSurface( V, faces, color_code )
%
% usage 3:
%     color_code = 4;
%     opt = [];   % reset
%     opt.faceAlpha = 0.5;
%     opt.sampleRatio = 0.1;
%     
%     plotSurface( V, phaseFaces, color_code, opt );
%
% input:
%   Argument color_code, and opt are optional.
%
%   V: Mesh nodes. It's a Nn-by-3 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y, z coordinates for that mesh node.
%     
%   faces: N-by-1 cell array. N is the number of phases.
%          faces{i} is Ne-by-3 array. Define triangular surface in one
%          phase.
%
%   color_code: Color code for selecting colormap.
%                Interger. Value: 0-10. Default value: 0
%                0: grayscale, 1: lines, 2: parula, 3: turbo, 4: jet, 5: hot
%                6: cool, 7: summer, 8: winter, 9: bone, 10: pink.
%
%   opt: a structure array. It is the extra options for plotSurface.
%        It stores extra parameter settings for plotSurface.
%
%   opt.faceAlpha: face transparency. It's a scalar value in range [0,1].
%              Default value: 1
%
%   opt.beta: brightness adjustment of colormap. Scalar value in range 
%             [-1, 1]. The colors brighten when beta >0. The colors darken 
%             when beta <0. The magnitude of the color change is 
%             proportional to the magnitude of beta.
%             Default value: 0
%
%   opt.sampleRatio: down-sampling ratio. 
%         'opt.sampleRatio = 1' means no down-sampling.
%         'opt.sampleRatio = 0.1' means roundly 10% down-sampling before 
%           plotting surface.
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
        color_code = [];
    end
    
    if nargin < 4
        opt = [];
    end

    % ---------------------------------------------------------------------
    % verify field names and set values for opt
    opt = setOption( opt );

    % ---------------------------------------------------------------------
    % If input is empty, assign defaualt value to input
    if isempty(color_code)
        color_code = 0;
    end

    %--------------------------------------------------------------------
    % check element type 
    ele_wid = size(V,2);

    if ele_wid == 3        
        % do nothing
    else
        error("V - wrong size. Must be an N-by-3 array.")
    end

    %--------------------------------------------------------------------
    num_phase = length( faces );
    
    %--------------------------------------------------------------------
    % set color
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
    % GLOBAL DOWNSAMPLING
    % Run the expensive reduction ONCE on the entire dataset
    if opt.sampleRatio < 1
        [faces, V] = global_fast_reduce(faces, V, opt.sampleRatio);
    end

    %--------------------------------------------------------------------
    % FAST PRE-ALLOCATION
    % Calculate the exact total number of faces to pre-allocate memory
    total_faces = sum(cellfun(@(x) size(x, 1), faces));
    
    all_F = zeros(total_faces, 3);
    all_C = zeros(total_faces, 3);
    
    % Populate the unified arrays linearly
    current_idx = 1;
    for i = 1: num_phase
        F_i = faces{i};
        num_f = size(F_i, 1);
        
        if num_f > 0
            all_F(current_idx : current_idx + num_f - 1, :) = F_i;
            all_C(current_idx : current_idx + num_f - 1, :) = repmat(colors(i, :), num_f, 1);
            current_idx = current_idx + num_f;
        end
    end

    %--------------------------------------------------------------------
    % plot mesh
    figure('GraphicsSmoothing', 'on', 'Renderer', 'opengl');
    hold on;
    axis image off;
    % camlight headlight;
    
    % Single patch using the unified, globally reduced vertices
    patch( ...
        'Faces', all_F, ...
        'Vertices', V, ...
        'FaceVertexCData', all_C, ... 
        'FaceColor', 'flat', ...      
        'FaceAlpha', opt.faceAlpha, ...
        'EdgeColor', 'none', ... 
        'FaceLighting', 'flat', ...       
        'BackFaceLighting', 'unlit', ...  
        'AmbientStrength', 0.5);

    view([45 30]);
    hold off
    %--------------------------------------------------------------------
    
end

function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    % initialize new_opt with default field names & value 
    new_opt.faceAlpha = 1;
    new_opt.beta = 0;
    new_opt.sampleRatio = 1;
    
    if isempty(opt)
        return
    end

    if ~isstruct(opt)
        error("opt is not a structure array. Not valid input.")
    end

    nameC = fieldnames(opt);

    for i = 1: length(nameC)
        if isfield( new_opt, nameC{i} )
            new_opt.(nameC{i}) = opt.(nameC{i});
        else
            error("Field name %s in opt is not correct.", nameC{i});
        end
    end
end

function [faces_cell, V] = global_fast_reduce(faces_cell, V, target_ratio)
% global_fast_reduce decimates a multi-phase mesh processing vertices only once
    
    target_vertices = size(V, 1) * target_ratio;
    grid_bins = max(2, round(sqrt(target_vertices / 6))); 
    
    min_V = min(V, [], 1);
    max_V = max(V, [], 1);
    max_range = max(max_V - min_V);
    
    if max_range == 0
        max_range = 1; 
    end
    
    % Quantize globally
    V_quantized = round(((V - min_V) ./ max_range) * grid_bins);
    
    % Find unique vertices ONE TIME for the entire mesh
    [~, unique_idx, vertex_map] = unique(V_quantized, 'rows');
    V = V(unique_idx, :);
    
    % Remap the faces for each phase
    for i = 1:length(faces_cell)
        F = faces_cell{i};
        if isempty(F)
            continue;
        end
        
        mapped_F = vertex_map(F);
        
        % Remove degenerate faces (triangles crushed into lines or points)
        valid_faces = (mapped_F(:,1) ~= mapped_F(:,2)) & ...
                      (mapped_F(:,2) ~= mapped_F(:,3)) & ...
                      (mapped_F(:,1) ~= mapped_F(:,3));
                  
        faces_cell{i} = mapped_F(valid_faces, :);
    end
end