function objs = regionpropsOS(Label, varargin)
% regionpropsOS: Measure properties of image regions
% Open-source version of function regionprops
%
%   STATS = REGIONPROPSOS(L,PROPERTIES) measures properties for each
%   labeled region in the label matrix L. Positive integer elements of L
%   correspond to different regions. For example, the set of elements of L
%   equal to 1 corresponds to region 1; the set of elements of L equal to 2
%   corresponds to region 2; and so on.
%
%   Note: This open-source implementation is specifically optimized for
%   extracting 'Image' and 'BoundingBox' properties without requiring the
%   Image Processing Toolbox. Other properties passed in the PROPERTIES
%   argument are currently ignored to maximize execution speed.
%       L = bwlabel(BW, 4);
%       stats = regionpropsOS(L, 'Image', 'BoundingBox');
%
%   Property Definitions
%   --------------------
%   'BoundingBox' - A 1-by-4 vector [x y width height] specifying the
%                   smallest bounding box containing the region. x and y
%                   are the spatial coordinates of the upper-left corner
%                   of the bounding box.
%
%   'Image'       - A binary image (logical) of the same size as the
%                   bounding box of the region. The "on" pixels correspond 
%                   to the region, and all other pixels are "off".
%
%   Class Support
%   -------------
%   L must be a numeric array containing positive integers. L can have any
%   numeric class. Negative-valued pixels or zeros are treated as 
%   background.
%
%   Example
%   -------
%       % Label a synthetic image and extract properties
%       BW = logical([1 1 0 0; 1 1 0 0; 0 0 1 1; 0 0 1 1]);
%       L = bwlabel(BW, 4);     % Assuming bwlabel is available
%       stats = regionpropsOS(L, 'Image', 'BoundingBox');
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh


    % Determine the total number of connected components (objects)
    % by finding the maximum integer label in the matrix.
    num_objs = max(Label(:));
    
    % Initialize the output structure array. 
    % Preallocating the array using repmat prevents dynamic memory allocation 
    % during the loop, significantly improving execution speed for large numbers 
    % of objects.
    objs = repmat(struct('Image', [], 'BoundingBox', []), num_objs, 1);
    
    % Iterate over each labeled object
    for k = 1:num_objs
        
        % Find the linear indices of all pixels belonging to the k-th object,
        % then convert them to row (r) and column (c) subscripts.
        [r, c] = find(Label == k);
        
        % Handle potential edge cases where a label index might be missing
        % (e.g., if the matrix contains labels 1, 2, and 4, but no 3).
        if isempty(r)
            objs(k).BoundingBox = [0.5, 0.5, 0, 0];
            objs(k).Image = logical([]);
            continue;
        end
        
        % Identify the extreme spatial coordinates of the object.
        minR = min(r);
        maxR = max(r);
        minC = min(c);
        maxC = max(c);
        
        % -----------------------------------------------------------------
        % Compute BoundingBox
        
        % Calculate the spatial dimensions of the bounding box.
        width = maxC - minC + 1;
        height = maxR - minR + 1;
        
        % The MATLAB standard for bounding boxes uses spatial coordinates
        % where the center of the first pixel is (1,1) and its upper-left
        % corner is (0.5, 0.5). Therefore, subtract 0.5 from the minimum
        % subscripts to align with MATLAB's default spatial coordinate system.
        % Format: [min_col - 0.5, min_row - 0.5, width, height]
        objs(k).BoundingBox = [minC - 0.5, minR - 0.5, width, height];
        
        % -----------------------------------------------------------------
        % Compute Image
        
        % Extract the submatrix from the original label matrix that is 
        % completely enclosed by the bounding box limits.
        subImage = Label(minR:maxR, minC:maxC);
        
        % Create a tight logical mask for the object. 
        % We use '== k' instead of just '> 0' to ensure that if a different, 
        % discontiguous object encroaches into this bounding box space, its 
        % pixels are correctly ignored (evaluated as false/0).
        objs(k).Image = (subImage == k);
        % -----------------------------------------------------------------
    end

end