function [L, numComponents] = bwlabelOS(BW, mode)
% bwlabelOS: Label connected components in 2-D binary image 
% Open-source version of function bwlabel
%
%   L = bwlabelOS(BW,N) returns a matrix L, of the same size as BW,
%   containing labels for the connected components in BW. N can have a
%   value of either 4 or 8, where 4 specifies 4-connected objects and 8
%   specifies 8-connected objects; if the argument is omitted, it defaults
%   to 8.
%
%   The elements of L are integer values greater than or equal to 0.  The
%   pixels labeled 0 are the background.  The pixels labeled 1 make up one
%   object, the pixels labeled 2 make up a second object, and so on.
%
%   [L,NUM] = bwlabelOS(BW,N) returns in NUM the number of connected objects
%   found in BW.
%
%   ------------------------------------------------------------------
%   To extract features from a binary image using REGIONPROPS using the
%   default connectivity, just pass BW directly into REGIONPROPS, i.e.,
%   REGIONPROPS(BW).
%
%   To compute a label matrix having a more memory-efficient data type
%   (e.g., uint8 versus double), use the LABELMATRIX function on the output
%   of BWCONNCOMP. 
%
%   Class Support
%   -------------
%   BW can be logical or numeric, and it must be real, 2-D, and nonsparse.
%   L is of class double.
%
%   Example
%   -------
%       BW = logical([1 1 1 0 0 0 0 0
%                     1 1 1 0 1 1 0 0
%                     1 1 1 0 1 1 0 0
%                     1 1 1 0 0 0 1 0
%                     1 1 1 0 0 0 1 0
%                     1 1 1 0 0 0 1 0
%                     1 1 1 0 0 1 1 0
%                     1 1 1 0 0 0 0 0]);
%       L = bwlabelOS(BW,4)
%       [r,c] = find(L == 2)
%
%   Notes
%   -------
%   To mimic the speed of compiled internal functions, this function uses 
%   a Breadth-First Search (BFS) algorithm combined with a zero-padding 
%   technique. By padding the image with zeros, we can compute linear 
%   indices for 4-connected or 8-connected neighbors without writing slow, 
%   repetitive boundary checks (e.g., checking if a pixel is on the edge 
%   of the image). It also uses vectorized logical indexing to evaluate all
%   neighbors simultaneously.
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh

    % 1. Input Validation
    if nargin < 2
        mode = 8;
    end
    
    if mode ~= 4 && mode ~= 8
        error('Mode must be either 4 or 8.');
    end
    
    % Ensure input is logical
    BW = BW ~= 0; 
    
    [rows, cols] = size(BW);
    if rows == 0 || cols == 0
        L = zeros(rows, cols);
        numComponents = 0;
        return;
    end
    
    % 2. Pad the image with zeros 
    % This is a crucial optimization. It ensures we can look at adjacent 
    % pixels without running into boundary limits, removing the need for 
    % slow "if x > 1" edge checks in the loop.
    padRows = rows + 2;
    padCols = cols + 2;
    BW_pad = false(padRows, padCols);
    BW_pad(2:end-1, 2:end-1) = BW;
    
    L_pad = zeros(padRows, padCols); % Output label matrix (padded)
    
    % 3. Precalculate linear index offsets for neighbors
    % MATLAB arrays are column-major. Moving down a row is +1, 
    % moving right a column is +padRows.
    if mode == 4
        offsets = [-1, 1, -padRows, padRows];
    else % mode == 8
        offsets = [-1, 1, -padRows, padRows, ...
                   -padRows-1, -padRows+1, padRows-1, padRows+1];
    end
    
    % 4. Initialization
    numComponents = 0;
    linearIndices = find(BW_pad); % Find all foreground pixels
    
    % Preallocate a queue for the Breadth-First Search (BFS)
    % Max possible size is the number of 1s in the image
    queue = zeros(length(linearIndices), 1); 
    
    % 5. Core BFS loop
    for k = 1:length(linearIndices)
        idx = linearIndices(k);
        
        % If this pixel hasn't been labeled yet, it's a new component
        if L_pad(idx) == 0 
            numComponents = numComponents + 1;
            
            % Initialize BFS queue
            head = 1;
            tail = 1;
            queue(tail) = idx;
            L_pad(idx) = numComponents;
            
            % Traverse the current connected component
            while head <= tail
                curr = queue(head);
                head = head + 1;
                
                % Calculate linear indices of all neighbors
                neighbors = curr + offsets;
                
                % Vectorized check: Find neighbors that are foreground (1) 
                % AND unvisited (0)
                valid_idx = BW_pad(neighbors) & (L_pad(neighbors) == 0);
                valid_nbrs = neighbors(valid_idx);
                
                num_valid = length(valid_nbrs);
                if num_valid > 0
                    % Label them immediately to prevent re-queueing
                    L_pad(valid_nbrs) = numComponents;
                    
                    % Add valid neighbors to the end of the queue
                    queue(tail+1 : tail+num_valid) = valid_nbrs;
                    tail = tail + num_valid;
                end
            end
        end
    end
    
    % 6. Extract the unpadded result
    L = L_pad(2:end-1, 2:end-1);
end