function B = bwboundariesOS(BW, conn)
% bwboundariesOS Trace boundaries in a binary image
% Open-source version of function bwboundaries
%
%   B = bwboundariesOS(BW, 8) traces the exterior boundary of a single
%   binary object in BW, as well as the boundaries of holes inside that
%   object. Nonzero pixels in BW are treated as object pixels, and zero
%   pixels are treated as background pixels.
%
%   B is a P-by-1 cell array, where P is the number of boundaries found.
%   The first cell, B{1}, contains the exterior object boundary. The
%   remaining cells, B{2}, B{3}, ..., contain hole boundaries if holes are
%   present.
%
%   Each cell contains a Q-by-2 matrix. Each row of the matrix contains the
%   row and column coordinates of one boundary pixel:
%
%       [row, column]
%
%   This function is designed as a lightweight open-source replacement for
%   the restricted use case:
%
%       B = bwboundaries(BW, 8);
%
%   Important assumptions:
%       1. BW is a real, nonsparse, 2-D numeric or logical array.
%       2. BW contains one foreground object.
%       3. The connectivity input conn is always 8.
%       4. Only the first output B is needed.
%       5. Label matrix, number of objects, and adjacency matrix outputs
%          are not implemented.
%
%   Connectivity behavior:
%       When conn = 8, object pixels are traced using 8-connectivity.
%       To match MATLAB's topology convention, holes are detected and
%       traced using the opposite background connectivity, namely
%       4-connectivity.
%
%   Class Support
%   -------------
%   BW can be logical or numeric. It must be real, 2-D, and nonsparse.
%   The output B is a cell array containing double coordinate matrices.
%
%   Example
%   -------
%       bw = [
%          0 0 0 0 0 0
%          0 0 1 1 1 0
%          0 1 1 1 1 0
%          0 1 1 1 1 0
%          0 1 0 0 1 0
%          0 1 0 0 1 0
%          0 0 1 1 1 0
%          0 0 0 0 0 0 ];
%
%       bw = logical(bw);
%       B = bwboundariesOS(bw, 8);
%
%       figure;
%       imagesc(bw);
%       axis image;
%       colormap gray;
%       hold on;
%
%       for k = 1:length(B)
%           boundary = B{k};
%           plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
%       end
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh


    % ---------------------------------------------------------------------
    % Input parsing and validation
    % ---------------------------------------------------------------------

    if nargin < 2
        conn = 8;
    end

    % This implementation is intentionally limited to the user's target
    % workflow, where bwboundaries is always called with conn = 8.
    if conn ~= 8
        error('bwboundariesOS:OnlyConn8', ...
              'This simplified version only supports conn = 8.');
    end

    % Match the basic input requirements of MATLAB bwboundaries:
    % numeric or logical, real, 2-D, nonsparse.
    if ~(isnumeric(BW) || islogical(BW)) || ~isreal(BW) || ...
            ~ismatrix(BW) || issparse(BW)
        error('bwboundariesOS:InvalidInput', ...
              'BW must be a real, nonsparse, 2-D numeric or logical array.');
    end

    % Convert input to logical.
    %
    % This matches MATLAB behavior:
    %   nonzero values are treated as foreground/object pixels.
    %
    % Note:
    %   NaN ~= 0 returns true, so NaN values are treated as foreground,
    %   consistent with MATLAB bwboundaries behavior.
    BW = (BW ~= 0);

    % Initialize output as an empty cell array.
    B = cell(0,1);

    % If there are no foreground pixels, there are no object boundaries.
    if ~any(BW(:))
        return;
    end

    % ---------------------------------------------------------------------
    % Step 1: Trace the exterior object boundary
    % ---------------------------------------------------------------------
    %
    % In bwboundaries(BW, 8), foreground object pixels are considered
    % 8-connected. Therefore, the exterior boundary is traced using
    % an 8-neighbor Moore boundary tracing algorithm.
    %
    % The resulting boundary is stored as B{1}.
    B{end+1,1} = localTraceBoundary8(BW);

    % ---------------------------------------------------------------------
    % Step 2: Detect holes inside the object
    % ---------------------------------------------------------------------
    %
    % MATLAB uses complementary connectivity to avoid topological ambiguity.
    %
    % If foreground objects use 8-connectivity, then background holes are
    % defined using 4-connectivity.
    %
    % The function localFindHoles4 finds background components that are not
    % connected to the image border by 4-connectivity. These components are
    % holes inside the foreground object.
    Holes = localFindHoles4(BW);

    % ---------------------------------------------------------------------
    % Step 3: Trace hole boundaries
    % ---------------------------------------------------------------------
    %
    % Each connected hole component is labeled using 4-connectivity.
    % Then each hole boundary is traced and appended to B.
    %
    % Therefore:
    %   B{1}     = exterior object boundary
    %   B{2:end} = hole boundaries
    if any(Holes(:))
        [Lh, numHoles] = localLabelComponents(Holes, 4);

        for k = 1:numHoles
            holeMask = (Lh == k);

            % Important:
            % Hole boundaries are traced with 4-connected behavior.
            % This avoids diagonal shortcuts across hole corners and better
            % matches MATLAB bwboundaries(BW, 8).
            B{end+1,1} = localTraceBoundary4(holeMask);
        end
    end
end


% -------------------------------------------------------------------------
function boundary = localTraceBoundary8(M)
%LOCALTRACEBOUNDARY8 Trace one 8-connected foreground boundary.
%
%   boundary = localTraceBoundary8(M) traces the boundary of one connected
%   foreground component in the logical mask M.
%
%   The returned boundary is a Q-by-2 matrix of row-column coordinates.
%
%   This function implements a Moore-neighbor tracing procedure. At each
%   boundary pixel, the algorithm searches the eight neighboring pixels in
%   a fixed clockwise order to find the next boundary pixel.
%
%   Neighbor search order:
%
%       W, NW, N, NE, E, SE, S, SW
%
%   This order was chosen to match MATLAB bwboundaries(BW, 8) behavior for
%   the intended single-object use case.

    [m, n] = size(M);

    % Find the first foreground pixel in column-major order, matching
    % MATLAB's find behavior.
    idx = find(M, 1, 'first');

    if isempty(idx)
        boundary = zeros(0,2);
        return;
    end

    [startR, startC] = ind2sub([m, n], idx);

    % Eight neighbor offsets.
    %
    % The order is important. It controls the direction of tracing and the
    % sequence of output boundary coordinates.
    dirs = [
         0, -1;   % W
        -1, -1;   % NW
        -1,  0;   % N
        -1,  1;   % NE
         0,  1;   % E
         1,  1;   % SE
         1,  0;   % S
         1, -1    % SW
    ];

    % ---------------------------------------------------------------------
    % Single-pixel edge case
    % ---------------------------------------------------------------------
    %
    % MATLAB bwboundaries returns a closed boundary even for a single
    % isolated pixel. Therefore, the coordinate is repeated:
    %
    %       [r c
    %        r c]
    %
    % Without this special case, the output would contain only one point
    % and would not match native bwboundaries.
    if ~localHasForegroundNeighbor(M, startR, startC, dirs)
        boundary = [
            startR, startC;
            startR, startC
        ];
        return;
    end

    % Upper bound on the number of tracing steps.
    %
    % For normal objects, the boundary length is much smaller than this.
    % The guard prevents infinite loops if an unexpected pathological case
    % is encountered.
    maxSteps = max(100, 8 * nnz(M) + 8);

    boundary = zeros(maxSteps + 1, 2);

    % Current boundary pixel.
    pR = startR;
    pC = startC;

    % Backtrack pixel.
    %
    % The initial backtrack pixel is chosen as the pixel immediately west
    % of the starting pixel. This defines the initial search direction.
    bR = startR;
    bC = startC - 1;

    % Store the starting boundary pixel.
    boundary(1,:) = [pR, pC];
    count = 1;

    for step = 1:maxSteps

        % Determine where the backtrack pixel lies relative to the current
        % boundary pixel.
        rel = [bR - pR, bC - pC];

        dirIdx = find(dirs(:,1) == rel(1) & dirs(:,2) == rel(2), 1);

        if isempty(dirIdx)
            error('bwboundariesOS:TraceError', ...
                  'Internal tracing error: invalid backtrack direction.');
        end

        found = false;

        % Search clockwise starting from the backtrack direction.
        %
        % The first foreground neighbor encountered becomes the next
        % boundary pixel.
        for offset = 0:7
            j = mod(dirIdx - 1 + offset, 8) + 1;

            qR = pR + dirs(j,1);
            qC = pC + dirs(j,2);

            if qR >= 1 && qR <= m && qC >= 1 && qC <= n && M(qR, qC)
                found = true;
                foundIdx = j;
                nextR = qR;
                nextC = qC;
                break;
            end
        end

        if ~found
            error('bwboundariesOS:TraceError', ...
                  'Unable to find next boundary pixel.');
        end

        % Update the backtrack pixel.
        %
        % The new backtrack pixel is the neighbor immediately before the
        % found foreground neighbor in the clockwise search order.
        prevIdx = mod(foundIdx - 2, 8) + 1;

        bR = pR + dirs(prevIdx,1);
        bC = pC + dirs(prevIdx,2);

        % Advance to the next boundary pixel.
        pR = nextR;
        pC = nextC;

        count = count + 1;
        boundary(count,:) = [pR, pC];

        % Stop once the trace returns to the starting pixel.
        %
        % The returned boundary is explicitly closed because the starting
        % coordinate appears again as the final coordinate.
        if pR == startR && pC == startC
            break;
        end
    end

    if count >= maxSteps
        error('bwboundariesOS:TraceError', ...
              'Boundary tracing did not close.');
    end

    % Trim unused preallocated rows.
    boundary = boundary(1:count,:);
end


% -------------------------------------------------------------------------
function boundary4 = localTraceBoundary4(M)
%LOCALTRACEBOUNDARY4 Trace a 4-connected hole boundary.
%
%   boundary4 = localTraceBoundary4(M) traces the boundary of one
%   4-connected hole component.
%
%   In MATLAB bwboundaries(BW, 8), foreground objects are 8-connected, but
%   holes are treated as 4-connected background components. This avoids
%   topological ambiguity between object and background connectivity.
%
%   Implementation strategy:
%
%       1. Trace the hole component using the 8-neighbor tracer.
%       2. Inspect each consecutive pair of boundary pixels.
%       3. If a diagonal step occurs, insert an intermediate pixel so that
%          the resulting boundary follows a 4-connected path.
%
%   This prevents diagonal shortcuts across hole corners.

    boundary8 = localTraceBoundary8(M);

    if isempty(boundary8)
        boundary4 = boundary8;
        return;
    end

    % Single-pixel hole edge case.
    %
    % localTraceBoundary8 already returns:
    %
    %       [r c
    %        r c]
    %
    % which is also the desired 4-connected closed boundary.
    if size(boundary8,1) == 2 && isequal(boundary8(1,:), boundary8(2,:))
        boundary4 = boundary8;
        return;
    end

    if size(boundary8,1) <= 1
        boundary4 = boundary8;
        return;
    end

    % The 4-connected trace may insert at most one extra point for each
    % diagonal step, so twice the 8-connected boundary length is sufficient.
    maxLen = 2 * size(boundary8,1) + 5;
    boundary4 = zeros(maxLen, 2);

    count = 1;
    boundary4(count,:) = boundary8(1,:);

    for k = 1:size(boundary8,1)-1

        p = boundary8(k,:);
        q = boundary8(k+1,:);

        dr = q(1) - p(1);
        dc = q(2) - p(2);

        % -----------------------------------------------------------------
        % Convert diagonal moves into 4-connected moves
        % -----------------------------------------------------------------
        %
        % A diagonal move from p to q changes both row and column by one:
        %
        %       abs(dr) == 1 and abs(dc) == 1
        %
        % For 4-connectivity, this must be replaced by two axis-aligned
        % moves through an intermediate pixel.
        if abs(dr) == 1 && abs(dc) == 1

            % Candidate intermediate pixels:
            %
            %   cand1: move in row first, then column
            %   cand2: move in column first, then row
            cand1 = [q(1), p(2)];
            cand2 = [p(1), q(2)];

            % Prefer the candidate that belongs to the hole mask.
            in1 = localPixelInMask(M, cand1(1), cand1(2));
            in2 = localPixelInMask(M, cand2(1), cand2(2));

            if in1 && ~in2
                mid = cand1;
            elseif in2 && ~in1
                mid = cand2;
            elseif in1 && in2
                % Ambiguous 2-by-2 case.
                %
                % Both candidates are valid foreground pixels of the hole.
                % This tie-break rule gives stable MATLAB-like ordering.
                if dr * dc > 0
                    mid = cand1;
                else
                    mid = cand2;
                end
            else
                % Rare fallback.
                %
                % If neither candidate is inside the mask, do not insert an
                % intermediate point. This should rarely occur for a true
                % 4-connected hole component.
                mid = [];
            end

            % Append intermediate point if it is valid and not a duplicate
            % of the most recently stored point.
            if ~isempty(mid)
                if ~isequal(boundary4(count,:), mid)
                    count = count + 1;
                    boundary4(count,:) = mid;
                end
            end
        end

        % Append the next boundary point if it is not a duplicate.
        if ~isequal(boundary4(count,:), q)
            count = count + 1;
            boundary4(count,:) = q;
        end
    end

    % Trim unused preallocated rows.
    boundary4 = boundary4(1:count,:);

    % Ensure explicit closure.
    %
    % MATLAB bwboundaries returns closed coordinate sequences, meaning that
    % the final boundary coordinate equals the first coordinate.
    if size(boundary4,1) > 1 && ~isequal(boundary4(1,:), boundary4(end,:))
        boundary4(end+1,:) = boundary4(1,:);
    end
end


% -------------------------------------------------------------------------
function tf = localPixelInMask(M, r, c)
%LOCALPIXELINMASK True if pixel location is inside image and foreground.
%
%   tf = localPixelInMask(M, r, c) returns true if the row-column location
%   (r,c) lies within M and M(r,c) is true.

    [m, n] = size(M);

    tf = ...
        r >= 1 && r <= m && ...
        c >= 1 && c <= n && ...
        M(r,c);
end


% -------------------------------------------------------------------------
function tf = localHasForegroundNeighbor(M, r, c, dirs)
%LOCALHASFOREGROUNDNEIGHBOR Test whether a pixel has a foreground neighbor.
%
%   tf = localHasForegroundNeighbor(M, r, c, dirs) returns true if the
%   pixel at row r, column c has at least one foreground neighbor using the
%   neighbor offsets specified by dirs.

    [m, n] = size(M);
    tf = false;

    for k = 1:size(dirs,1)
        rr = r + dirs(k,1);
        cc = c + dirs(k,2);

        if rr >= 1 && rr <= m && cc >= 1 && cc <= n && M(rr, cc)
            tf = true;
            return;
        end
    end
end


% -------------------------------------------------------------------------
function Holes = localFindHoles4(BW)
%LOCALFINDHOLES4 Find 4-connected background holes.
%
%   Holes = localFindHoles4(BW) returns a logical mask containing
%   background components that are not connected to the image border by
%   4-connectivity.
%
%   This is equivalent in purpose to:
%
%       BWcomplement = imcomplement(BW);
%       BWholes = imclearborder(BWcomplement, 4);
%
%   but it does not require Image Processing Toolbox.
%
%   Method:
%       1. Compute the background mask BG = ~BW.
%       2. Flood-fill all background pixels connected to the image border
%          using 4-connectivity.
%       3. Any remaining background pixels are holes.

    [m, n] = size(BW);

    BG = ~BW;
    outside = false(m, n);

    % Queue arrays for flood fill.
    %
    % qR and qC store row and column coordinates of background pixels that
    % still need to be processed.
    qR = zeros(numel(BW), 1);
    qC = zeros(numel(BW), 1);

    head = 1;
    tail = 0;

    % Collect all image border coordinates.
    %
    % Duplicate corner coordinates are harmless because each pixel is added
    % to the queue only once.
    border = [
        ones(n,1),       (1:n)';
        m * ones(n,1),   (1:n)';
        (1:m)',          ones(m,1);
        (1:m)',          n * ones(m,1)
    ];

    % Initialize queue with all border background pixels.
    for k = 1:size(border,1)
        r = border(k,1);
        c = border(k,2);

        if BG(r,c) && ~outside(r,c)
            tail = tail + 1;
            qR(tail) = r;
            qC(tail) = c;
            outside(r,c) = true;
        end
    end

    % Four-neighbor offsets.
    nbrs = [
        -1,  0;
         1,  0;
         0, -1;
         0,  1
    ];

    % Flood-fill the outside background.
    while head <= tail
        r = qR(head);
        c = qC(head);
        head = head + 1;

        for k = 1:4
            rr = r + nbrs(k,1);
            cc = c + nbrs(k,2);

            if rr >= 1 && rr <= m && cc >= 1 && cc <= n
                if BG(rr,cc) && ~outside(rr,cc)
                    tail = tail + 1;
                    qR(tail) = rr;
                    qC(tail) = cc;
                    outside(rr,cc) = true;
                end
            end
        end
    end

    % Background pixels that are not connected to the border are holes.
    Holes = BG & ~outside;
end


% -------------------------------------------------------------------------
function [L, num] = localLabelComponents(M, conn)
%LOCALLABELCOMPONENTS Label connected components in a binary mask.
%
%   [L, num] = localLabelComponents(M, conn) labels connected foreground
%   components in the logical mask M.
%
%   conn can be either 4 or 8.
%
%   L is a label matrix of the same size as M. Background pixels are zero.
%   The k-th connected component has label k.
%
%   num is the number of connected components.
%
%   This helper is a lightweight replacement for bwlabel in the limited
%   internal use case of labeling hole components.

    [m, n] = size(M);

    L = zeros(m, n);
    num = 0;

    % Select neighbor offsets based on requested connectivity.
    if conn == 4
        nbrs = [
            -1,  0;
             1,  0;
             0, -1;
             0,  1
        ];
    elseif conn == 8
        nbrs = [
            -1, -1;
            -1,  0;
            -1,  1;
             0, -1;
             0,  1;
             1, -1;
             1,  0;
             1,  1
        ];
    else
        error('bwboundariesOS:BadConnectivity', ...
              'Connectivity must be 4 or 8.');
    end

    % Queue arrays used for breadth-first component labeling.
    qR = zeros(numel(M), 1);
    qC = zeros(numel(M), 1);

    % Scan image in column-major order, matching MATLAB's linear indexing.
    for c = 1:n
        for r = 1:m

            % A new component starts at any foreground pixel that has not
            % already been labeled.
            if M(r,c) && L(r,c) == 0

                num = num + 1;

                head = 1;
                tail = 1;

                qR(tail) = r;
                qC(tail) = c;
                L(r,c) = num;

                % Breadth-first flood fill for this component.
                while head <= tail
                    cr = qR(head);
                    cc = qC(head);
                    head = head + 1;

                    for k = 1:size(nbrs,1)
                        rr = cr + nbrs(k,1);
                        col = cc + nbrs(k,2);

                        if rr >= 1 && rr <= m && col >= 1 && col <= n
                            if M(rr,col) && L(rr,col) == 0
                                tail = tail + 1;
                                qR(tail) = rr;
                                qC(tail) = col;
                                L(rr,col) = num;
                            end
                        end
                    end
                end
            end
        end
    end
end