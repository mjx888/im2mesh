function [lat, lon] = polymergeOS(lat, lon)
% polymergeOS: Merge line segments with matching endpoints
% Open-source version of function polymerge
%
% Assume The inputs LAT and LON can be column or row vectors with 
% NaN-separated parts
% Assume nargin < 3
% Assume tol = 0
%
%   [latMerged, lonMerged] = polymergeOS(LAT, LON) accepts a multipart line
%   in latitude-longitude with vertices stored in arrays LAT and LON, and
%   merges the parts wherever a pair of end points coincide. For this
%   purpose, an end point can be either the first or last vertex in a given
%   part. When a pair of parts are merged, they are combined into a single
%   part and the duplicate common vertex is removed.  If two first vertices
%   coincide or two last vertices coincide, then the vertex order of one of
%   the parts will be reversed.  A merge is applied anywhere that the end
%   points of exactly two distinct parts coincide, so that an indefinite
%   number of parts can be chained together in a single call to polymergeOS.
%   If three or more distinct parts share a common end point, however, the
%   choice of which parts to merge is ambiguous and therefore none of the
%   corresponding parts are connected at that common point.
%
%   The inputs LAT and LON can be column or row vectors with NaN-separated
%   parts (and identical NaN locations in each array) or they can be cell
%   arrays with each part in a separate cell. The form of the output
%   arrays, latMerged and lonMerged, matches the inputs in this regard.
%
%   See also POLYSPLIT, POLYJOIN.
% Copyright 1996-2012 The MathWorks, Inc.
%
% example:
%     lat = [3 2 NaN 1 2 NaN 5 6 NaN 3 4]';
%     lon = [13 12 NaN 11 12 NaN 15 16 NaN 13 14]';
%     [latm, lonm] = polymergeOS(lat,lon)

    % Assume The inputs LAT and LON can be column or row vectors with 
    % NaN-separated parts
    % Assume nargin < 3
    % Assume tol = 0
    
    % Parse and validate inputs.
    checklatlonOS(lat, lon, mfilename, 'LAT', 'LON', 1, 2)
    
    tol = 0;
    
    % Work with column vectors.
    rowVectorInput = (size(lat,2) > 1);
    lat = lat(:);
    lon = lon(:);
    
    % Ensure non-empty, NaN-terminated input arrays.
    [lat, lon] = removeExtraNanSeparatorsOS(lat, lon);
    if isempty(lat)
        lat = NaN;
        lon = NaN;
    elseif ~isnan(lat(end))
        lat(end+1,1) = NaN;
        lon(end+1,1) = NaN;
    end
    
    % Attempt to merge if the line has two or more parts.
    numberOfParts = sum(isnan(lat(:)));
    if numberOfParts > 1
        [lat, lon] = mergeMultipart(lat, lon, ...
            @(lat1, lon1, lat2, lon2) hypot(lat1 - lat2, lon1 - lon2), ...
            @midpoint, tol);
    end
    
    % Convert back to row vectors, if necessary.
    if rowVectorInput
        lat = lat';
        lon = lon';
    end


%--------------------------------------------------------------------------

function [x, y] = midpoint(x1, y1, x2, y2)

x = mean([x1 x2], 2);
y = mean([y1 y2], 2);

%--------------------------------------------------------------------------

function [c1Merged, c2Merged] ...
    = mergeMultipart(c1, c2, distanceFcn, midpointFcn, tol)
% Merge a multipart line. Inputs c1 and c2 are NaN-delimited (and
% terminated) column vectors.  Input distanceFcn is a handle to a function
% that computes an array of distances given two pairs of coordinate arrays.
% Input midPointFcn is a handle to a function that computes the point
% halfway between a pair of points. Input TOL is the maximum distance at
% which two points are considered to be coincident.

% The first step is to determine the indices in the vertex arrays c1 and c2
% that correspond to the starting points ("first") and ending points
% ("last") of each part. From there, we compute a "connection" vector with
% two elements per part.

[first, last, connection] = identifyConnections(c1, c2, distanceFcn, tol);

% Before tracing through the connections, bring pairs of vertices that are
% separated by less than TOL into perfect coincidence by replacing each
% member of the pair with their midpoint. 

if tol > 0
    % Note: The adjustment process here and the loop below, in which the
    % full set of vertex coordinates are copied to the output, both involve
    % tracing through the connection vector.  Beyond that similarity,
    % however, the two processes are very different.  Here we make sure to
    % adjust each end point that is part of a pair to be merged -- so we
    % "visit" all of the non-zero elements of connection.  Later we visit
    % each of the "free" end points, which are indicated by values of 0 in
    % the connection vector.  The paths through the connection vector will
    % be very different in these two cases.  This is one reason why it
    % works well to perform the adjustments via a pre-processing step.
    [c1, c2] = adjustCoincidentVertices( ...
        c1, c2, first, last, connection, midpointFcn);
end

% Excluding identities, loops (parts that start and end in the same place),
% and junctions at which more than two ending or starting points coincide,
% each element of the connection vector corresponds to either the starting
% point or ending point of one of the parts, arranged in the following
% order:
%
%                   [start1 end1 start2 end2 ...]
%
% If one part connects to another part (with a pair-wise coincidence of
% starting/ending vertices), then the element corresponding to the first
% part will contain the index of the element corresponding to the second
% part. On the other hand, if a starting or ending point does not coincide
% with the starting or ending point of another part, then value of its
% "connection" element will be 0.  Thus, the essential step in merging the
% input parts is to trace these self-references through the connection
% vector, starting at a zero-valued element, and finishing up and starting
% over with another new output part whenever another 0 is encountered.
%
% For example, if the value of the k-th element of the connections vector
% is 7 and k happens to be even, this means that the ending point of the
% p-th part (where p = k / 2) coincides with the starting point of the 4-th
% part (because 4 == (7 + 1) / 2).  Connections like this can be traced
% from one part to another.  Continuing with this example, we can move on
% to the ending point 4-th part by looking the 8-th element of the
% connection vector (because 8 == 7 + 1 = 2*4 - 1).  Suppose its value is
% 12; this means that it coincides with the end point of the 6-th part.
% Arriving at an ending point (as opposed to a starting point) simply means
% that we have to reverse the vertices during the merge.  So in this
% example we move forward as we copy the vertices from 4-th part to the
% output, but then copy the 6-th part in reverse order. To see where to go
% from starting point of the 6-th part, we look at the 11-th element (11 =
% 2*6 - 1 = 12 - 1) of the connection vector.  Suppose that there is no
% starting or ending vertex that coincides with this point; in this case,
% the 11-th element of connections will contain value 0.
%
% In some cases, one or more input part will stand apart and not need to be
% merged with any other part. If this is the case for the 10th part, for
% example, then elements 19 and 20 of the connection vector will contain 0.

% Initialize the output vertex arrays, filling them with NaN. The size of
% the inputs provides a safe upper bound on the size of the outputs,
% because there is no way for the number of elements to increase during the
% merge process: each part-to-part connection results in the loss of one
% NaN-separator and one set of vertex coordinates.
c1Merged = NaN(size(c1));
c2Merged = c1Merged;

% Iterate through the connection vector. Use the variable k to denote its
% k-th element.  Use p to denote the index of the corresponding part. 
k = 0;
nextOutputIndex = 1;
endpoints = find(connection == 0);
visitedEnd = false(size(endpoints));
visitedAll = false(size(connection));
while ~all(visitedEnd)
    if k == 0
        % Start a new output part.
        k = endpoints;
        k(visitedEnd) = [];
        k = k(1);
    end
    visitedEnd(endpoints == k) = true;
    visitedAll(connection == k) = true;

    % Update k and get indices needed for copy.    
    [k, n1, n2, m1, m2, step] ...
        = traversePartGettingIndices(k, first, last, nextOutputIndex);
    
    % Copy vertices.
    c1Merged(n1:n2) = c1(m1:step:m2);
    c2Merged(n1:n2) = c2(m1:step:m2);
    
    visitedAll(connection == k) = true;    
    if connection(k) == 0
        % End the current output part.
        visitedEnd(endpoints == k) = true;
        k = 0;
        
        % Add 2 in order to move past the vertex c1Merged(n2), c2Merged(n2)
        % and the NaN-separator that follows it.
        nextOutputIndex = n2 + 2;
    else
        % Append to the current output part.
        k = connection(k);
        
        % Allow c1Merged(n2), c2Merged(n2) to be overwritten by the
        % coincident vertex from the next input part.
        nextOutputIndex = n2;
    end
end

% At this stage, all parts with a free endpoint, or merged to a part
% with a free endpoint, have been visited.  But one or more closed loops
% (either single part or multipart) could remain.  They need to be
% traversed also.

k = 0;
visitedAll = visitedAll | (connection == 0);
while ~all(visitedAll)
    if k == 0
        k = connection;
        k(visitedAll) = [];
        k = k(1);
    end
    visitedAll(connection == k) = true;
    
    % Update k and get indices needed for copy.
    [k, n1, n2, m1, m2, step] ...
        = traversePartGettingIndices(k, first, last, nextOutputIndex);
    
    % Copy vertices.
    c1Merged(n1:n2) = c1(m1:step:m2);
    c2Merged(n1:n2) = c2(m1:step:m2);
    
    visitedAll(connection == k) = true;    
    if visitedAll(connection == connection(k))
        % End the current output part
        k = 0;
        
        % Add 2 in order to move past the vertex c1Merged(n2), c2Merged(n2)
        % and the NaN-separator that follows it.
        nextOutputIndex = n2 + 2;
    else
        % Append to the current output part.
        k = connection(k);
        
        % Allow c1Merged(n2), c2Merged(n2) to be overwritten by the
        % coincident vertex from the next input part.
        nextOutputIndex = n2;
    end
end

% After all the input vertices have been copied, it's very likely that
% there are extra NaN-separators that need to be removed from the end of
% the output vertex arrays.
[c1Merged, c2Merged] = removeExtraNanSeparatorsOS(c1Merged, c2Merged);

%--------------------------------------------------------------------------

function [k, n1, n2, m1, m2, step] ...
    = traversePartGettingIndices(k, first, last, nextOutputIndex)
% Traverse the current part, updating the endpoint index k, and getting the
% indices needed to copy the vertices from the input (n1:n2) to the output
% (m1:step:m2) vertex arrays.

% Having arrived at one end of a part as identified by the endpoint
% number k, identify the input part number p and then reset k to
% correspond to the opposite end of the part. Get the indices needed to
% copy the vertices from the p-th part to the output vertex arrays,
% taking care to move in the correct direction (from first to last if
% the k is initially odd, and from last to first if k is initially
% even).
if isodd(k)
    % Copy this part as-is (from first vertex to last vertex).
    k = k + 1;
    p = k / 2;
    m1 = first(p);
    m2 = last(p);
    step = 1;
else
    % k is even; reverse vertex direction when copying this part.
    p = k / 2;
    k = k - 1;
    m1 = last(p);
    m2 = first(p);
    step = -1;
end
n1 = nextOutputIndex;
n2 = nextOutputIndex + last(p) - first(p);

%--------------------------------------------------------------------------

function [first, last, connection] ...
    = identifyConnections(c1, c2, distanceFcn, tol)
% Determine the indices of the starting points ("first") and ending points
% ("last") in the vertex arrays c1 and c2.  Compute a "connection" vector
% with one element for each starting point and one element for each ending
% point, ordered as interleaved first-last pairs.  Thus connection(k), for
% even k, corresponds to the last vertex in the p-th part, where p = k/2.
% If this vertex coincides (within tolerance TOL) in a unique, pair-wise
% sense with the first or last vertex of another part, then the value of
% connection(k) is the index (in connection itself) of the other vertex in
% the pair. Otherwise it is zero.  For odd k, connection(k) corresponds to
% the first vertex in part p = (k + 1)/2.
%
% A "unique, pair-wise" coincidence is ensured by imposing three types of
% exclusions:
%
% 1. Identities (every endpoint coincides with itself, but these
%    coincidences are not relevant when merging distinct parts)
%
% 2. Loops: parts that start and end in the same place (there are not
%    relevant either, when merging distinct parts)
%
% 3. Junctions at which more than two ending or starting points coincide
%    (resulting in ambiguities which would require additional information
%    to resolve)

n = isnan(c1(:));
first = find([true; n(1:end-1)]);
last  = find([n(2:end);  false]);

% Store the coordinates of the starting and ending vertices, ordered
% pair-wise ([start1 end1 start2 end2 ...]) in vectors e1 and e2.
endpoints = [first'; last'];
endpoints = endpoints(:);
e1 = c1(endpoints);
e2 = c2(endpoints);

% Compute a "coincidence matrix", imposing the exclusions described above.
% The matrix C is 2N-by-2N, where N is the number of parts in the multipart
% line defined by c1 and c2.
C = coincidenceMatrix(e1, e2, distanceFcn, tol);

% C is symmetric and consists mostly of 0s, with at most one 1 per column
% and row.  We can transform it into a 2N-by-1 vector of self-referential
% indices without loss of information.
[rows, cols] = find(C);
connection = zeros(1,size(C,1));
connection(cols) = rows;

%--------------------------------------------------------------------------

function C = coincidenceMatrix(e1, e2, distanceFcn, tol)
% For pair-wise ordered endpoint coordinate vectors e1 and e2, each having
% length M, compute an M-by-M matrix C such that C is true for each pair of
% coincident points, excluding identities, loops,  and junctions of three
% or more parts.  (M == 2*N where N is the number of parts.)

M = numel(e1);

% Replicate columns.
e1 = e1(:,ones(1,M));
e2 = e2(:,ones(1,M));

% Compute raw coincidence matrix.
C = (distanceFcn(e1, e2, e1', e2') <= tol);

% Apply identity and loop exclusions: Set 2-by-2 submatrices along the
% diagonal to false. This includes the full diagonal itself and, for odd k,
% the pairs e1(k),e2(k) and e1(k+1),e2(k+1).
for k = 1:2:M
    C([k k+1],[k k+1]) = false;
end

% Apply multi-junction exclusion: In each row and column for which more
% than one element is true, set all elements to false.
sum1 = sum(C,1);
sum2 = sum(C,2);
C(:, sum1 > 1) = false;
C(sum2 > 1, :) = false;

%--------------------------------------------------------------------------

function [c1, c2] = adjustCoincidentVertices( ...
    c1, c2, first, last, connection, midpointFcn)
% Take the midpoint of coincident starting and ending vertex positions.

% Iterate through the connection vector, keeping track of status via the
% logical vector "visited".  Initialize visited to true for non-coincident
% endpoints, which are coded by the value 0. For each non-zero element of
% "connection" (with index k1), see if the element it points to (with index
% k2) is also non-zero (not an endpoint, that is). If so, then use the
% "first" and "last" arrays to determine the corresponding indices into the
% input vertex arrays c1 and c2 (i1 and i2, respectively), compute the
% midpoint location, and assign its coordinates back into c1 and c2.
visited = (connection == 0);
while ~all(visited)
    k1 = find(~visited);
    k1 = k1(1);
    visited(k1) = true;
    k2 = connection(k1);
    if k2 ~= 0
        visited(k2) = true;
        if isodd(k1)
            i1 = first((k1 + 1) / 2);
        else
            i1 = last(k1 / 2);
        end
        if isodd(k2)
            i2 = first((k2 + 1) / 2);
        else
            i2 = last(k2 / 2);
        end
        [c1Mid, c2Mid] = midpointFcn(c1(i1), c2(i1), c1(i2), c2(i2));
        c1(i1) = c1Mid;
        c2(i1) = c2Mid;
        c1(i2) = c1Mid;
        c2(i2) = c2Mid;
    end
end

%--------------------------------------------------------------------------

function tf = isodd(x)
% Return true if and only if x is an odd integer.

tf = (mod(x,2) == 1);
