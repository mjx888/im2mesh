function Bs = getExactBoundsOS( bw )
% getExactBoundsOS: get the exact boundaries (polygonal) of binary image
% sub-function: outerboundary.m and holeboundary.m
% getExactBoundsOS is the fully open-source version of getExactBounds
%
% Use bwlabelOS.m, regionpropsOS.m, bwboundariesOS.m
% 
% usage:
%	Bs = getExactBoundsOS( bw );
%
% input bw - binary image
% ouput Bs - stores column (x) and row (y) coordinates of exact polygonal 
%           boundaries. Bs is a p-by-1 cell array, where p is the number of
%           (4 connected) objects and (4 connected) holes. Each cell in the
%           cell array contains a q-by-2 matrix. Each row in the matrix 
%           contains the column and row coordinates of a boundary corner. q
%           is the number of boundary corner for the corresponding region.
%           Size of Bs{i} is (1+number_of_vertices)-by-2.
%
% Comment:
% There're huge different between function bwboundaries and getExactBounds.
% The output of bwboundaries and getExactBounds may have different length.
%
% 1. bwboundaries(BW,conn,options)
% Return row (y) and column (x) coordinates of boundary pixel locations.
% Only inner boundaries are obatined.
% Parameter conn is very important for function bwboundaries. According to 
% bwboundaries.m in Matlab, 
% if conn = 4, objects are 4 connected and holes are 8 connected; 
% if conn = 8, objects are 8 connected and holes are 4 connected.
% This is to avoid topological errors.
%
% 2. getExactBounds( bw )
% Return column (x) and row (y) coordinates of the corner of boundary
% pixels.
% However, function getExactBounds don't have connectivity parameter. Both 
% objects and holes are 4 connected, since outputs are exact boundary. 
%
% Example:
%     bw = imread('bw.tif');
%     Bs = getExactBounds( bw );
% 
%     hold on
%     imshow(bw,'InitialMagnification','fit');
%     for k = 1:length(Bs)
%        plot( Bs{k}(:,1), Bs{k}(:,2), 's-g','MarkerSize',5,'LineWidth',2 );
%     end
%     hold off
% 
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%   

    % initialize
    Bs = cell(0);
    % label objects using 4-connected neighborhood
    [ Label, num_object ] = bwlabelOS( bw, 4 );
    objs = regionpropsOS( Label, 'Image', 'BoundingBox' );

    % get boundary of each 4-connected object, and update Bs
    for i = 1: num_object
        % one object at a time
        % B_temp is a cell & inner boundary
        B_temp = bwboundariesOS( (objs(i).Image)', 8 );

        % convert back to global coordinate
        for j = 1: length(B_temp)
            B_temp{j}(:,1) = B_temp{j}(:,1) + objs(i).BoundingBox(1) -0.5;  
            B_temp{j}(:,2) = B_temp{j}(:,2) + objs(i).BoundingBox(2) -0.5;
        end

        % convert B_temp (inner boundary) to B (exact boundary)
        if length(B_temp) == 1       % no holes
            B = outerboundary( B_temp(1) );
        elseif length(B_temp) > 1    % exist holes
            B = [ outerboundary( B_temp(1) ); 
                holeboundary( B_temp(2:end) ) ];
        else
            error("boundary error")
        end

        % update Bs
        Bs = [ Bs; B ];
    end

end


function cellout = outerboundary( cellin )
% outerboundary: get exact boundary of outer contour of objects
% sub-function of getExactBounds.m
% A faster version of the original outerboundary function.
%
% INPUT:
%   cellin : N-by-1 cell. Each cell is a polygon, M-by-2, in clockwise order.
%
% OUTPUT:
%   cellout: N-by-1 cell. Each cell is the "outer boundary" polygon for the
%            corresponding input polygon.
%
% an example:
%     bw = [0 0 0 0 0;
%           0 0 1 1 0;
%           0 1 1 1 0;
%           0 0 0 1 0;
%           0 0 0 0 0];
%     bw = logical( bw );
%     B_in = bwboundaries( bw', 8, 'noholes' );
%     B_ex = outerboundary( B_in );
% 
%     hold on
%     imshow(bw,'InitialMagnification','fit');
%     for i = 1: length(B_ex)
%         plot( B_in{i}(:,1), B_in{i}(:,2), 'O-k','MarkerSize',5,'LineWidth',2 );
%         plot( B_ex{i}(:,1), B_ex{i}(:,2), 'v-g','MarkerSize',5,'LineWidth',2 );
%     end
%     hold off
%     legend( 'result of bwboundaries(bw'',8,''noholes'' )' ,...
%             'outerboundary( bwboundaries() )');
%
%
% Feb 2025 update note:
% Preallocation of arrays so that we do not grow or shrink arrays repeatedly inside loops.
% Single-pass insertions rather than [A(1:k), newValues, A(k+1:end)] in a loop. 
% Repeated concatenation or slicing is notoriously costly in MATLAB. 
% Instead, we collect any newly needed points in one pass, then add them all at once.
% Minor vectorization and reduction in repeated checks where possible.
% Removed a few redundant operations in places.
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Feb 2025

nPolys = numel(cellin);
cellout = cell(nPolys, 1);

for j = 1 : nPolys
    currBoundary = cellin{j};
    nPoints      = size(currBoundary, 1);
    
    % Quick sanity check
    if nPoints < 2
        error('size error')
    elseif nPoints == 2
        % Single pixel boundary
        [cx, cy]        = meshgrid( currBoundary(1,1)+[-.5 .5], ...
                                    currBoundary(1,2)+[-.5 .5] );
        pnt_list        = zeros(5,2);
        pnt_list(1:4,:) = [cx(:), cy(:)];
        pnt_list(3:4,:) = flip( pnt_list(3:4,:), 1 );
        pnt_list(5,:)   = pnt_list(1,:);
        new_pnt_list    = pnt_list;
    else
        % More than one pixel
        % -----------------------------------------------------------
        % STEP 1 & STEP 2
        % -----------------------------------------------------------
        % We'll generate a set of "midpoints" (offset by +/- 0.5) and
        % also add special corner vertices. We do so in a single pass.
        
        % Add an extra row to close it easily
        poly = [currBoundary; currBoundary(2,:)];

        % The maximum number of points we can generate in step 1 & 2:
        %   - We get 1 new "midpoint" for each edge (nPoints-1 edges)
        %   - We might add up to 2 new corner points for each edge
        % So let's be generous and preallocate for 3*(nPoints-1) + 1 
        % (the +1 for closure).
        
        maxSz   = 3*(nPoints-1) + 1;
        pnt_list = zeros(maxSz, 2);
        count    = 0;

        for k = 1 : (size(poly,1) - 1)
            p1 = poly(k,:);
            p2 = poly(k+1,:);

            % STEP 1: midpoints offset by +/- 0.5
            if p1(1) == p2(1)
                % vertical line
                tempY = 0.5*(p1(2)+p2(2));
                if p2(2) > p1(2)
                    tempX = p1(1) - 0.5;
                else
                    tempX = p1(1) + 0.5;
                end
                midpt = [tempX, tempY];

            elseif p1(2) == p2(2)
                % horizontal line
                tempX = 0.5*(p1(1)+p2(1));
                if p2(1) > p1(1)
                    tempY = p1(2) + 0.5;
                else
                    tempY = p1(2) - 0.5;
                end
                midpt = [tempX, tempY];

            else
                % skew line => just store midpoint for now
                midpt = 0.5*(p1 + p2);
            end

            % STEP 2: handle "special corner" additions if the newly
            % created point + previously created point form a
            % horizontal or vertical line, we might add corners.
            if count > 0
                pLast = pnt_list(count,:);
                pCurr = midpt;
                
                % vertical check
                if abs(pLast(1)-pCurr(1)) < 1e-12
                    % vertical line from pLast to midpt
                    % We check poly(k,:) to see if we need corners
                    if pCurr(2) > pLast(2)
                        if p1(1) < pLast(1)
                            corner = [ p1 + [-0.5 -0.5];
                                       p1 + [-0.5 +0.5] ];
                        else
                            corner = [];
                        end
                    else
                        if p1(1) > pLast(1)
                            corner = [ p1 + [+0.5 +0.5];
                                       p1 + [+0.5 -0.5] ];
                        else
                            corner = [];
                        end
                    end

                    nC = size(corner,1);
                    if nC>0
                        pnt_list(count+1:count+nC,:) = corner;
                        count = count + nC;
                    end

                % horizontal check
                elseif abs(pLast(2)-pCurr(2)) < 1e-12
                    % horizontal line from pLast to midpt
                    if pCurr(1) > pLast(1)
                        if p1(2) > pLast(2)
                            corner = [ p1 + [-0.5 +0.5];
                                       p1 + [+0.5 +0.5] ];
                        else
                            corner = [];
                        end
                    else
                        if p1(2) < pLast(2)
                            corner = [ p1 + [+0.5 -0.5];
                                       p1 + [-0.5 -0.5] ];
                        else
                            corner = [];
                        end
                    end

                    nC = size(corner,1);
                    if nC>0
                        pnt_list(count+1:count+nC,:) = corner;
                        count = count + nC;
                    end
                end
            end

            % Put in the midpoint
            count = count + 1;
            pnt_list(count,:) = midpt;
        end

        % Trim unused space
        pnt_list = pnt_list(1:count,:);
        
        % Check closure
        if ~isequal(pnt_list(1,:), pnt_list(end,:))
            error('polygon not closed after step 1/2.')
        end

        % -----------------------------------------------------------
        % STEP 3: For every pair of consecutive points in pnt_list,
        % if they form a skew line (both x and y differ), we insert
        % one extra corner point. We'll do a single pass and store
        % new points in a second list, then combine them.
        % -----------------------------------------------------------
        nSegs         = size(pnt_list,1)-1;  % number of edges
        cornerFlags   = false(nSegs,1);
        cornerPoints  = zeros(nSegs,2);

        for k = 1 : nSegs
            p1 = pnt_list(k,:);
            p2 = pnt_list(k+1,:);

            % If both x and y differ => skew
            if (abs(p1(1)-p2(1)) > 1e-12) && (abs(p1(2)-p2(2)) > 1e-12)
                cornerFlags(k) = true;
                % figure out which corner
                if (p2(1) < p1(1)) && (p2(2) > p1(2))
                    cornerPoints(k,:) = [ p2(1), p1(2) ];
                elseif (p2(1) > p1(1)) && (p2(2) > p1(2))
                    cornerPoints(k,:) = [ p1(1), p2(2) ];
                elseif (p2(1) > p1(1)) && (p2(2) < p1(2))
                    cornerPoints(k,:) = [ p2(1), p1(2) ];
                elseif (p2(1) < p1(1)) && (p2(2) < p1(2))
                    cornerPoints(k,:) = [ p1(1), p2(2) ];
                end
            end
        end

        % Now build the final new_pnt_list with inserted corners
        % We'll have at most nSegs extra points
        new_pnt_size = size(pnt_list,1) + sum(cornerFlags);
        new_pnt_list = zeros(new_pnt_size,2);

        idxNew = 1;
        for k = 1 : nSegs
            new_pnt_list(idxNew,:) = pnt_list(k,:);
            idxNew = idxNew + 1;
            if cornerFlags(k)
                new_pnt_list(idxNew,:) = cornerPoints(k,:);
                idxNew = idxNew + 1;
            end
        end
        % close the polygon
        new_pnt_list(idxNew,:) = pnt_list(end,:);

        % final sanity check
        if ~isequal(new_pnt_list(1,:), new_pnt_list(end,:))
            error('polygon not closed after step 3.')
        end
    end

    cellout{j} = new_pnt_list;
end

end

function cellout = holeboundary( cellin )
% holeboundary: get exact boundary of the holes
% sub-function of getExactBounds.m
% input - cellin is a cell (N*1), cellin{j} is clockwise polygonal boundary
%
% an example:
%     bw = [0 1 0 0 0;
%           1 0 1 1 0;
%           1 0 0 1 0;
%           1 0 0 1 0;
%           1 1 1 0 0];
%     bw = logical( bw );
%     B_in = bwboundaries( bw', 8 );    
%     % B_in{1} is boundary of the object
%     % B_in{2} is boundary of hole inside the object
%     B_ex = holeboundary( B_in(2) );
%     % B_ex is the exact boundary of hole inside the object
% 
%     hold on
%     imshow(bw,'InitialMagnification','fit');
%     plot( B_in{2}(:,1), B_in{2}(:,2), 'O-r','MarkerSize',5,'LineWidth',2 );
%     plot( B_ex{1}(:,1), B_ex{1}(:,2), 'v-g','MarkerSize',5,'LineWidth',2 );
%     hold off
%     legend( 'result of bwboundaries( )' ,...
%             'holeboundary( bwboundaries() )');
%
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, May 2020 

cellout = cell( length(cellin), 1 );

for j = 1: length(cellin)
    if size( cellin{j}, 1 ) < 2
        error('size error')
    elseif size( cellin{j}, 1 ) == 2
        % one pixel
        % four corner
        [cx,cy] = meshgrid( cellin{j}(1,1)+[-.5 .5], cellin{j}(1,2)+[-.5 .5] );
        pnt_list = zeros( 5, 2 );
        pnt_list( 1:4, : ) = [ cx(:), cy(:) ];
        pnt_list( 3:4, : ) = flip( pnt_list(3:4,:), 1 );
        pnt_list( 5, : ) = pnt_list( 1, : );
    else
        % two or more pixels
        poly = [ cellin{j}; cellin{j}(2,:) ];   % to make pnt_list be close
        len_pnt_list = ( size( cellin{j}, 1 ) - 1 )*4;
        pnt_list = zeros( len_pnt_list, 2 );    % initailize
        count = 0;
        
        % get exact boundary from inner boundary
        for k = 1: length(poly)-2
            % vector
            vec_1 = [ poly(k+1,1)-poly(k,1), poly(k+1,2)-poly(k,2) ];
            vec_2 = [ poly(k+2,1)-poly(k+1,1), poly(k+2,2)-poly(k+1,2) ];
            xy = [ poly( k, : ); poly( k+1, : ); poly( k+2, : ) ];
            % generate points on the corner of a pixel
            temp = getCorner( vec_1, vec_2, xy );
            
            % delete repeated vertex
            % compare last of pnt_list and first of temp
            if count>0 && isequal( pnt_list(count,:), temp(1,:) )
                temp(1,:) = [];
            end
            % put temp into pnt_list
            if ~isempty( temp )
                len_temp = size( temp, 1 );
                pnt_list( count+(1:len_temp), : ) = temp;
                count = count + len_temp;
            end
        end
        % delete redudant zeros in pnt_list
        pnt_list( count+1:end, : ) = [];
    end
    
    if ~isequal( pnt_list(1,:), pnt_list(end,:) )
        error('polygon not close')
    end
    
    cellout{j} = pnt_list;
end

end

function temp = getCorner( vec_1, vec_2, xy )
% generate points on the corner of a pixel
%
%     % vector
%     vec_1 = [ poly(k+1,1)-poly(k,1), poly(k+1,2)-poly(k,2) ];
%     vec_2 = [ poly(k+2,1)-poly(k+1,1), poly(k+2,2)-poly(k+1,2) ];
%     xy = [ poly( k, : ); poly( k+1, : ); poly( k+2, : ) ];
%     temp = getCorner( vec_1, vec_2, xy );
%

    % rotate xy, make vec_1 to [1 0] direction
    angl = angle( vec_1(1) + vec_1(2)*1i );
    Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];
    rot = round( Rmat( -angl ) );
    xyr = xy * rot;     % size = (n,2) * (2,2) = (n,2)

    % generate points on the corner of a pixel
    if isequal( vec_1, vec_2 )                       % case 1
        tempr = [ ( xyr(1,1) + xyr(2,1) )/2, xyr(1,2) + 0.5;
                  ( xyr(2,1) + xyr(3,1) )/2, xyr(1,2) + 0.5 ];

    elseif isequal( vec_1, -vec_2 )                  % case 2
        tempr = [ ( xyr(1,1) + xyr(2,1) )/2, xyr(1,2) + 0.5;
                  ( xyr(1,1) + xyr(2,1) )/2 + 1, xyr(1,2) + 0.5;
                  ( xyr(1,1) + xyr(2,1) )/2 + 1, xyr(1,2) - 0.5;
                  ( xyr(1,1) + xyr(2,1) )/2, xyr(1,2) - 0.5 ];

    elseif dot( vec_1, vec_2 ) == 0 && ...            % case 3
            vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1) > 0 % cross product
        tempr = [ ( xyr(1,1) + xyr(2,1) )/2, xyr(1,2) + 0.5 ];

    elseif dot( vec_1, vec_2 ) == 0 && ...            % case4
            vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1) < 0
        tempr = [ ( xyr(1,1) + xyr(2,1) )/2, xyr(1,2) + 0.5;
                  ( xyr(1,1) + xyr(2,1) )/2 + 1, xyr(1,2) + 0.5;
                  ( xyr(1,1) + xyr(2,1) )/2 + 1, xyr(1,2) - 0.5 ];
    else
        error('unexpected cases');
    end
    
    % rotate backward
    temp = tempr * rot';
end