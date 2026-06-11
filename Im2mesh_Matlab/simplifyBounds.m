function new_bounds = simplifyBounds( bounds, tolerance, thresh_num_vert )
% simplifyBounds: simplify polygonal boundaries using Douglas-Peucker 
% Polyline Simplification (dpsimplify.m)
%
% usage:
%   new_bounds = simplifyBounds( bounds, tolerance, thresh_num_vert );
%   new_bounds = simplifyBounds( bounds, tolerance );
%   new_bounds = simplifyBounds( bounds, 0 );   % no simplification
%
% input:
%   bounds - a nested cell array of 2d polygonal boundaries.
%            Polygons in bounds{i} belong to the i-th part or phase.
%            bounds{i}{j} is one of the polygons in the i-th part.
%            bounds{i}{j} is a n-by-2 array for x y coordinates of vertices
%            in a polygon. You can use 
%            plot( bounds{i}{j}(:,1), bounds{i}{j}(:,2) ) to view the
%            polygon. Use plotBounds( bounds ) to view all polygons.
%   bounds{i}{j} is a polygon boundary with (NaN,NaN). (NaN,NaN) is the 
%   label for control point.
% 
%   tolerance - parameter for polygon or polyline simplification.
%               Function: dpsimplify.m (Douglas–Peucker algorithm)
% 
%   thresh_num_vert - threshold for the number of vertices.
%                     If the number of vertices in a ppolyline is not 
%                     larger than this threshold, do not perform 
%                     simplification.
%                     It can be set as an integer or an array with 
%                     two elements. See section 4 in Tutorial.pdf
%
% output:
%   new_bounds - cell array. new_bounds{i}{j} is without (NaN,NaN)
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%
    
    % check the number of inputs
    if nargin == 2
        thresh_num_vert = 0;
    elseif nargin == 3
        % normal case
    else
        error("check the number of inputs");
    end
    
    new_bounds = bounds;
    
    % create a mapping vector for the actual number of tolerance
    tolMapVec = createMapVec( thresh_num_vert, tolerance );

    % simplify each polygonal boundary
    for i = 1: length(bounds)
        for j = 1: length(bounds{i})
            % bounds{i}{j} is a N-by-2 coordinate array of a polygon boundary with (NaN,NaN)
            
            % convert a polygon to cell array that consists of polylines
            % x, y are N-by-1 cell arrays with one polygon segment per cell
            [x, y] = polysplitOS( bounds{i}{j}(:,1), bounds{i}{j}(:,2) );
            
            % simplify each polyline in the polygon
            for k = 1: numel(x)
                % poly_O is the k-th polyline in the polygon
                poly_O = [ x{k}, y{k} ];  % N-by-2
                
                if poly_O(1,:) == poly_O(end,:)
                    num_vert = length(poly_O) - 1;
                else
                    num_vert = length(poly_O);
                end

                % the actual number of tolerance
                real_tol = getRealValue( num_vert, tolMapVec );
                
                % pre-process polyline
                % dpsimplify() is sensitive to the orientation of 
                % polyline, so reorient first
                poly_O = reorient( poly_O );

                % simplify polyline
                poly_temp = dpsimplify( poly_O, real_tol );
                
                % update
                x{k} = poly_temp(:,1);
                y{k} = poly_temp(:,2);
            end
            
            new_bounds{i}{j} = [];
            [ new_bounds{i}{j}(:,1), new_bounds{i}{j}(:,2) ] = polyjoinOS(x, y);
        end
    end
    
    new_bounds = mergeBounds( new_bounds );     % delete (NaN,NaN)
end

function polyline = reorient( polyline )
% reorient: make the starting point of polyline the left-bottom. 
% If polyline is a polygon, make it to counter clockwise.

    x = polyline(:,1);
    y = polyline(:,2);

    if x(1) == x(end) && y(1) == y(end)
        % the polyline is a polygon, make it to counter clockwise
        if ispolycwOS(x, y)
            x = x( end:-1:1 );  % convert to counter clockwise
            y = y( end:-1:1 );
        end
    else
        % the polyline is not a polygon
        % make the starting point of polyline the left-bottom
        if x(1) < x(end) || ...
                ( x(1) == x(end) && y(1) <= y(end) )
            % starting point of polyline already at left-bottom
            % so do nothing
        elseif ( x(1) == x(end) && y(1) > y(end) ) ||...
                x(1) > x(end)
            % starting point not at left-bottom
            % reverse
            x = x(end:-1:1);
            y = y(end:-1:1);
        else
            error('other cases')
        end
    end
    
    % update polyline
    polyline(:,1) = x;
    polyline(:,2) = y;
            
end

function bounds = mergeBounds( bounds )
% delete (NaN,NaN) in bounds{i}{j} using polymerge

    for i = 1: length(bounds)
        for j = 1: length(bounds{i})
            poly = [];
            [ poly(:,1), poly(:,2) ] = polymergeOS( ...
                                    bounds{i}{j}(:,1), bounds{i}{j}(:,2) );
            
            if isnan( poly(end,1) )
                poly = poly( 1:end-1, : );  % delete last NaN
            else
                error('not ending with NaN');
            end
            bounds{i}{j} = poly;
        end
    end
end

function real_value = getRealValue( num_vert, mapVec )
% getRealValue: get the real number (piecewise)
% 
% real_value is an integer. It's the real number when the
% number of vertices is num_vert.
%

    if num_vert > numel(mapVec)
        real_value = mapVec(end);
    else
        real_value = mapVec( num_vert );
    end
end

function mapVec = createMapVec( thresh_bound, ymax )
% createMapVec: create a mapping vector
    
    % ---------------------------------------------------------------------
    % chekc input
    t = thresh_bound;
    
    if numel(t) == 1
        t(2) = t(1);
    end

    if numel(t) == 2 && t(1) > t(2)
        % switch value to make t(1) < t(2)
        temp = t(2);
        t(2) = t(1);
        t(1) = temp;
    end

    if numel(t) > 2
        error('Num of elements in thresh_num_vert should not larger than 2.');
    end
    
    % ---------------------------------------------------------------------
    % initialize
    if t(2) > 0
        len = ceil( t(2)+1 );
    else
        len = 2;
    end

    mapVec = ymax * ones(len,1);
    
    % ---------------------------------------------------------------------
    % single threshold
    if numel(t) == 1 || t(1) == t(2)
        mapVec( 1: ceil(t(1)-1) ) = 0;
        return
    end

    % ---------------------------------------------------------------------
    % two thresholds
    % t = round(t);     % !!!
    
    if t(2) <= 0
        return
    end
    
    if t(1) > 0
        mapVec( 1: t(1) ) = 0;
        k = ymax/(t(2)-t(1));
        x = t(1) : t(2);
        x = x(:);
        y = k * ( x - t(1));
        mapVec( x ) = y;    % !!!
        
    elseif t(1) <= 0 && t(2) > 0
        k = ymax/(t(2)-t(1));
        x = t(1) : t(2);
        x = x(:);

        idx = x<=0 ;
        x(idx) = [];

        y = k * ( x - t(1));
        mapVec( x ) = y;    % !!!
    end

    % ---------------------------------------------------------------------
end