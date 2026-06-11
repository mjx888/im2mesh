function new_bounds = smoothBounds( bounds, lambda, mu, iters, thresh_num_turn, thresh_num_vert )
% smoothBounds: smooth polygonal boundaries using 2d Taubin Smoothing 
% (taubinSmooth.m)
%
% usage:
%   new_bounds = smoothBounds( bounds, lambda, mu, iters, thresh_num_turn, thresh_num_vert );
%   new_bounds = smoothBounds( bounds, lambda, mu, iters, thresh_num_turn );
%   new_bounds = smoothBounds( bounds, lambda, mu, iters );
%   new_bounds = smoothBounds( bounds, lambda, mu, 0 );     % no smoothing
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
%   lambda, mu, iters - parameters for Taubin Smoothing (taubinSmooth.m)
%
%   lambda: How far each node is moved toward the average position of its
%           neighbours during every second iteration. 0<lambda<1
%
%   mu: How far each node is moved opposite the direction of the average 
%       position of its neighbours during every second iteration. -1<mu<0
%
%   iters: Number of smoothing steps. Non negative interger
%          If iters == 0, no smooth.
%
%   thresh_num_turn:  Threshlod for the number of knick-point 
%                     vertices or turning points. 
%                     If the number of knick-point vertices or 
%                     turning points in a polyline is smaller than 
%                     this threshold, do not perform smoothing. 
%
%   thresh_num_vert: Threshlod for the number of vertices. 
%                    If the number of vertices in a polyline is smaller
%                    than this threshold, do not perform smoothing.
%                    It can be set as an integer or an array with 
%                    two elements. See section 4 in Tutorial.pdf
%
% output:
%   new_bounds - cell array. new_bounds{i}{j} is a polygon boundary with 
%               (NaN,NaN). (NaN,NaN) is the label for control point.
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%
    
    % check the number of inputs
    if nargin == 4
        thresh_num_turn = 0;
        thresh_num_vert = 0;
    elseif nargin == 5
        thresh_num_vert = 0;
    elseif nargin == 6
        % normal case
    else
        error('check the number of inputs');
    end

    new_bounds = bounds;

    if iters == 0
        % no smoothing
        return
    end

    % create a mapping vector for the actual number of iterations
    iterMapVec = createMapVec( thresh_num_vert, iters );
    
    % smooth each polygonal boundary
    for i = 1: length(bounds)
        for j = 1: length(bounds{i})
            % bounds{i}{j} is a N-by-2 coordinate array of a polygon boundary with (NaN,NaN)

            % convert a polygon to cell array that consists of polylines
            % x, y are N-by-1 cell arrays with one polygon segment per cell
            [x, y] = polysplitOS( bounds{i}{j}(:,1), bounds{i}{j}(:,2) );
            
            % smooth polyline in the polygon
            for k = 1: numel(x)
                % poly_O is the k-th polyline in the polygon
                poly_O = [ x{k}, y{k} ];  % N-by-2
                
                % use function dpsimplify() to remove non-knick/non-turning points
                poly_turn_pnt = dpsimplify( poly_O, eps );

                % check the number of turning points in a polyline
                if length(poly_turn_pnt)-1 < thresh_num_turn
                    % If the number of turning-point vertices in a polyline
                    % is smaller than thresh_num_turn, don't smooth.
                    continue
                else
                    if poly_O(1,:) == poly_O(end,:)
                        num_vert = length(poly_O) - 1;
                    else
                        num_vert = length(poly_O);
                    end

                    % the actual number of iterations
                    real_iter = getRealValue( num_vert, iterMapVec );
                    
                    % smooth polyline using taubinSmooth function
                    poly_temp = taubinSmooth( poly_O, lambda, mu, real_iter );
                    
                    % update
                    x{k} = poly_temp(:,1);
                    y{k} = poly_temp(:,2);
                end
            end
            
            [ new_bounds{i}{j}(:,1), new_bounds{i}{j}(:,2) ] = polyjoinOS(x, y);
        end
    end

    % rounding the vertex coordinates to n digits
    % purpose: avoid numeric error of taubinSmooth
    n_digit_decimal = 2;    % N digits to the right of the decimal point
    new_bounds = roundBounds( new_bounds, n_digit_decimal );

end


function bounds = roundBounds( bounds, n_digit_decimal )
% round the vertex coordinates in polygonal boundaries to n digits
% input:
%   n_digit_decimal: N digits to the right of the decimal point

    btemp = bounds;

    for i = 1: length(btemp)
        for j = 1: length(btemp{i})
             btemp{i}{j} = round( btemp{i}{j}, n_digit_decimal );
        end
    end
    
    bounds = btemp;

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
    t = round(t);   % !!!
    
    if t(2) <= 0
        return
    end
    
    if t(1) > 0
        mapVec( 1: t(1) ) = 0;
        k = ymax/(t(2)-t(1));
        x = t(1) : t(2);
        x = x(:);
        y = k * ( x - t(1));
        mapVec( x ) = round( y );
        
    elseif t(1) <= 0 && t(2) > 0
        k = ymax/(t(2)-t(1));
        x = t(1) : t(2);
        x = x(:);

        idx = x<=0 ;
        x(idx) = [];

        y = k * ( x - t(1));
        mapVec( x ) = round( y );
    end

    % ---------------------------------------------------------------------
end









