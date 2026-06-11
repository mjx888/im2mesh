function bounds = delZeroAreaPoly( bounds )
% delZeroAreaPoly: delete polygon with zero area
%
% usage:
%	bounds = delZeroAreaPoly( bounds );
%
%   bounds - a nested cell array of 2d polygonal boundaries.
%            Polygons in bounds{i} belong to the i-th part or phase.
%            bounds{i}{j} is one of the polygons in the i-th part.
%            bounds{i}{j} is a n-by-2 array for x y coordinates of vertices
%            in a polygon. You can use 
%            plot( bounds{i}{j}(:,1), bounds{i}{j}(:,2) ) to view the
%            polygon. Use plotBounds( bounds ) to view all polygons.
%
%   Jiexian Ma, mjx0799@gmail.com, Feb 2025

    %----------------------------------------------------------------------
    % find sub-polygon in bounds{i}{j} with zero area
    
    for i = 1: length(bounds)
        for j = 1: length(bounds{i})
            % check bounds{i}{j}

            % count NaN
            if sum(isnan( bounds{i}{j}(:,1) )) == 0
                % normal case, no NaN
                continue   
            end

            % exist NaN
            % something wierd happen
            % need to delete wierd sub-polygon with zero area

            % convert a polygon to cell array that consists of polylines
            % x, y are N-by-1 cell arrays with one polygon segment per cell
            [x, y] = polysplitOS( bounds{i}{j}(:,1), bounds{i}{j}(:,2) );
            
            mark_empty_segment = false( numel(x), 1 );

            for k = 1: numel(x)
                if polyarea( x{k}, y{k} ) == 0
                    mark_empty_segment(k) = true;
                end
            end

            % delete wierd sub-polygon with zero area
            x(mark_empty_segment) = [];
            y(mark_empty_segment) = [];

            if numel(x) ~= 1
                error("BAD input")
            end

            % update bounds{i}{j}
            bounds{i}{j} = [ x{1}, y{1} ];
        end
    end

    %---------------------------------------------------------------------
    % find polygon bounds{i}{j} with zero area

    mark_empty_bounds = false( length(bounds), 1 );
    
    for i = 1: length(bounds)
        mark_zero_poly = false( length(bounds{i}), 1 );
        % check bounds{i}{j}
        for j = 1: length(bounds{i})
            if polyarea( bounds{i}{j}(:,1), bounds{i}{j}(:,2) ) == 0
                mark_zero_poly(j) = true;
            end
        end      
        bounds{i}( mark_zero_poly ) = [];   % delete bounds{i}{j}
        
        % check bounds{i}
        if isempty(bounds{i})
            mark_empty_bounds(i) = true;
        end
    end
    
    bounds( mark_empty_bounds ) = [];       % delete bounds{i}
    %----------------------------------------------------------------------
end