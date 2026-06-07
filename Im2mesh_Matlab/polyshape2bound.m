function bounds = polyshape2bound( pC )
% polyshape2bound: Convert a cell array of polyshape objects to a cell 
% array of polygonal boundaries
%
% usage:
%	bounds = polyshape2bound( pC );
%
% input:
% pC - cell array of polyshape objects. pC{i} is a polyshape object, which 
%      corresponds to polygons in bounds{i}
%
% output:
%   bounds - a nested cell array of 2d polygonal boundaries.
%            Polygons in bounds{i} belong to the i-th part or phase.
%            bounds{i}{j} is one of the polygons in the i-th part.
%            bounds{i}{j} is a n-by-2 array for x y coordinates of vertices
%            in a polygon. You can use 
%            plot( bounds{i}{j}(:,1), bounds{i}{j}(:,2) ) to view the
%            polygon. Use plotBounds( bounds ) to view all polygons.
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % delete empty element in cell array
    pC = pC( ~cellfun('isempty', pC) );

    bounds = cell( length(pC), 1 );
    
    for i = 1: length(bounds)
        [ x, y ] = boundary(pC{i});
        [ xC, yC ] = polysplitOS( x, y );

        num_face = numel(xC);
        bounds{i} = cell( num_face, 1 );

        for j = 1: num_face
            bounds{i}{j} = [ xC{j}, yC{j} ];
        end
    end

    % post-process
    bounds = bounds2CCW( bounds );       % counterclockwise ordering
end


function bounds = bounds2CCW( bounds )
% Convert polygon contour to counterclockwise vertex ordering
% 
    for i = 1: length(bounds)
        for j = 1: length(bounds{i})
            [ bounds{i}{j}(:,1), bounds{i}{j}(:,2) ] = poly2ccwOS( ...
                                    bounds{i}{j}(:,1), bounds{i}{j}(:,2) );
            
        end
    end
end










