function [x2, y2] = poly2ccwOS(x1, y1)
% poly2ccwOS: Convert polygon contour to counterclockwise vertex ordering
% Open-source version of function poly2ccw
%
% Assume x1, y1 is vectors (not consist of NaN)
%
%   [X2, Y2] = poly2ccwOS(X1, Y1) arranges the vertices in the polygonal
%   contour (X1, Y1) in counterclockwise order, returning the result in X2
%   and Y2. If X1 and Y1 can contain multiple contours, represented either
%   as NaN-separated vectors or as cell arrays, then each contour is
%   converted to clockwise ordering.  X2 and Y2 have the same format
%   (NaN-separated vectors or cell arrays) as X1 and Y1.
%
%   Example
%   -------
%   Convert a clockwise-ordered square to counterclockwise ordering.
%
%       x1 = [0 0 1 1 0];
%       y1 = [0 1 1 0 0];
%       ispolycwOS(x1, y1)
%       [x2, y2] = poly2ccwOS(x1, y1);
%       ispolycwOS(x2, y2)
%
%   See also ISPOLYCW, POLY2CW, POLYSHAPE
% Copyright 2004-2017 The MathWorks, Inc.

    % Assume x1, y1 is vectors (not consist of NaN)

    input_is_cell = iscell(x1);
    
    if ~input_is_cell
       input_is_row = (size(x1, 1) == 1);
       [x1, y1] = polysplitOS(x1, y1);
    end
    
    x2 = x1;
    y2 = y1;
    
    for k = 1:numel(x1)
       if ispolycwOS(x2{k}, y2{k})
          x2{k} = x2{k}(end:-1:1);
          y2{k} = y2{k}(end:-1:1);
       end
    end
    
    if ~input_is_cell
       [x2, y2] = polyjoinOS(x2, y2);
       if input_is_row
          x2 = x2';
          y2 = y2';
       end
    end
