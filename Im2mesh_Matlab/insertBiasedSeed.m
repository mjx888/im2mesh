function xyNew = insertBiasedSeed( xy, iters, ratio )
% insertBiasedSeed: inserts biased points between vertices of an 2d edge
%
% Take the coordinates x, y of a 2D edge, and returns new coordinates
% xNew, yNew. The new coordinates will have biased points inserted on each
% edge
%
% Biased points are inserted using the following equation:
%         xNew = ratio * x(1) + (1-ratio) * x(2);
%         yNew = ratio * y(1) + (1-ratio) * y(2);
%
% usage:
%	xyNew = insertBiasedSeed( xy, iters, ratio );
%
% input:
%	xy is a 2-by-2 array. Each row is a vertex for a 2D edge
%
%	iters is the number of iterations
%
%	ratio is a numeric, -1 < ratio < 1. Biased ratio.
%         if ratio is positive, normal case
%         if ratio is negative, flip biased direction
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % ---------------------------------------------------------------------
    % check the number of inputs
    if nargin == 1
        iters = 1;
        ratio = 0.5;
    elseif nargin == 2
        ratio = 0.5;
    elseif nargin == 3
        % all input arguments are given
    else
        error("Check the number of inputs");
    end

    if size(xy,1) ~= 2
        error('Input should be a 2-by-2 array. Each row is a vertex.')
    end
    
    % ---------------------------------------------------------------------
    % if ratio is negative, flip rows of xy
    if ratio < 0
        xy = flip(xy);
        ratio = -ratio;
    end

    % ---------------------------------------------------------------------
    % insert

    for i = 1: iters
        x = xy(:,1);
        y = xy(:,2);
    
        % number of vertex
        n = length(x);
    
        % Number of edges is (n-1). 
        % We add two vertices per edge + final closure.
        m = n + 1;
        xNew = zeros(m, 1);
        yNew = zeros(m, 1);
    
        % Fill new arrays
        xNew(1) = x(1);
        yNew(1) = y(1);
        
        xNew(2) = ratio * x(1) + (1-ratio) * x(2);
        yNew(2) = ratio * y(1) + (1-ratio) * y(2);
        
        % last vertex
        xNew(3:end) = x(2:end);
        yNew(3:end) = y(2:end);
    
        xyNew = [xNew, yNew];

        % update xy for next iteration
        xy = xyNew;
    end
    % ---------------------------------------------------------------------
end










