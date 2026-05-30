function poly = taubinSmooth(poly, lambda, mu, iters)
% taubinSmooth: Taubin smoothing (2D) without shrinkage.
% detail: https://pypi.org/project/shapelysmooth/#taubin
% 		  https://www.cse.wustl.edu/~taoju/cse554/lectures/lect06_Fairing.pdf
%
% usage:
%	poly = taubinSmooth(poly, lambda, mu, iters)
%
% inputs:
%    poly: N-by-2 array of (x,y)
%    lambda    : forward Laplacian factor   (0 < lambda < 1)
%    mu        : reverse Laplacian factor   (-1 < mu < 0)
%    iters     : number of iterations

    % Check if polyline is closed
    isClosed = isequal(poly(1,:), poly(end,:));
    
    % Repeat the pair of (lambda, mu) steps 'iters' times
    for k = 1: iters
        % Forward Laplacian
        poly = smoothStep(poly, lambda, isClosed);
        
        % Reverse Laplacian
        poly = smoothStep(poly, mu, isClosed);
    end
end

%--------------------------------------------------------------------------
function newLS = smoothStep(oldLS, factor, isClosed)
% Performs one Laplacian-like smoothing pass on oldLS with the given factor.
% For closed curves, replicates the endpoint-handling of the original code.

    n = size(oldLS, 1);
    newLS = oldLS; % Preallocate the updated curve

    if n < 3
        % Nothing to do for fewer than 3 points
        return;
    end
    
    % =========== Interior Points (Same for Open/Closed) ===========
    % Indices 2..(n-1)
    i = 2 : (n - 1);
    % Average neighbors for each interior point
    neighborAvg = 0.5 * (oldLS(i - 1, :) + oldLS(i + 1, :));
    % Laplacian step
    newLS(i, :) = oldLS(i, :) + factor .* (neighborAvg - oldLS(i, :));

    % =========== Endpoint Handling for Closed Polylines ===========
    if isClosed
        % The original code updates the first point separately using
        % the average of the second and second-last points, then copies
        % it to the last point. This is *not* a standard wrap-around
        % approach, so we replicate that exactly.

        % 2nd point:  oldLS(2, :),    second-last point: oldLS(end-1, :)
        endPoints = [oldLS(2, :); oldLS(n - 1, :)];
        avgEnds   = 0.5 * (endPoints(1, :) + endPoints(2, :));

        % Move first point toward that average
        newLS(1, :) = oldLS(1, :) + factor .* (avgEnds - oldLS(1, :));
        % Make last point identical
        newLS(end, :) = newLS(1, :);
    end
end
