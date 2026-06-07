function tf = ispolycwOS(x, y)
% ispolycwOS: True if polygon vertices are in clockwise order
% Open-source version of function ispolycw
%
%   TF = ispolycwOS(X, Y) returns true if the closed polygonal contour 
%   represented by X and Y is ordered in the clockwise direction. 
%
%   Assumptions Enforced:
%   1. X and Y are numeric vectors without NaN values.
%   2. X and Y represent a single polygon.
%   3. The polygon is closed (x(1) == x(end) && y(1) == y(end)).

    % ---------------------------------------------------------------------
    % Check inputs
    % ---------------------------------------------------------------------
    % 1. Check if inputs are numeric vectors
    if ~isnumeric(x) || ~isnumeric(y)
        error('ispolycwOS:InvalidDataType', ...
            'Inputs X and Y must be numeric.');
    end
    if ~isvector(x) || ~isvector(y)
        error('ispolycwOS:NotAVector', 'Inputs X and Y must be vectors.');
    end
    
    % 2. Check if X and Y match in size
    if numel(x) ~= numel(y)
        error('ispolycwOS:SizeMismatch', ...
            'Inputs X and Y must have the same number of elements.');
    end

    % Handle empty arrays gracefully (matches original MATLAB behavior)
    if isempty(x)
        tf = true;
        return;
    end

    % 3. Check for NaN values (Ensures it is a single, continuous polygon)
    if any(isnan(x)) || any(isnan(y))
        error('ispolycwOS:ContainsNaN', ...
            ['Inputs X and Y cannot contain NaN values. ' ...
            'This function only supports a single polygon.']);
    end

    % 4. Check if the polygon is closed
    if x(1) ~= x(end) || y(1) ~= y(end)
        error('ispolycwOS:NotClosed', ...
            ['The polygon must be closed. The first and last vertices' ...
            ' must be identical (x(1)==x(end) and y(1)==y(end)).']);
    end

    % ---------------------------------------------------------------------
    % Core logic
    % ---------------------------------------------------------------------
    % Ensure inputs are column vectors to avoid dimension mismatches
    x = x(:);
    y = y(:);

    % Polygons with 2 or fewer unique vertices (<= 3 elements since it is 
    % closed) do not have a well-defined orientation. The original function
    % returns true.
    if numel(x) <= 3
        tf = true;
        return;
    end

    % Calculate the signed area using the modified Shoelace Formula.
    signedArea = sum(x(2:end) .* y(1:end-1) - x(1:end-1) .* y(2:end));

    % A result >= 0 indicates a clockwise orientation
    tf = (signedArea >= 0);

    % ---------------------------------------------------------------------
end