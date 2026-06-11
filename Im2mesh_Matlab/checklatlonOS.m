function checklatlonOS(lat, lon, func_name, lat_var_name, lon_var_name, lat_pos, lon_pos)
% checklatlonOS: Validate pair of latitude-longitude arrays
% Open-source version of function checklatlon
%
%   checklatlonOS(LAT, LON, FUNC_NAME, LAT_VAR_NAME, LON_VAR_NAME,
%   LAT_POS, LON_POS) ensures that LAT and LON are real vectors, arrays, 
%   or cell arrays of matching size, having class 'double' or 'single', 
%   and that for every NaN-valued element of LAT there is identically-positioned
%   NaN-valued element in LON, and vice versa.

    % 1. Handle Cell Arrays
    if iscell(lat) || iscell(lon)
        % Make sure both are cell arrays if one is
        assert(iscell(lat) && iscell(lon), ...
            "map:" + func_name + ":latlonCellMismatch", ...
            ['Function %s expected its %s and %s input arguments,\n', ...
            '%s and %s, to both be cell arrays or both be numeric arrays.'], ...
            upper(func_name), num2ordinal_local(lat_pos), num2ordinal_local(lon_pos), ...
            lat_var_name, lon_var_name)
        
        % Make sure the cell arrays are the same size
        assert(isequal(size(lat), size(lon)), ...
            "map:" + func_name + ":latlonSizeMismatch", ...
            ['Function %s expected its %s and %s input arguments,\n', ...
            '%s and %s, to match in size.'], ...
            upper(func_name), num2ordinal_local(lat_pos), num2ordinal_local(lon_pos), ...
            lat_var_name, lon_var_name)
        
        % Recursively check each cell contents
        for i = 1:numel(lat)
            check_numeric_latlon(lat{i}, lon{i}, func_name, lat_var_name, lon_var_name, lat_pos, lon_pos);
        end
        
    % 2. Handle standard Numeric Arrays / Vectors
    else
        check_numeric_latlon(lat, lon, func_name, lat_var_name, lon_var_name, lat_pos, lon_pos);
    end

end

% --- Helper Functions ---

function check_numeric_latlon(lat, lon, func_name, lat_var_name, lon_var_name, lat_pos, lon_pos)
    % Real-valued double or single?
    validateattributes( ...
        lat, {'double','single'} ,{'real'}, func_name, lat_var_name, lat_pos)
    validateattributes( ...
        lon, {'double','single'} ,{'real'}, func_name, lon_var_name, lon_pos)

    % Sizes match?
    assert(isequal(size(lat),size(lon)), ...
        "map:" + func_name + ":latlonSizeMismatch", ...
        ['Function %s expected its %s and %s input arguments,\n', ...
        '%s and %s, to match in size.'], ...
        upper(func_name), num2ordinal_local(lat_pos), num2ordinal_local(lon_pos), ...
        lat_var_name, lon_var_name)

    % NaN positions correspond?
    assert(isequal(isnan(lat),isnan(lon)), ...
        "map:" + func_name + ":latlonNaNMismatch", ...
        ['Function %s expected its %s and %s input arguments,\n', ...
        '%s and %s, to have NaN-separators in corresponding positions.'], ...
        upper(func_name), num2ordinal_local(lat_pos), num2ordinal_local(lon_pos), ...
        lat_var_name, lon_var_name)
end

function ordStr = num2ordinal_local(n)
    % Replicates Mapping Toolbox's num2ordinal function
    if ~isnumeric(n) || ~isscalar(n)
        ordStr = num2str(n); % Fallback
        return;
    end
    
    n = abs(round(n));
    rem100 = rem(n, 100);
    rem10  = rem(n, 10);
    
    % Determine the correct suffix (st, nd, rd, th)
    if rem100 >= 11 && rem100 <= 13
        suffix = 'th';
    elseif rem10 == 1
        suffix = 'st';
    elseif rem10 == 2
        suffix = 'nd';
    elseif rem10 == 3
        suffix = 'rd';
    else
        suffix = 'th';
    end
    
    ordStr = sprintf('%d%s', n, suffix);
end