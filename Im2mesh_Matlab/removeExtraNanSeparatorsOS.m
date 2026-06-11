function [xdata, ydata, zdata] = ...
    removeExtraNanSeparatorsOS(xdata, ydata, zdata)
% removeExtraNanSeparatorsOS:  Clean up NaN separators in polygons and lines
% Open-source version of function removeExtraNanSeparators
%
%   [XDATA, YDATA] = removeExtraNanSeparatorsOS(XDATA, YDATA) removes NaNs
%   from the vectors XDATA and YDATA, leaving only isolated NaN separators.
%   If present, one or more leading NaNs are removed entirely.  If present,
%   a single trailing NaN is preserved.  NaNs are removed, but never added,
%   so if the input lacks a trailing NaN, so will the output.  XDATA and
%   YDATA must match in size and have identical NaN locations.
%
%   [XDATA, YDATA, ZDATA] = removeExtraNanSeparatorsOS(XDATA, YDATA, ZDATA)
%   removes NaNs from the vectors XDATA, YDATA, and ZDATA, leaving only
%   isolated NaN separators and optionally, consistent with the input, a
%   single trailing NaN.
%
%   Examples
%   --------
%   xin = [NaN NaN 1:3 NaN 4:5 NaN NaN NaN 6:9 NaN NaN]
%   yin = xin;
%   [xout, yout] = removeExtraNanSeparatorsOS(xin, yin);
%   xout
%
%   xin = [NaN 1:3 NaN NaN 4:5 NaN NaN NaN 6:9]'
%   yin = xin;
%   zin = xin;
%   [xout, yout, zout] = removeExtraNanSeparatorsOS(xin, yin, zin);
%   xout
%
% Copyright 2005-2008 The MathWorks, Inc.

if nargin < 3
    if ~isequal(isnan(xdata), isnan(ydata))
        error(['map:' mfilename ':inconsistentXY'], ...
            'XDATA and YDATA mismatch in size or NaN locations.')
    end
else
    if ~isequal(isnan(xdata), isnan(ydata), isnan(zdata))
        error(['map:' mfilename ':inconsistentXYZ'], ...
            'XDATA, YDATA (or ZDATA) mismatch in size or NaN locations.')
    end
end

if ~isempty(xdata) && ~isvector(xdata)
    error(['map:' mfilename ':nonVectorInput'], ...
        'XDATA, YDATA (and ZDATA) must be vectors.')
end

% Locate extra NaNs.
n = isnan(xdata(:));
firstOrPrecededByNaN = [true; n(1:end-1)];
extraNaN = n & firstOrPrecededByNaN;

% Remove extra NaNs.
xdata(extraNaN) = [];
ydata(extraNaN) = [];
if nargin >= 3
    zdata(extraNaN) = [];
end
