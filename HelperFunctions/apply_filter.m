% %%> Author: Mathias Rehm 
%> Email: mathias.rehm@tum.de
%> Date: 2023-08-24

%% apply_filter
% Requirements:
%   * A vector (or array) of data that you want to smooth or filter.
%   * One or more filter methods can be specified via 'filtermethod' and 
%     related arguments in varargin.
%
% Overview:
%   * The function applies one or multiple filters (movmean, movmedian,
%     gaussian, sgolay, etc.) to the input data, handling NaNs 
%     appropriately.
%   * The length of the filter window can be specified as a fixed number 
%     ('framelenfixed', N) or as a fraction of the data length 
%     ('framelenfactor', f).
%   * Some optional steps like outlier handling ('filloutlier') can be 
%     activated.
%
% Setting Options (in varargin):
%   * 'filtermethod': ('movmean', 'movmedian', 'gaussian', 'sgolay') – 
%       multiple calls can be used to layer filters in sequence.
%   * 'framelenfixed' or 'framelenfactor': sets the window size. E.g. 
%       'framelenfixed',50 or 'framelenfactor',0.005.
%   * 'order': for sgolay filter (default = 3).
%   * 'repeat': number of times each filter is repeated (default = 1).
%   * 'filloutlier': if set to 1 or true, uses MATLAB's filloutliers(...) to 
%       replace outliers with 'linear' interpolation before filtering.
%
% Function calculates:
%   * An output vector (or array) 'output' that is the filtered version of 
%     the input data.
%
% Example:
%   yFiltered = apply_filter(mySignal, 'filtermethod','sgolay', ...
%                               'framelenfixed', 30, 'repeat', 2);
%   % Applies a moving average of window size 30, repeated twice.

function [output] = apply_filter(data, varargin)

% Default values for filter settings
defaultVal.framelen.movmean   = 10;
defaultVal.framelen.movmedian = 10;
defaultVal.framelen.gaussian  = 10;
defaultVal.framelen.sgolay    = 50;
defaultVal.repeatFilter.num   = 1;      % Default number of filter repeats
defaultVal.order.sgolay       = 3;      % Default order for sgolay

dispAppliedFilters = 0;  % If 1, displays which filters were applied

% Preallocate space
filtermethod = {};
framelen     = {};
order        = {};
repeatFilter = {};

fillOutlier = 0;

% Convert numeric to string in varargin
for i = 1 : length(varargin)
    if isnumeric(varargin{i})
        varargin{i} = num2str(varargin{i});
    end
end

% Attempt to lowercase all varargin
try 
    varargin = lower(varargin);
catch
    disp('Please use only strings as input for varargin for apply_filter!')
end

% Remove empty entries in varargin
for i = length(varargin) : -1 : 1
    if isempty(varargin{i})
        try
            varargin(i) = [];
            varargin(i-1) = [];
        catch
            disp('First entry in varargin cannot be empty!')
        end
    end
end

% Parse varargin
for i = 1 : 2 : length(varargin)
    switch true
        % Define filtermethod
        case (contains(varargin{i},'filtermethod','IgnoreCase',true) == 1)
            filterStr = varargin{i+1};
            switch true
                case (contains(filterStr,'mean','IgnoreCase',true) == 1)
                    filtermethod{end+1} = 'movmean';
                    framelen{length(filtermethod)} = defaultVal.framelen.movmean;
                case (contains(filterStr,'median','IgnoreCase',true) == 1)
                    filtermethod{end+1} = 'movmedian';
                    framelen{length(filtermethod)} = defaultVal.framelen.movmedian;                    
                case (contains(filterStr,'gaus','IgnoreCase',true) == 1 || ...
                      contains('gauß',filterStr,'IgnoreCase',true) == 1)
                    filtermethod{end+1} = 'gaussian';
                    framelen{length(filtermethod)} = defaultVal.framelen.gaussian;                      
                case (contains(filterStr,'sgola','IgnoreCase',true) == 1 || ...
                      contains('golay',filterStr,'IgnoreCase',true) == 1 || ...
                      contains('savitzk',filterStr,'IgnoreCase',true) == 1)
                    filtermethod{end+1} = 'sgolay';
                    framelen{length(filtermethod)} = defaultVal.framelen.sgolay;
                    order{length(filtermethod)}    = defaultVal.order.sgolay;
            end

        % Define length of filter window
        case (contains(varargin{i}, 'framel', 'IgnoreCase',true) == 1)
            framelenStr = varargin{i};
            if (contains(framelenStr, 'fix', 'IgnoreCase',true) == 1)
                framelen{length(filtermethod)} = round(str2double(varargin{i+1}));
            elseif (contains(framelenStr, 'fact', 'IgnoreCase',true) == 1)
                framelen{length(filtermethod)} = round(str2double(varargin{i+1}) * length(data));
            else
                disp(['Recommended to use either "framelenfactor" or "framelenfixed". ', ...
                      'Otherwise, it is automatically selected.'])
                % Decide automatically based on value
                if str2double(varargin{i+1}) < 1
                    framelen{length(filtermethod)} = round(str2double(varargin{i+1}) * length(data));
                elseif str2double(varargin{i+1}) > 1
                    framelen{length(filtermethod)} = round(str2double(varargin{i+1}));
                else
                    disp('Automatic selection between fixed framelength and factor did not work!')
                end
            end

        % Sgolay order
        case (contains(varargin{i}, 'order', 'IgnoreCase',true) == 1)
            order{length(filtermethod)} = round(str2double(varargin{i+1}));

        % Option to apply filter multiple times
        case (contains(varargin{i}, 'repeat', 'IgnoreCase',true) == 1)
            repeatFilter{end+1} = round(str2double(varargin{i+1}));

        % Fill outliers
        case ( (contains(varargin{i}, 'fill', 'IgnoreCase', true)  || ...
                contains(varargin{i}, 'delete', 'IgnoreCase', true)) && ...
                 contains(varargin{i}, 'outlier', 'IgnoreCase', true) )
            fillOutlier = varargin{i + 1};
            if ischar(fillOutlier) || isstring(fillOutlier)
                fillOutlier = logical(str2double(fillOutlier));
            elseif isnumeric(fillOutlier)
                fillOutlier = logical(fillOutlier);
            end
    end
end

% If no filtermethod is defined, default to 'sgolay'
if isempty(filtermethod)
    filtermethod{end+1} = 'sgolay';
    framelen{length(filtermethod)} = defaultVal.framelen.sgolay;
    order{length(filtermethod)} = defaultVal.order.sgolay;
    warning('No filtermethod defined. Using Savitzky-Golay as default.')
end

% Check repeatFilter usage
if length(repeatFilter) < length(filtermethod)-1 && length(repeatFilter) > 1
    disp('Define repeatFilter for every used filtermethod or for none (default is repeat = 1).')
end

% If needed, fill up repeatFilter so it matches the number of filters
while length(repeatFilter) < length(filtermethod)
    repeatFilter{end+1} = defaultVal.repeatFilter.num;
end

% If requested, fill outliers before applying filters
if fillOutlier
    try
        data = filloutliers(data, 'linear', 'movmedian', 3, 'ThresholdFactor', 80);
    catch
        warning('filloutliers could not be performed.')
    end
end

% --- Perform filtering ---
for i = 1 : length(filtermethod)
    idx = isnan(data);
    dataNan = data(idx);
    dataNotNan = data(~idx);

    switch filtermethod{i}
        case {'movmean', 'movmedian', 'gaussian'}
            for ii = 1 : repeatFilter{i}
                dataNotNan = smoothdata(dataNotNan, filtermethod{i}, framelen{i});
            end

        case 'sgolay'
            framelenNow = max(floor(framelen{i} / 2) * 2 + 1, order{i} + 1);
            if mod(framelenNow, 2) == 0
                framelenNow = framelenNow + 1;
            end
            for ii = 1 : repeatFilter{i}
                try
                    dataNotNan = sgolayfilt(dataNotNan, order{i}, framelenNow);
                catch
                    warning('Filtering with Savitzky-Golay did not work!')
                end
            end
    end

    % Recombine the NaN and filtered portions
    try
        reconnectedData = data; 
        reconnectedData(idx) = dataNan;
        reconnectedData(~idx) = dataNotNan;
        data = reconnectedData;
    catch
        data = [dataNan(:); dataNotNan(:)];
    end
end

% Prepare output
strFiltermethods = [];
for idx = 1 : length(filtermethod)
    if ~isempty(order) && idx <= length(order) && ~isempty(order{idx})
        strFiltermethods = [strFiltermethods, filtermethod{idx}, ...
                            '(framelen: ', num2str(framelen{idx}), ...
                            '; order: ', num2str(order{idx}), ') '];
    else
        strFiltermethods = [strFiltermethods, filtermethod{idx}, ...
                            '(framelen: ', num2str(framelen{idx}), ...
                            '; order: - ) '];
    end
end

if exist("dispAppliedFilters", "var") && dispAppliedFilters == 1
    disp(['Filters applied: ', strFiltermethods]);
end

output = data;
end
