function [cellIdentifier_keys, cellIdentifier_values] = extract_cell_identifier(tableFilter)
%> Authors : Mathias Rehm
%> E-mail  : mathias.rehm@tum.de
%> Date    : 2025-09-12
%
% extract_cell_identifier
% Parse one or multiple {"Key","Value"} filters.
% Accepts:
%   * 1x2 cell: {'Key','Value'}
%   * Nx2 cell: {'Key1','Val1'; 'Key2','Val2'; ...}
%   * struct : fields are Keys, field values are Values
% Returns:
%   * cellIdentifier_keys   - 1xN cell of char
%   * cellIdentifier_values - 1xN cell of char
%
% If empty input, returns {} and {}.

%% 1) Defaults
cellIdentifier_keys   = {};
cellIdentifier_values = {};

%% 2) Empty → nothing to filter
if isempty(tableFilter)
    return;
end

%% 3) Cell input
if iscell(tableFilter)
    % 3a) Single pair in a vector cell: {'Key','Value'}
    if isvector(tableFilter) && numel(tableFilter) == 2
        k = tableFilter{1}; v = tableFilter{2};
        if isempty(k) || isempty(v)
            error('TableFilter must include both Key and Value (non-empty).');
        end
        cellIdentifier_keys   = {char(k)};
        cellIdentifier_values = {char(v)};
        return;
    end
    % 3b) Multiple pairs in an Nx2 cell
    if size(tableFilter,2) == 2 && size(tableFilter,1) >= 1
        K = tableFilter(:,1);
        V = tableFilter(:,2);
        for i = 1:numel(K)
            if isempty(K{i}) || isempty(V{i})
                error('All Key/Value pairs must be non-empty (row %d).', i);
            end
        end
        cellIdentifier_keys   = cellfun(@char, K(:).', 'UniformOutput', false);
        cellIdentifier_values = cellfun(@char, V(:).', 'UniformOutput', false);
        return;
    end
    error('TableFilter cell input must be 1x2 or Nx2 (Key/Value) pairs.');
end

%% 4) Struct input
if isstruct(tableFilter)
    f = fieldnames(tableFilter);
    if isempty(f), return; end
    cellIdentifier_keys = cellfun(@char, f(:).', 'UniformOutput', false);
    vals = struct2cell(tableFilter);
    cellIdentifier_values = cellfun(@char, vals(:).', 'UniformOutput', false);
    return;
end

%% 5) Anything else → error
error('Unsupported TableFilter type. Use 1x2/Nx2 cell or struct.');
end
