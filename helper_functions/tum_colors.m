function colors = tum_colors()
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-10
%
% tum_colors
% Return a TUM-aligned color palette and default plotting order.
%
% INPUTS:
%   none
%
% OUTPUTS:
%   colors - struct with named RGB triplets in [0, 1] and a `colorOrder`
%            field for axes styling

colors.black       = [0 0 0];
colors.tumBlueDark = [7 33 64]./255;         % #072140
colors.tumBlue     = [48 112 179]./255;      % #3070B3 (TUM main blue)
colors.darkBlue    = [0 82 147]./255;        % #005293
colors.lightBlue   = [100 160 200]./255;     % #64A0C8
colors.lighterBlue = [152 198 234]./255;     % #98C6EA
colors.green       = [162 173 0]./255;       % #A2AD00
colors.orange      = [227 114 34]./255;      % #E37222
colors.gray        = [153 153 153]./255;     % #999999
colors.mediumGray  = [106 117 126]./255;     % #6A757E (TUM gray 4)
colors.darkGray    = [71 80 88]./255;        % #475058 (TUM gray 3)
colors.lightGray   = [218 215 203]./255;     % #DAD7CB

colors.colorOrder = [
    colors.black;
    colors.tumBlueDark;
    colors.tumBlue;
    colors.darkBlue;
    colors.lightBlue;
    colors.lighterBlue;
    colors.green;
    colors.orange;
    colors.gray;
    colors.mediumGray;
    colors.darkGray;
    colors.lightGray
];
end
