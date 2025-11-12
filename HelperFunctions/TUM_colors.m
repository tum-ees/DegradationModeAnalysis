%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-10
%
% TUM_colors
% Sets the default axes color order using a TUM-aligned palette inspired by
% https://gist.github.com/lnksz/51e3566af2df5c7aa678cd4dfc8305f7.

black       = [0 0 0];
tumBlueDark = [7 33 64]./255;         % #072140
tumBlue     = [48 112 179]./255;      % #3070B3 (TUM main blue)
darkBlue    = [0 82 147]./255;        % #005293
lightBlue   = [100 160 200]./255;     % #64A0C8
lighterBlue = [152 198 234]./255;     % #98C6EA
green       = [162 173 0]./255;       % #A2AD00
orange      = [227 114 34]./255;      % #E37222
gray        = [153 153 153]./255;     % #999999
mediumGray  = [106 117 126]./255;     % #6A757E (tum-grey-4)
darkGray    = [71 80 88]./255;        % #475058 (tum-grey-3)
lightGray   = [218 215 203]./255;     % #DAD7CB

tumColorOrder = [
    black;
    tumBlueDark;
    tumBlue;
    darkBlue;
    lightBlue;
    lighterBlue;
    green;
    orange;
    gray;
    mediumGray;
    darkGray;
    lightGray
];

set(groot, 'defaultAxesColorOrder', tumColorOrder);
