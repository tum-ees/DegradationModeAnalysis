function [savePath_model, savePathDMA, savePathBlend2Content] = handle_paths(settings)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function:
%   * adds all necessary paths to conduct the DMA
%   * creates all save paths needed
%
%% ---------- ADD PATH TO PLOTTING FUNCTIONS ------------------------------
    % get aging study path from settings
    pathSave = settings.pathSaveResults;

%% --------------- CREATE SAVE PATHS --------------------------------------
    % Decide common path parts once, then reuse ---------------------------
    parts = {'01_Results'};                                   
    mk = @(leaf) fullfile(pathSave, parts{:}, leaf);  % tiny path builder
    
    % Final save paths
    savePath_model          = mk('Results_Parameter');
    savePathDMA             = mk('Results_DMA');
    savePathBlend2Content   = mk('Results_Blend2Content');
    
    % ensure directories exist
    for p = {savePath_model, savePathDMA, savePathBlend2Content}
        if ~isfolder(p{1}), mkdir(p{1}); end
    end

end

