% TEMPLATE FOR CUSTOM FEATURE
% -------------------------------
% 
% This script is a template for a custom feature.
% Please notice that the path to this feature has to be provided to the GUI
% before it is computed.

% --- computation feature

% Perform the computation of the feature. The coordinates of the
% trajectories are stored in a column vector cell array, which can be 
% addressed as obj(ii).coords{jj} with jj the index of the trajectory. The 
% output is a matrix with [1 .. # localizations, 1 .. #dimensions]. The 
% corresponding times can be retrieved as obj(ii).time{jj} with jj the 
% index of the trajectory. The output is a column vector with the frames.
% The unit in which the coordinates are stored as the same as the imported
% input file with trajectories. This is usually in pixels (coords) and
% frames (time).
% 
% Please do not use the following variables in this script:
% - preserveFlag
% - featureNamesFlag
% - filepathFeatures
% - wb
% - nObj
% - nDim
% - T
% - do

% For example, we compute the numer of frames the trajectory spans.
featureNameData = nan(size(obj(ii).time,1),1);
for jj = 1:size(obj(ii).time,1)
    featureNameData(jj) = (obj(ii).time{jj}(end) - obj(ii).time{jj}(1)) + 1;
end

% Alternatively, let's compute the center of mass of the coordinates. We can
% speed up the computation by using 'cellfun'.
% featureNameData = cellfun(@(x) mean(x,1),obj(ii).coords);

% --- write to table

% Save the feature to the table with features. Make sure that
% featureNameData is a numeric vector. The featureName is the name of the
% feature as it is stored in the table.
T.featureName = featureNameData;

% Store the units of the feature. Remove this line if the feature has no
% units. The units have to be given in terms of the 'pixelsize', which is
% the size of a single pixel in the trajectory coordinates, and in terms of
% 'dt', which is the time between two subsequent frames. For instance for
% the units of the diffusion constant one would write 'pixelsize.^2./dt'.
T.Properties.VariableUnits{end} = 'pixelsize.^2./dt';

% Add a short description of the feature between ' '. Description may
% contain spaces.
T.Properties.VariableDescriptions{end} = 'Short description of feature';

% multiple feature can be computed in a single script when they are pasted
% below each other.