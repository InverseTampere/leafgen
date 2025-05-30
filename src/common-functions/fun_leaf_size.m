% This file is part of LeafGen
% 
% LeafGen is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% LeafGen is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with LeafGen.  If not, see <https://www.gnu.org/licenses/>.

function leafScaleFactors = fun_leaf_size(targetArea,baseArea, ...
                                          dType,dParameters)

% Initialize number of generated leafs
nLeaves = 0;

% Sample leaves from the size distribution
if strcmp(dType,'constant')
    % Find the amount of leaves needed for target area (dParameters first
    % element contains the constant area of a single leaf)
    constantArea = dParameters(1);
    nLeaves = round(targetArea/constantArea);
    % Store leaf scaling factor (same for all dimensions)
    leafScaleFactors = sqrt(constantArea/baseArea)*ones(nLeaves,3);
else
    % Calculate the estimated size for the output vector
    switch dType
        case 'uniform'
            nLeafEstimate = ceil(2*targetArea/mean(dParameters));
        case 'normal'
            nLeafEstimate = ceil(2*targetArea/dParameters(1));
    end
    % Initialize the output vector
    leafScaleFactors = zeros(nLeafEstimate,3);
    % Sample leaves until target area is reached
    areaAdded = 0;
    while areaAdded < targetArea
        switch dType
            case 'uniform'
                sampledArea = (dParameters(2)-dParameters(1))*rand(1) ...
                              + dParameters(1);
            case 'normal'
                sampledArea = sqrt(dParameters(2))*randn(1) ...
                              + dParameters(1);
                % Resample negative and zero values
                while sampledArea <= 0
                    sampledArea = sqrt(dParameters(2))*randn(1) ...
                                  + dParameters(1);
                end
        end
        areaAdded = areaAdded + sampledArea;
        % Store leaf scaling factor (same for all dimensions)
        leafScaleFactors(nLeaves+1,:) = sqrt(sampledArea/baseArea)...
                                        *ones(1,3);
        nLeaves = nLeaves + 1;
    end
    % Remove empty rows from the end of output vector
    leafScaleFactors = leafScaleFactors(1:nLeaves,:);
end

end