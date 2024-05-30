function leafScaleFactors = fun_leaf_size(targetArea,baseArea, ...
                                          dType,dParameters)

% Initialize number of generated leafs
nLeaves = 0;

% Sample leaves from the size distribution
if strcmp(dType,'constant')
    % Find the amount of leaves needed for target area (dParameters
    % contains the constant area of a single leaf)
    nLeaves = round(targetArea/dParameters);
    % Store leaf scaling factor (same for all dimensions)
    leafScaleFactors = (constantArea/baseArea)*ones(nLeaves,3);
else
    % Calculate the estimated size for the output vector
    switch dType
        case 'uniform'
            nLeafEstimate = ceil(1.2*targetArea/mean(dParameters));
        case 'normal'
            nLeafEstimate = ceil(1.2*targetArea/dParameters(1));
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