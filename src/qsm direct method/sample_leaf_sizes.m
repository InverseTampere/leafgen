function [leafScaleFactors,leafParent] = sample_leaf_sizes( ...
                                              Leaves, ...
                                              CylinderParameters, ...
                                              TargetDistributions, ...
                                              cylinderCandidateLeafArea)
% Leaf base area
baseArea = Leaves.base_area;

% LSD type and parameter function
dType = TargetDistributions.dTypeLSD;
fun_param = TargetDistributions.fun_size_params;

% Total number of cylinders
nCyl = length(cylinderCandidateLeafArea);

% Cylinder structural variables
relH = CylinderParameters.relative_height;
relD = CylinderParameters.relative_distance_along_subbranch;
cDir = CylinderParameters.compass_direction;

% Leaf areas to be generated on cylinders
areaMissingPerCyl = cylinderCandidateLeafArea;

% Initialize output vectors
leafScaleFactors = zeros(10000,3);
leafParent = zeros(10000,1);

% Initialize loop variables
iLeaf = 0;
iCyl = 1;
targetReached = 0;
while targetReached == 0
    % Leaf size distribution parameters for cylinder
    params = fun_param(relH(iCyl),relD(iCyl),cDir(iCyl));
    % Mean leaf
    switch dType
        case 'constant'
            meanArea = params;
        case 'uniform'
            meanArea = (params(2)-params(1))/2;
        case 'normal'
            meanArea = params(1);
    end

    % Add leaf if necessary
    if areaMissingPerCyl(iCyl) >= meanArea
        % Increment leaf count
        iLeaf = iLeaf + 1;
        % Sample leaf size
        switch dType
            case 'constant'
                sampledArea = params;
            case 'uniform'
                sampledArea = (params(2)-params(1))*rand(1) + params(1);
            case 'normal'
                sampledArea = sqrt(params(2))*randn(1) + params(1);
        end
        leafScaleFactors(iLeaf,:) = sqrt(sampledArea/baseArea)*ones(1,3);
        leafParent(iLeaf) = iCyl;
        areaMissingPerCyl(iCyl) = areaMissingPerCyl(iCyl) - sampledArea;
    else
        iCyl = iCyl + 1;
    end

    % Check if target area is reached
    if iCyl == nCyl && areaMissingPerCyl(end) < meanArea
        targetReached = 1;
        nLeaves = iLeaf;
    end
    
    % Extend output variables if necessary
    if iLeaf >= size(leafScaleFactors,1) && targetReached == 0
        leafScaleFactors = [leafScaleFactors; zeros(10000,3)]; %#ok<AGROW>
        leafParent       = [leafParent;       zeros(10000,1)]; %#ok<AGROW>
    end
end
% Remove empty rows from the end of output vector
leafScaleFactors = leafScaleFactors(1:nLeaves,:);
leafParent       = leafParent(1:nLeaves,:);

end