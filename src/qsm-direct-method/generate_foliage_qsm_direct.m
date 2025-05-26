function [Leaves,QSM] = generate_foliage_qsm_direct(qsm, ...
                                              TargetDistributions, ...
                                              LeafProperties, ...
                                              totalLeafArea, ...
                                              varargin)

%% Initialize struct-based QSM as QSMBCylindrical-class object

if isa(qsm,'struct')
    QSM = QSMBCylindrical(qsm);
elseif isa(qsm,'QSMBCylindrical')
    QSM = qsm;
else
    error("Invalid QSM type. Input QSM should be a struct or " + ...
          "a QSMBCylindrical-class object.")
end

%% Check the correctness of inputs

if isfield(LeafProperties,'twigLengthLimits')
    LeafProperties.petioleLengthLimits = LeafProperties.twigLengthLimits;
    LeafProperties = rmfield(LeafProperties,'twigLengthLimits');
end

QSM = check_inputs_direct(QSM,TargetDistributions,LeafProperties, ...
                          totalLeafArea);

%% Optional input default values
intersectionPrevention = true;
overSamplingFactor = 2;
PetioleDirectionDistribution.flag = false;
Phyllotaxis.flag = false;

%% Read optional inputs

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg

    if ischar(varargin{i})

        switch lower(varargin{i})

            case 'preventintersections'
                assert(i < NArg && isa(varargin{i+1},'logical') && ...
                       isscalar(varargin{i+1}), ...
                       'Argument following ''PreventIntersections'' should be a boolean.')
                intersectionPrevention = varargin{i+1};
                i = i + 1;

            case 'oversamplingfactor'
                assert(i < NArg && isnumeric(varargin{i+1}) && ...
                       isscalar(varargin{i+1}) && varargin{i+1}>1, ...
                       'Argument following ''OverSamplingFactor'' should be a scalar above the value of 1.');
                overSamplingFactor = varargin{i+1};
                i = i + 1;

            case 'petioledirectiondistribution'
                assert(i < NArg && isa(varargin{i+1},'function_handle'),...
                       'Argument following ''PetioleDirectionDistribution'' should be a function handle.');
                PetioleDirectionDistribution.flag = true;
                PetioleDirectionDistribution.dist_fun = varargin{i+1};
                i = i + 1;

            case 'phyllotaxis'
                assert(i < NArg && isa(varargin{i+1},'struct'), ...
                       'Argument following ''Phyllotaxis'' should be a struct.')
                Phyllotaxis = varargin{i+1};
                Phyllotaxis.flag = true;
                if PetioleDirectionDistribution.flag == true
                    warning('Petiole direction distribution cannot be used simultaneously with phyllotaxis enabled')
                end

            otherwise
                warning(['Skipping unknown parameters: ''' varargin{i} '''']);
        end
    end
    i = i + 1;
end

%% Extract QSM properties
ParLabels = {'relative_height',...
             'relative_position',...
             'branch_order',...
             'axis',...
             'length',...
             'radius',...
             'is_last',...
             'start_point',...
             'branch_index',...
             'index_in_branch'};
CylinderParameters = QSM.get_block_properties(ParLabels);

% Cut cylinders to shorter ones to increase the discretization on relative
% distance on subbranch (when phyllotaxis is not in use)
if Phyllotaxis.flag == false
    [CylinderParameters,originalIndex] = preprocess_cylinders( ...
                                            CylinderParameters, ...
                                            [0 0.1], ...
                                            "relative branch");
else
    originalIndex = cumsum(ones(QSM.block_count,1));
end

% Calculate the relative distance of cylinders along subbranch
CylinderParameters.relative_distance_along_subbranch = ...
    calculate_relative_distance_along_subbranch(CylinderParameters);

% Calculate compass direction of cylinders with respect to mean value of
% horizontal cylinder location
CylinderParameters.compass_direction = ...
    calculate_compass_direction(CylinderParameters);

%% Determining leaf area budgets for cylinders using LADD

% Relative leaf area budgets
relativeCylinderLeafArea = fun_leaf_area_density(CylinderParameters, ...
                                                 TargetDistributions);

%% Initialize leaf base area and candidate leaf area
LeavesInit = LeafModelTriangle(LeafProperties.vertices, ...
                           LeafProperties.triangles);
candidateArea = overSamplingFactor*totalLeafArea;
cylinderCandidateLeafArea = candidateArea*relativeCylinderLeafArea;

%% Sample leaf sizes
[leafScaleFactors,leafParentPP] = sample_leaf_sizes(LeavesInit, ...
                                                CylinderParameters, ...
                                                TargetDistributions, ...
                                                cylinderCandidateLeafArea);

%% Sample leaf orientations
[leafDir,leafNormal,petioleStart,petioleEnd] = sample_leaf_orientations(...
                                             CylinderParameters,...
                                             leafParentPP,...
                                             TargetDistributions,...
                                             LeafProperties,...
                                             PetioleDirectionDistribution,...
                                             Phyllotaxis);

%% Randomize the order of leaves

randOrder = randperm(size(leafScaleFactors,1));
leafScaleFactors = leafScaleFactors(randOrder,:);
leafParent       = leafParentPP(randOrder,:);
leafDir          = leafDir(randOrder,:);
leafNormal       = leafNormal(randOrder,:);
petioleStart        = petioleStart(randOrder,:);
petioleEnd          = petioleEnd(randOrder,:);

%% Add leaves on QSM

% Average leaf area
avgAr = mean(LeavesInit.base_area*(leafScaleFactors(:,1).^2));
% Estimate on total leaf count
leafCountEst = int64(1.1*round(totalLeafArea/avgAr));
% Initialize leaf object
Leaves = LeafModelTriangle(LeafProperties.vertices, ...
                           LeafProperties.triangles, ...
                           0,leafCountEst);

% Set the correct leaf parent cylinder indices (i.e. indices before
% pre-processing of cylinders)
leafParent = zeros(size(leafParentPP));
for iLeaf = 1:size(leafScaleFactors,1)
    leafParent(iLeaf) = originalIndex(leafParentPP(iLeaf));
end

% Add leaves to the model
if intersectionPrevention == true
    Leaves = add_leaves_qsm(QSM,Leaves,leafScaleFactors,leafParent,...
                            leafDir,leafNormal,petioleStart,petioleEnd,...
                            totalLeafArea);
else
    iLeaf = 1;
    totalArea = 0;
    nLeaves = size(leafScaleFactors,1);
    while totalArea < totalLeafArea && iLeaf < nLeaves
        Leaves.add_leaf(petioleEnd(iLeaf,:),leafDir(iLeaf,:),...
                        leafNormal(iLeaf,:),leafScaleFactors(iLeaf,:), ...
                        1,petioleStart(iLeaf,:));
        totalArea = totalArea ...
                    + (leafScaleFactors(iLeaf,1)^2)*Leaves.base_area;
        iLeaf = iLeaf + 1;
    end
end

% Trim excess rows
Leaves.trim_slack;

end