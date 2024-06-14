function Leaves = generate_foliage_qsm_direct(QSM, ...
                                              TargetDistributions, ...
                                              LeafProperties, ...
                                              totalLeafArea, ...
                                              varargin)

%% Check the correctness of inputs
QSM = check_inputs_direct(QSM,TargetDistributions,LeafProperties, ...
                          totalLeafArea);

%% Optional input default values
intersectionPrevention = true;
overSamplingFactor = 2;
TwigDirectionDistribution.flag = false;
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

            case 'twigdirectiondistribution'
                assert(i < NArg && isa(varargin{i+1},'function_handle'),...
                       'Argument following ''TwigDirectionDistribution'' should be a function handle.');
                TwigDirectionDistribution.flag = true;
                TwigDirectionDistribution.dist_fun = varargin{i+1};
                i = i + 1;

            case 'phyllotaxis'
                assert(i < NArg && isa(varargin{i+1},'struct'), ...
                       'Argument following ''Phyllotaxis'' should be a struct.')
                Phyllotaxis = varargin{i+1};
                Phyllotaxis.flag = true;
                if TwigDirectionDistribution.flag == true
                    warning('Twig direction distribution cannot be used simultaneously with phyllotaxis enabled')
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

% Scale by target area to get cylinder leaf area budgets
cylinderLeafArea = totalLeafArea*relativeCylinderLeafArea;

%% Initialize leaf object and candidate leaf area
Leaves = LeafModelTriangle(LeafProperties.vertices, ...
                           LeafProperties.triangles);
baseArea = Leaves.base_area;
candidateArea = overSamplingFactor*totalLeafArea;
cylinderCandidateLeafArea = candidateArea*relativeCylinderLeafArea;

%% Sample leaf sizes
[leafScaleFactors,leafParent] = sample_leaf_sizes(Leaves, ...
                                                CylinderParameters, ...
                                                TargetDistributions, ...
                                                cylinderCandidateLeafArea);

%% Sample leaf orientations
[leafDir,leafNormal,twigStart,twigEnd] = sample_leaf_orientations(...
                                             CylinderParameters,...
                                             leafParent,...
                                             TargetDistributions,...
                                             LeafProperties,...
                                             TwigDirectionDistribution,...
                                             Phyllotaxis);

%% Randomize the order of leaves
randOrder = randperm(size(leafScaleFactors,1));
leafScaleFactors = leafScaleFactors(randOrder,:);
leafParent       = leafParent(randOrder,:);
leafDir          = leafDir(randOrder,:);
leafNormal       = leafNormal(randOrder,:);
twigStart        = twigStart(randOrder,:);
twigEnd          = twigEnd(randOrder,:);

%% Add leaves on QSM
if intersectionPrevention == true
    Leaves = add_leaves_qsm(QSM,Leaves,leafScaleFactors,leafParent,...
                            leafDir,leafNormal,twigStart,twigEnd,...
                            totalLeafArea);
else
    iLeaf = 1;
    totalArea = 0;
    nLeaves = size(leafScaleFactors,1);
    while totalArea < totalLeafArea && iLeaf < nLeaves
        Leaves.add_leaf(twigEnd(iLeaf,:),leafDir(iLeaf,:),...
                        leafNormal(iLeaf,:),leafScaleFactors(iLeaf,:), ...
                        1,twigStart(iLeaf,:));
        totalArea = totalArea ...
                    + (leafScaleFactors(iLeaf,1)^2)*Leaves.base_area;
        iLeaf = iLeaf + 1;
    end
end

end