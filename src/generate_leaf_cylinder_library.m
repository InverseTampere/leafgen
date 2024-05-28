function LeafCylLib = generate_leaf_cylinder_library(Nodes, ...
                                                     LibraryDistributions, ...
                                                     LeafProperties, ...
                                                     varargin)

%% Default values

nLeafObjectsPerNode = 3;
intersectionPrevention = false;
overSamplingFactor = 2;
TwigDirectionDistribution.flag = false;
Phyllotaxis.flag = false;

if LibraryDistributions.dTypeLODinc == "uniform"
    Nodes.pLODinc1 = 0;
    Nodes.pLODinc2 = 0;
end

if LibraryDistributions.dTypeLODaz == "uniform"
    Nodes.pLODaz1 = 0;
    Nodes.pLODaz2 = 0;
end

%% Read inputs

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg

    if ischar(varargin{i})

        switch lower(varargin{i})

            case 'nleafobjectspernode'
                assert(i < NArg && isnumeric(varargin{i+1}) && ...
                       isscalar(varargin{i+1}) && varargin{i+1} > 0, ...
                       'Argument following ''NLeafObjectsPerNode'' should be a positive integer.');
                nLeafObjectsPerNode = varargin{i+1};
                i = i + 1;

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

%% Library metadata

% Leaf distributions of the library
LeafCylLib.LeafDistributions = LibraryDistributions;

% Library nodes
LeafCylLib.Nodes = Nodes;

% Twig length limits
LeafCylLib.twigLengthLimits = LeafProperties.twigLengthLimits;

% Leaf base information
LeafCylLib.LeafBaseModel.vertices = LeafProperties.vertices;
LeafCylLib.LeafBaseModel.tris     = LeafProperties.triangles;

% Number of leaf objects generated for each node
LeafCylLib.nLeafObjectsPerNode = nLeafObjectsPerNode;

%% Leaf-cylinder generation

% Ranges of library variables
nLODinc1 = length(Nodes.pLODinc1);
nLODinc2 = length(Nodes.pLODinc2);
nLODaz1  = length(Nodes.pLODaz1);
nLODaz2  = length(Nodes.pLODaz2);
nLSD1    = length(Nodes.pLSD1);
nLSD2    = length(Nodes.pLSD2);

nLen = length(Nodes.cylinderLength);
nRad = length(Nodes.cylinderRadius);
nInc = length(Nodes.cylinderInclinationAngle);
nAz  = length(Nodes.cylinderAzimuthAngle);
nAr  = length(Nodes.cylinderLeafArea);

% Number of nodes per library variable
nNodesPerLibVar = [nLen,nRad,nInc,nAz,nAr,nLODinc1,nLODinc2,nLODaz1,...
                   nLODaz2,nLSD1,nLSD2,nLeafObjectsPerNode];

% Initializing the cell for leaf-cylinder library (empty cell element added
% to dimensions with only one node to ensure 12-dimensional size for the 
% cell variable)
LeavesObjects = cell(nNodesPerLibVar + (nNodesPerLibVar==1));

% Total number of iterations
nIter = prod(nNodesPerLibVar);
iIter = 0;
wb = waitbar(iIter/nIter,num2str(iIter/nIter),...
             'Name','Progress bar');

% Create Leaves objects for the library
for  iLODinc1 = 1:nLODinc1
    for iLODinc2 = 1:nLODinc2
    for iLODaz1  = 1:nLODaz1
    for iLODaz2  = 1:nLODaz2
    for iLSD1    = 1:nLSD1
    for iLSD2    = 1:nLSD2

    for iCyl = 1:nLeafObjectsPerNode % loop over the desired amount of 
                                       % leaf-cylinders for each node

    % Loop over cylinder attributes
    for iLen = 1:nLen
    len = Nodes.cylinderLength(iLen);

    for iRad = 1:nRad
    rad = Nodes.cylinderRadius(iRad);

    for iInc = 1:nInc
    inc = Nodes.cylinderInclinationAngle(iInc);

    for iAz  = 1:nAz
    az  = Nodes.cylinderAzimuthAngle(iAz);

    for iAr  = 1:nAr
    ar  = Nodes.cylinderLeafArea(iAr);

    % Waitbar update
    iIter = iIter + 1;
    waitbar(iIter/nIter,wb, ...
            sprintf('Leaf cylinder generation in progress\n%d/%d', ...
                    iIter,nIter) ...
           );

    % If cylinder inclination is vertical pick already made leaf cylinder
    % if possible (cylinder azimuth loses significance for vertical 
    % cylinder)
    if (iInc == 1 || iInc == nInc) && iAz > 1
        Leaves = LeavesObjects{iLen,iRad,iInc,1,iAr,iLODinc1,iLODinc2, ...
                               iLODaz1,iLODaz2,iLSD1,iLSD2,iCyl};
        LeavesObjects{iLen,iRad,iInc,iAz,iAr,iLODinc1,iLODinc2,iLODaz1, ...
                      iLODaz2,iLSD1,iLSD2,iCyl} = Leaves;
        clearvars Leaves;
        % Skip to next iteration of the active for-loop
        continue
    end


    % Initialize Leaves object
    Leaves = LeafModelTriangle(LeafProperties.vertices, ...
                               LeafProperties.triangles);

    % Sample leaves from leaf size function
    [leafScaleFactors,nLeaves,maxLeafSize] = fun_leaf_size( ...
        overSamplingFactor*ar, ...
        Leaves.base_area, ...
        LibraryDistributions.dTypeLSD, ...
        [Nodes.pLSD1(iLSD1), Nodes.pLSD2(iLSD2)] ...
        );

    % Attach the leaves to the cylinder with leaf orientation
    % distribution
    [leafDir,leafNormal,twigStart,twigEnd] = fun_leaf_orientation( ...
        len,rad,inc,az, ...
        nLeaves, ...
        LeafProperties.twigLengthLimits, ...
        LibraryDistributions.dTypeLODinc, ...
        LibraryDistributions.dTypeLODaz, ...
        [Nodes.LODinc1(iLODinc1), Nodes.LODinc2(iLODinc2)], ...
        [Nodes.LODaz1(iLODaz1), Nodes.LODaz2(iLODaz2)], ...
        TwigDirectionDistribution, ...
        Phyllotaxis ...
        );

    % Add leaves to the model
    if intersectionPrevention == true
        Leaves = add_leaves(Leaves,len,rad,inc,az,ar,twigStart,twigEnd, ...
                            leafDir,leafNormal,leafScaleFactors, ...
                            maxLeafSize, ...
                            max(LeafProperties.twigLengthLimits));
    else
        iLeaf = 1;
        totalArea = 0;
        while totalArea < ar && iLeaf < nLeaves
            Leaves.add_leaf(twigEnd(iLeaf,:), ...
                leafDir(iLeaf,:), ...
                leafNormal(iLeaf,:), ...
                leafScaleFactors(iLeaf,:), ...
                1, ...
                twigStart(iLeaf,:) ...
                );
            totalArea = leafScaleFactors(iLeaf,1)*leafScaleFactors(iLeaf,2) ...
                        *Leaves.base_area;
            iLeaf = iLeaf + 1;
        end
    end

    % Add Leaves object to the node cell
    LeavesObjects{iLen,iRad,iInc,iAz,iAr,iLODinc1,iLODinc2,iLODaz1, ...
                  iLODaz2,iLSD1,iLSD2,iCyl} = Leaves;
    clearvars Leaves;

    end % cylinderLeafArea
    end % cylinderAzimuthAngle
    end % cylinderInclinationAngle
    end % cylinderRadius
    end % cylinderLenght

    end % nLeafCylindersPerNode

    end % LSD2
    end % LSD1
    end % LODaz2
    end % LODaz1
    end % LODinc2
end     % LODinc1

LeafCylLib.LeavesObjects = LeavesObjects;

% Delete waitbar
delete(wb)

%
end


