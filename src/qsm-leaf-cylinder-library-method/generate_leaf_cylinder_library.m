function LeafCylLib = generate_leaf_cylinder_library(...
                                        LibraryDistributions, ...
                                        Nodes, ...
                                        LeafProperties, ...
                                        varargin)

%% Default values

nLeafObjectsPerNode = 3;
intersectionPrevention = true;
overSamplingFactor = 2;
PetioleDirectionDistribution.flag = false;
Phyllotaxis.flag = false;

if LibraryDistributions.dTypeLODinc == "uniform"
    Nodes.pLODinc1 = 0;
    Nodes.pLODinc2 = 0;
end
if LibraryDistributions.dTypeLODaz == "uniform"
    Nodes.pLODaz1 = 0;
    Nodes.pLODaz2 = 0;
end
if LibraryDistributions.dTypeLSD == "constant"
    Nodes.pLSD2 = 0;
end

%% Read inputs

% Check the correctness of inputs

if isfield(LeafProperties,'twigLengthLimits')
    LeafProperties.petioleLengthLimits = LeafProperties.twigLengthLimits;
    LeafProperties = rmfield(LeafProperties,'twigLengthLimits');
end

check_inputs_library(LibraryDistributions,Nodes,LeafProperties);

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg

    if ischar(varargin{i})

        switch lower(varargin{i})

            case 'nleafobjectspernode'
                assert(i < NArg && isnumeric(varargin{i+1}) && ...
                       isscalar(varargin{i+1}) && varargin{i+1} > 0, ...
                       "Argument following ''NLeafObjectsPerNode''"...
                       +" should be a positive integer.");
                nLeafObjectsPerNode = varargin{i+1};
                i = i + 1;

            case 'preventintersections'
                assert(i < NArg && isa(varargin{i+1},'logical') && ...
                       isscalar(varargin{i+1}), ...
                       "Argument following ''PreventIntersections''"...
                       +" should be a boolean.")
                intersectionPrevention = varargin{i+1};
                i = i + 1;

            case 'oversamplingfactor'
                assert(i < NArg && isnumeric(varargin{i+1}) && ...
                       isscalar(varargin{i+1}) && varargin{i+1}>1, ...
                       "Argument following ''OverSamplingFactor''"...
                       +" should be a scalar above the value of 1.");
                overSamplingFactor = varargin{i+1};
                i = i + 1;

            case 'petioledirectiondistribution'
                assert(i < NArg && isa(varargin{i+1},'function_handle'),...
                       "Argument following"...
                       +" ''PetioleDirectionDistribution'' should be a"...
                       +" function handle.");
                PetioleDirectionDistribution.flag = true;
                PetioleDirectionDistribution.dist_fun = varargin{i+1};
                i = i + 1;

            case 'phyllotaxis'
                assert(i < NArg && isa(varargin{i+1},'struct'), ...
                       "Argument following ''Phyllotaxis'' should be a"...
                       +" struct.")
                Phyllotaxis = varargin{i+1};
                Phyllotaxis.flag = true;
                if PetioleDirectionDistribution.flag == true
                    warning("Petiole direction distribution cannot be"...
                            +" used simultaneously with phyllotaxis"...
                            +" enabled")
                end

            otherwise
                warning("Skipping unknown parameters:"...
                        +" ''"+varargin{i}+"''");
        end
    end
    i = i + 1;
end

%% Library metadata

% Leaf distributions of the library
LeafCylLib.LeafDistributions = LibraryDistributions;

% Library nodes
LeafCylLib.Nodes = Nodes;

% Leaf base information
LeafCylLib.LeafBaseModel.vertices  = LeafProperties.vertices;
LeafCylLib.LeafBaseModel.triangles = LeafProperties.triangles;

% Other properties
LeafCylLib.Properties.petioleLengthLimits = ...
    LeafProperties.petioleLengthLimits;
LeafCylLib.Properties.nLeafObjectsPerNode = nLeafObjectsPerNode;
LeafCylLib.Properties.intersectionPrevention = intersectionPrevention;
LeafCylLib.Properties.overSamplingFactor = overSamplingFactor;
LeafCylLib.Properties.PetioleDirectionDistribution = ...
    PetioleDirectionDistribution;
LeafCylLib.Properties.Phyllotaxis = Phyllotaxis;


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

% Total number of nodes
LeafCylLib.totalNodes = prod(nNodesPerLibVar);

% Initializing the cell for leaf-cylinder library (empty cell element added
% to dimensions with only one node to ensure 12-dimensional size for the 
% cell variable)
LeafObjects = cell(nNodesPerLibVar + (nNodesPerLibVar==1));

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
    if all([(inc < 1e-6 || abs(inc-pi) < 1e-6), ...
            iAz > 1, ...
            PetioleDirectionDistribution.flag == false, ...
            Phyllotaxis.flag == false])
        Leaves = LeafObjects{iLen,iRad,iInc,1,iAr,iLODinc1,iLODinc2, ...
                               iLODaz1,iLODaz2,iLSD1,iLSD2,iCyl};
        LeafObjects{iLen,iRad,iInc,iAz,iAr,iLODinc1,iLODinc2,iLODaz1, ...
                      iLODaz2,iLSD1,iLSD2,iCyl} = Leaves;
        clearvars Leaves;
        % Skip to next iteration of the active for-loop
        continue
    end


    % Initialize Leaves object
    Leaves = LeafModelTriangle(LeafProperties.vertices, ...
                               LeafProperties.triangles);

    % Sample leaves from leaf size function
    leafScaleFactors = fun_leaf_size(overSamplingFactor*ar, ...
                                     Leaves.base_area, ...
                                     LibraryDistributions.dTypeLSD, ...
                                     [Nodes.pLSD1(iLSD1), ...
                                      Nodes.pLSD2(iLSD2)]);
    nLeaves = size(leafScaleFactors,1);
    maxLeafSize = max(max(leafScaleFactors))*max(Leaves.base_dimensions);


    % Attach the leaves to the cylinder with leaf orientation
    % distribution
    [leafDir,leafNormal,petioleStart,petioleEnd] = fun_leaf_orientation(...
        len,rad,inc,az, ...
        nLeaves, ...
        LeafProperties.petioleLengthLimits, ...
        LibraryDistributions.dTypeLODinc, ...
        LibraryDistributions.dTypeLODaz, ...
        [Nodes.pLODinc1(iLODinc1), Nodes.pLODinc2(iLODinc2)], ...
        [Nodes.pLODaz1(iLODaz1), Nodes.pLODaz2(iLODaz2)], ...
        PetioleDirectionDistribution, ...
        Phyllotaxis ...
        );

    % Add leaves to the model
    if intersectionPrevention == true
        Leaves = add_leaves(Leaves,len,rad,inc,az,ar,petioleStart, ...
                            petioleEnd,leafDir,leafNormal, ...
                            leafScaleFactors,maxLeafSize, ...
                            max(LeafProperties.petioleLengthLimits));
    else
        iLeaf = 1;
        totalArea = 0;
        while totalArea < ar && iLeaf < nLeaves
            Leaves.add_leaf(petioleEnd(iLeaf,:), ...
                leafDir(iLeaf,:), ...
                leafNormal(iLeaf,:), ...
                leafScaleFactors(iLeaf,:), ...
                1, ...
                petioleStart(iLeaf,:) ...
                );
            totalArea = totalArea ...
                        + (leafScaleFactors(iLeaf,1)^2)*Leaves.base_area;
            iLeaf = iLeaf + 1;
        end
    end

    % Add Leaves object to the node cell
    LeafObjects{iLen,iRad,iInc,iAz,iAr,iLODinc1,iLODinc2,iLODaz1,iLODaz2,...
               iLSD1,iLSD2,iCyl} = Leaves;
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

LeafCylLib.LeafObjects = LeafObjects;

% Delete waitbar
delete(wb)

%
end


