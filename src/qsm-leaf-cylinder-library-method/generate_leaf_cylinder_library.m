function LeafCylLib = generate_leaf_cylinder_library(...
                                        LibraryDistributions, ...
                                        Nodes, ...
                                        LeafProperties, ...
                                        varargin)

%% Default values

nLeafObjectsPerNode = 3;
intersectionPrevention = true;
overSamplingFactor = 10;
PetioleDirectionDistribution.flag = false;
Phyllotaxis.flag = false;
defineParallelWorkers = false;

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
                       isscalar(varargin{i+1}) && varargin{i+1} > 0 && ...
                       round(varargin{i+1}) == varargin{i+1}, ...
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
                            +" enabled.")
                end

            case 'parallelworkers'
                assert(i < NArg &&  isnumeric(varargin{i+1}) && ...
                       isscalar(varargin{i+1}) && varargin{i+1} > 0 && ...
                       round(varargin{i+1}) == varargin{i+1}, ...
                       "Argument following ''ParallelWorkers''"...
                       +" should be a positive integer.")
                nWorkers = varargin{i+1};
                defineParallelWorkers = true;

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

% Save the order of library varialbes
LeafCylLib.Properties.VariableDescription = ...
    ["(1) cylinder length", ...
     "(2) cylinder radius", ...
     "(3) cylinder inclination angle",...
     "(4) cylinder azimuth angle",...
     "(5) cylinder leaf area",...
     "(6) LOD inclination angle distribution parameter 1",...
     "(7) LOD inclination angle distribution parameter 2",...
     "(8) LOD azimuth angle distribution parameter 1",...
     "(9) LOD azimuth angle distribution parameter 2",...
     "(10) LSD distribution parameter 1",...
     "(11) LSD distribution parameter 2",...
     "(12) number of Leaves objects per library node"];

% Number of nodes per library variable
nNodesPerLibVar = [nLen,nRad,nInc,nAz,nAr,nLODinc1,nLODinc2,nLODaz1,...
                   nLODaz2,nLSD1,nLSD2,nLeafObjectsPerNode];
LeafCylLib.Properties.nNodesPerLibVar = nNodesPerLibVar;

% Total number of nodes
totalNodes = prod(nNodesPerLibVar);
LeafCylLib.Properties.totalNodes = totalNodes;

% Initializing the struct for leaf objects of the library
LeafObjects(totalNodes).Leaves = [];

% Manually define the number of parallel workers if supplied by user
if defineParallelWorkers == true
    parpool(nWorkers);
else 
    parpool
end

% Initialize a waitbar for the parfor loop
dq = parallel.pool.DataQueue;
wb = waitbar(0/totalNodes,num2str(0/totalNodes),...
             'Name','Progress bar');
wb.UserData = [0 totalNodes];
afterEach(dq,@(varargin) waitbarUpdate(wb))

% Create Leaves objects for the library
parfor iNode = 1:totalNodes

    % Make broadcast variables local to speed up parfor loop
    NodesLocal = Nodes;
    LeafPropertiesLocal = LeafProperties;
    LibraryDistributionsLocal = LibraryDistributions;

    % Find the library variable indices
    [iLen,iRad,iInc,iAz,iAr,iLODinc1,iLODinc2,iLODaz1,iLODaz2,iLSD1, ...
     iLSD2,~] = ind2sub(nNodesPerLibVar,iNode);

    % Cylinder attribute values
    len = NodesLocal.cylinderLength(iLen);
    rad = NodesLocal.cylinderRadius(iRad);
    inc = NodesLocal.cylinderInclinationAngle(iInc);
    az  = NodesLocal.cylinderAzimuthAngle(iAz);
    ar  = NodesLocal.cylinderLeafArea(iAr);

    % Initialize Leaves object
    Leaves = LeafModelTriangle(LeafPropertiesLocal.vertices, ...
                               LeafPropertiesLocal.triangles);

    % Sample leaves from leaf size function
    leafScaleFactors = fun_leaf_size(overSamplingFactor*ar, ...
                                     Leaves.base_area, ...
                                     LibraryDistributionsLocal.dTypeLSD,...
                                     [NodesLocal.pLSD1(iLSD1), ...
                                      NodesLocal.pLSD2(iLSD2)]);
    nLeaves = size(leafScaleFactors,1);
    maxLeafSize = max(max(leafScaleFactors))*max(Leaves.base_dimensions);


    % Attach the leaves to the cylinder with leaf orientation
    % distribution
    [leafDir,leafNormal,petioleStart,petioleEnd] = fun_leaf_orientation(...
        len,rad,inc,az, ...
        nLeaves, ...
        LeafPropertiesLocal.petioleLengthLimits, ...
        LibraryDistributionsLocal.dTypeLODinc, ...
        LibraryDistributionsLocal.dTypeLODaz, ...
        [NodesLocal.pLODinc1(iLODinc1), NodesLocal.pLODinc2(iLODinc2)], ...
        [NodesLocal.pLODaz1(iLODaz1), NodesLocal.pLODaz2(iLODaz2)], ...
        PetioleDirectionDistribution, ...
        Phyllotaxis ...
        );

    % Add leaves to the model
    if intersectionPrevention == true
        Leaves = add_leaves(Leaves,len,rad,inc,az,ar,petioleStart, ...
                            petioleEnd,leafDir,leafNormal, ...
                            leafScaleFactors,maxLeafSize, ...
                            max(LeafPropertiesLocal.petioleLengthLimits));
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

    % Add Leaves object to the library struct
    LeafObjects(iNode).Leaves = Leaves;
    Leaves = [];

    % Update waitbar
    send(dq,iNode);

end

% Close waitbar
close(wb);

% Closing the parallel pool
delete(gcp('nocreate'));

% Assigning the leaf object struct to the leaf cylinder library
LeafCylLib.LeafObjects = LeafObjects;

% Function for waitbar update of parfor loop
    function waitbarUpdate(wb)
        iters = wb.UserData;
        iters(1) = iters(1) + 1;
        waitbar(iters(1)/iters(2),wb, ...
                sprintf('Leaf cylinder generation in progress\n%d/%d', ...
                iters(1),iters(2)));
        wb.UserData = iters;
    end

end


