function LeafCylLib = generate_leaf_cylinder_library(Nodes, ...
                                                     LeafDistributions, ...
                                                     twigLengthLimits, ...
                                                     vertices, ...
                                                     tris, ...
                                                     varargin)

%% Default values

nLeafObjectsPerNode = 3;
intersectionPrevention = false;
overSamplingFactor = 2;
TwigDirectionDistribution.flag = false;
Phyllotaxis.flag = false;

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
LeafCylLib.LeafDistributions = LeafDistributions;

% Library nodes
LeafCylLib.Nodes = Nodes;

% Twig length limits
LeafCylLib.twigLengthLimits = twigLengthLimits;

% Leaf base information
LeafCylLib.LeafBaseModel.vertices = vertices;
LeafCylLib.LeafBaseModel.tris = tris;

% Number of leaf objects generated for each node
LeafCylLib.nLeafObjectsPerNode = nLeafObjectsPerNode;

%% Leaf-cylinder generation

% Ranges of library variables
nLodInc1 = length(Nodes.LodInc1);
nLodInc2 = length(Nodes.LodInc2);
nLodAz1  = length(Nodes.LodAz1);
nLodAz2  = length(Nodes.LodAz2);
nLsd1    = length(Nodes.Lsd1);
nLsd2    = length(Nodes.Lsd2);

nLen = length(Nodes.cylinderLength);
nRad = length(Nodes.cylinderRadius);
nInc = length(Nodes.cylinderInclinationAngle);
nAz  = length(Nodes.cylinderAzimuthAngle);
nAr  = length(Nodes.cylinderLeafArea);

% Number of nodes per library variable
nNodesPerLibVar = [nLen,nRad,nInc,nAz,nAr,nLodInc1,nLodInc2,nLodAz1,...
                   nLodAz2,nLsd1,nLsd2,nLeafObjectsPerNode];

% Initializing the cell for leaf-cylinder library
LeavesObjects = cell(nNodesPerLibVar);

% Total number of iterations
nIter = prod(nNodesPerLibVar);
iIter = 0;
wb = waitbar(iIter/nIter,num2str(iIter/nIter),...
             'Name','Progress bar');

% Create Leaves objects for the library
for  iLodInc1 = 1:nLodInc1
    for iLodInc2 = 1:nLodInc2
    for iLodAz1  = 1:nLodAz1
    for iLodAz2  = 1:nLodAz2
    for iLsd1    = 1:nLsd1
    for iLsd2    = 1:nLsd2

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
        Leaves = LeavesObjects{iLen,iRad,iInc,1,iAr,iLodInc1,iLodInc2, ...
                               iLodAz1,iLodAz2,iLsd1,iLsd2,iCyl};
        LeavesObjects{iLen,iRad,iInc,iAz,iAr,iLodInc1,iLodInc2,iLodAz1, ...
                      iLodAz2,iLsd1,iLsd2,iCyl} = Leaves;
        clearvars Leaves;
        % Skip to next iteration of the active for-loop
        continue
    end


    % Initialize Leaves object
    Leaves = LeafModelTriangle(vertices,tris);

    % Sample leaves from leaf size function
    [leafScaleFactors,nLeaves,maxLeafSize] = fun_leaf_size( ...
        overSamplingFactor*ar, ...
        Leaves.base_area, ...
        LeafDistributions.dTypeLsd, ...
        [Nodes.Lsd1(iLsd1), Nodes.Lsd2(iLsd2)] ...
        );

    % Attach the leaves to the cylinder with leaf orientation
    % distribution
    [leafDir,leafNormal,twigStart,twigEnd] = fun_leaf_orientation( ...
        len,rad,inc,az, ...
        nLeaves, ...
        twigLengthLimits, ...
        LeafDistributions.dTypeLodInc, ...
        LeafDistributions.dTypeLodAz, ...
        [Nodes.LodInc1(iLodInc1), Nodes.LodInc2(iLodInc2)], ...
        [Nodes.LodAz1(iLodAz1), Nodes.LodAz2(iLodAz2)], ...
        TwigDirectionDistribution, ...
        Phyllotaxis ...
        );

    % Add leaves to the model
    if intersectionPrevention == true
        Leaves = add_leaves(Leaves,len,rad,inc,az,ar,twigStart,twigEnd, ...
                            leafDir,leafNormal,leafScaleFactors, ...
                            maxLeafSize,max(twigLengthLimits));
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
    LeavesObjects{iLen,iRad,iInc,iAz,iAr,iLodInc1,iLodInc2,iLodAz1, ...
                  iLodAz2,iLsd1,iLsd2,iCyl} = Leaves;
    clearvars Leaves;

    end % cylinderLeafArea
    end % cylinderAzimuthAngle
    end % cylinderInclinationAngle
    end % cylinderRadius
    end % cylinderLenght

    end % nLeafCylindersPerNode

    end % Lsd2
    end % Lsd1
    end % LodAz2
    end % LodAz1
    end % LodInc2
end     % LodInc1

LeafCylLib.LeavesObjects = LeavesObjects;

% Delete waitbar
delete(wb)

%
end


