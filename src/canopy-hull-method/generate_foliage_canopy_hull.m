function [Leaves,aShape] = generate_foliage_canopy_hull(treePointCloud, ...
                                                   TargetDistributions, ...
                                                   LeafProperties, ...
                                                   totalLeafArea, ...
                                                   varargin ...
                                                   )
%% Check the correctness of inputs
check_inputs_canopy_hull(TargetDistributions, ...
                         LeafProperties, ...
                         totalLeafArea);

%% Initialize values
alpha = 1;
stemCoordinates = [0 0 0; 0 0 max(treePointCloud(:,3))];
pcSamplingWeights = 0;
voxelEdge = 0.15;
voxelThreshold = 5;
leafAttLabels = true(size(treePointCloud,1),1);
intersectionPrevention = true;
fDist_h = []; fDist_d = []; fDist_c = [];
maxfDist_h = []; maxfDist_d = []; maxfDist_c = [];

%% Read inputs

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'alpha'
                alpha = varargin{i+1};
            case 'stemcoordinates'
                stemCoordinates = varargin{i+1};
            case 'pcsamplingweights'
                pcSamplingWeights = varargin{i+1};
            case 'voxelsize'
                voxelEdge = varargin{i+1};
            case 'voxelthreshold'
                voxelThreshold = varargin{i+1};
            case 'leafattractorlabels'
                leafAttLabels = logical(varargin{i+1});
            case 'intersectionprevention'
                intersectionPrevention = varargin{i+1};
        end
    end
    i = i + 1;
end
%% Function definitions
fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_weibull = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));
fun_dewit = @(x,a,b) (1 + a*cos(b*x))/(pi/2+(a/b)*sin(b*pi/2));
%% Leaf area density distribution functions
% Distribution function and parameters for relative height
dTypeLADDh = TargetDistributions.dTypeLADDh;
pLADDh    = TargetDistributions.pLADDh;
switch dTypeLADDh
    case 'uniform'
        fDist_h = @(h) 1*ones(size(h));
        maxfDist_h = 1;
    case 'polynomial'
        fDist_h = @(h) polyval(pLADDh,h);
        maxfDist_h = max([0,polynomial_upper_limit(pLADDh)]);
    case 'polynomialmixture'
        nP = (length(pLADDh)-1)/2; % number of polynomial coefficients
        p1 = pLADDh(1:nP); % coefficients of the first polynomial
        p2 = pLADDh((nP+1):(2*nP)); % coefficients of the second polynom.
        w = pLADDh(end); % mixture model weight
        fDist_h = @(h) w*polyval(p1,h) + (1-w)*polyval(p2,h);
        maxfDist_h = max([0,polynomial_upper_limit(w*p1+(1-w)*p2)]);
    case 'weibull'
        l = pLADDh(1); % scale parameter
        k = pLADDh(2); % shape parameter
        fDist_h = @(h) fun_weibull(h,l,k);
        maxfDist_h = weibull_upper_limit(l,k);
    case 'weibullmixture'
        l1 = pLADDh(1); k1 = pLADDh(2); % parameters of the first dist.
        l2 = pLADDh(3); k2 = pLADDh(4); % parameters of the second dist.
        w = pLADDh(5); % mixture model weight
        fDist_h = @(h) w*fun_weibull(h,l1,k1) + (1-w)*fun_weibull(h,l2,k2);
        maxfDist_h = weibull_mm_upper_limit(pLADDh);
    case 'beta'
        a = pLADDh(1);
        b = pLADDh(2);
        fDist_h = @(h) fun_beta(h,a,b);
        maxfDist_h = beta_upper_limit(a,b);
    case 'betamixture'
        a1 = pLADDh(1); b1 = pLADDh(2); % parameters of the first dist.
        a2 = pLADDh(3); b2 = pLADDh(4); % parameters of the second dist.
        w = pLADDh(5); % mixture model weight
        fDist_h = @(h) w*fun_beta(h,a1,b1) + (1-w)*fun_beta(h,a2,b2);
        maxfDist_h = beta_mm_upper_limit(pLADDh);
end
% Distribution function and parameters for relative distance along
% sub-branch
dTypeLADDd = TargetDistributions.dTypeLADDd;
pLADDd    = TargetDistributions.pLADDd;
switch dTypeLADDd
    case 'uniform'
        fDist_d = @(d) 1*ones(size(d));
        maxfDist_d = 1;
    case 'polynomial'
        fDist_d = @(d) polyval(pLADDd,d);
        maxfDist_d = max([0,polynomial_upper_limit(pLADDd)]);
    case 'polynomialmixture'
        nP = (length(pLADDd)-1)/2; % order of polynomial
        p1 = pLADDd(1:nP); % coefficients of the first polynomial
        p2 = pLADDd((nP+1):(2*nP)); % coefficients of the second polynom.
        w = pLADDd(end); % mixture model weight
        fDist_d = @(d) w*polyval(p1,d) + (1-w)*polyval(p2,d);
        maxfDist_d = max([0,polynomial_upper_limit(w*p1+(1-w)*p2)]);
    case 'weibull'
        l = pLADDd(1); % scale parameter
        k = pLADDd(2); % shape parameter
        fDist_d = @(d) fun_weibull(d,l,k);
        maxfDist_d = weibull_upper_limit(l,k);
    case 'weibullmixture'
        l1 = pLADDd(1); k1 = pLADDd(2); % parameters of the first dist.
        l2 = pLADDd(3); k2 = pLADDd(4); % parameters of the second dist.
        w = pLADDd(5); % mixture model weight
        fDist_d = @(d) w*fun_weibull(d,l1,k1) + (1-w)*fun_weibull(d,l2,k2);
        maxfDist_d = weibull_mm_upper_limit(pLADDd);
    case 'beta'
        a = pLADDd(1);
        b = pLADDd(2);
        fDist_d = @(d) fun_beta(d,a,b);
        maxfDist_d = beta_upper_limit(a,b);
    case 'betamixture'
        a1 = pLADDd(1); b1 = pLADDd(2); % parameters of the first dist.
        a2 = pLADDd(3); b2 = pLADDd(4); % parameters of the second dist.
        w = pLADDd(5); % mixture model weight
        fDist_d = @(d) w*fun_beta(d,a1,b1) + (1-w)*fun_beta(d,a2,b2);
        maxfDist_d = beta_mm_upper_limit(pLADDd);
end
% Distribution function and parameters for compass direction
dTypeLADDc = TargetDistributions.dTypeLADDc;
pLADDc    = TargetDistributions.pLADDc;
switch dTypeLADDc
    case 'uniform'
        fDist_c = @(c) 1/(2*pi)*ones(size(c));
        maxfDist_c = 1/(2*pi);
    case 'vonmises'
        m = pLADDc(1); % mean
        k = pLADDc(2); % measure of concentration
        fDist_c = @(c) fun_vonmises(c,m,k);
        maxfDist_c = fDist_c(m);
    case 'vonmisesmixture'
        m1 = pLADDc(1); k1 = pLADDc(2); % parameters of the first dist.
        m2 = pLADDc(3); k2 = pLADDc(4); % parameters of the second dist.
        w = pLADDc(5); % mixture model weight
        fDist_c = @(c) w*fun_vonmises(c,m1,k1) ...
                       + (1-w)*fun_vonmises(c,m2,k2);
        maxfDist_c = vonmises_mm_upper_limit(pLADDc);
end

%% Leaf orientation distribution functions
% Distribution funtion and parameter value function for leaf inclination
% angle distribution
dTypeLODinc    = TargetDistributions.dTypeLODinc;
fun_inc_params = TargetDistributions.fun_inc_params;
switch dTypeLODinc
    case 'uniform'
        % Uniform distribution
        f_inc = @(x,p) 2/pi*ones(size(x));
        max_f_inc = @(p) 2/pi;
        % Set sampling type to rejection sampling
        lod_inc_sampling = 'rejection sampling';
    case 'spherical'
        % Spherical distribution function
        f_inc = @(x,p) sin(x);
        % Inverse of cumulative density function
        F_inc_inv = @(y,p) acos(1-y);
        % Set sampling type to inverse sampling
        lod_inc_sampling = 'inverse sampling';
    case 'dewit'
        % Generalized de Wit's distribution function
        f_inc = @(x,p) fun_dewit(x,p(1),p(2));
        % Upper limit for the distribution
        max_f_inc = @(p) max(f_inc([0 1 2 3]*(pi/p(2)),p));
        % Set sampling type to rejection sampling
        lod_inc_sampling = 'rejection sampling';
    case 'beta'
        % Beta distribution density function
        f_inc = @(x,p) fun_beta(2*x/pi,p(1),p(2));
        % Inverse of cumulative density function
        F_inc_inv = @(y,p) (pi/2)*betaincinv(y,p);
        % Set sampling type to inverse sampling
        lod_inc_sampling = 'inverse sampling';
    case 'constant'
        % Set sampling type to constant
        lod_inc_sampling = 'constant';
end
% Distribution function and parameter value function for leaf azimuth angle
% distribution
dTypeLODaz    = TargetDistributions.dTypeLODaz;
fun_az_params = TargetDistributions.fun_az_params;
switch dTypeLODaz
    case 'uniform'
        % Uniform distribution
        f_az = @(x,p) 1/(2*pi)*ones(size(x));
        max_f_inc = @(p) 1/(2*pi);
        % Set sampling type to rejection sampling
        lod_az_sampling = 'rejection sampling';
    case 'vonmises'
        % Von Mises distribution density function
        f_az = @(x,p) fun_vonmises(x,p(1),p(2));
        % Upper limit for the distribution
        max_f_az = @(p) fun_vonmises(p(1),p(1),p(2));
        % Set sampling type to rejection sampling
        lod_az_sampling = 'rejection sampling';
    case 'constant'
        % Set sampling type to constant
        lod_az_sampling = 'constant';
end
%% Leaf size distriubtion functions
dTypeLSD        = TargetDistributions.dTypeLSD;
fun_size_params = TargetDistributions.fun_size_params;
%% Generate alpha shape on the point cloud
disp('---------------------------------------')
disp('Generating alphashape on point cloud')
tic
aShape = alphaShape(treePointCloud,alpha);
toc
%% Extreme points of point cloud
% Maximum horizontal distance from origin to point cloud member
horzDistances = sqrt(treePointCloud(:,1).^2 + treePointCloud(:,2).^2);
maxHorzDist = max(horzDistances);

% Maximum height of the point cloud members
maxHeight = max(treePointCloud(:,3));

% Axiswise extreme points
xMin = min(treePointCloud(:,1));
xMax = max(treePointCloud(:,1));
yMin = min(treePointCloud(:,2));
yMax = max(treePointCloud(:,2));
zMin = min(treePointCloud(:,3));
zMax = max(treePointCloud(:,3));

%% Caculating voxel array and finding which voxels contain points
disp('Calculating point cloud voxelization')

tic

% Voxel parameters
nx = ceil((xMax-xMin)/voxelEdge);
ny = ceil((yMax-yMin)/voxelEdge);
nz = ceil((zMax-zMin)/voxelEdge);
xEdges = [xMin xMin+cumsum(voxelEdge*ones(1,nx))];
yEdges = [yMin yMin+cumsum(voxelEdge*ones(1,ny))];
zEdges = [zMin zMin+cumsum(voxelEdge*ones(1,nz))];
xCenters = 0.5*(xEdges(1:end-1)+xEdges(2:end));
yCenters = 0.5*(yEdges(1:end-1)+yEdges(2:end));
zCenters = 0.5*(zEdges(1:end-1)+zEdges(2:end));

% Initialize visualizaiton of point cloud voxelization
if any(pcSamplingWeights) == true
    figure, clf, hold on, grid on, axis equal, view(3)
    plot3(treePointCloud(:,1),treePointCloud(:,2),treePointCloud(:,3), ...
        'g.','MarkerSize',0.5)
    protoVertices = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1; 
                     1 0 1];
    voxelFaces = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
end

% Pick the leaf attractor points from point cloud (if not supplied by user, 
% all point cloud points act as leaf attractors)
pcRemaining = treePointCloud(leafAttLabels,:);

% Initialize loop variables
voxelInd = cell(nx,ny,nz);
voxelCenDir = zeros(nx,ny,nz);
pcVoxels = zeros(nx,ny,nz);
for ix = 1:nx
    % Slice of point cloud corrseponding to the x voxel coordinates
    inXslice = all([(pcRemaining(:,1) < xEdges(ix+1)) ...
                    (pcRemaining(:,1) >= xEdges(ix))],2);
    pcXslice = pcRemaining(inXslice,:);
    for iy = 1:ny
        % Column of point cloud corresponding to the x and y voxel
        % coordinates
        inXYcolumn = all([(pcXslice(:,2) < yEdges(iy+1)) ...
                          (pcXslice(:,2) >= yEdges(iy))],2);
        pcXYcolumn = pcXslice(inXYcolumn,:);
        for iz = 1:nz
            % Cell for voxel indices
            voxelInd{ix,iy,iz} = [ix iy iz];
            % Compass direction of voxel center point wrt. stem
            voxelCenter = [xCenters(ix) yCenters(iy) zCenters(iz)];
            voxelCenDir(ix,iy,iz) = xyz_to_compass_dir(voxelCenter, ...
                                                       stemCoordinates);
            % Find the number of points inside the voxel
            inXYZvoxel = all([(pcXYcolumn(:,3) < zEdges(iz+1)) ...
                              (pcXYcolumn(:,3) >= zEdges(iz))],2);
            % If the number of points surpasses voxel point threshold, set
            % voxel value to 1, otherwise 0
            if sum(inXYZvoxel) >= voxelThreshold
                pcVoxels(ix,iy,iz) = 1;
                if any(pcSamplingWeights) == true
                    % Visualize the voxel
                    voxelVertices = voxelEdge*protoVertices ...
                                    + [xEdges(ix) yEdges(iy) zEdges(iz)];
                    patch('vertices',voxelVertices,'faces',voxelFaces, ...
                          'facecolor','r','facealpha',0.0);
                end
            end
        end
    end
end
% Set point cloud sampling and voxelization information into a struct
PCSampling.binEdges = linspace(0,1,length(pcSamplingWeights)+1);
PCSampling.weights  = pcSamplingWeights;
PCSampling.xEdges = xEdges;
PCSampling.yEdges = yEdges;
PCSampling.zEdges = zEdges;
PCSampling.voxels = pcVoxels;
PCSampling.voxelIndex = voxelInd;
PCSampling.voxelCenterDirection = voxelCenDir;
PCSampling.nPCSampled = zeros(size(pcSamplingWeights));
PCSampling.nTotalSampled = zeros(size(pcSamplingWeights));
PCSampling.ratio = zeros(size(pcSamplingWeights));
if sum(pcSamplingWeights) > 0
    PCSampling.flag = true;
else
    PCSampling.flag = false;
end

toc 

%% Initialize leaf object and candidate leaf area
Leaves = LeafModelTriangle(LeafProperties.vertices, ...
                           LeafProperties.triangles);
baseArea = Leaves.base_area;

if intersectionPrevention == true
    candidateArea = 2*totalLeafArea;
else
    candidateArea = totalLeafArea;
end

%% Initializing leaf attribute variables
nInit = 10000;
leafStartPoints  = zeros(nInit,3);
incAngles        = zeros(nInit,1);
azAngles         = zeros(nInit,1);
leafNormal       = zeros(nInit,3);
leafDir          = zeros(nInit,3);
leafScaleFactors = zeros(nInit,3);

%% Sample leaf positions, orientations and sizes
disp('Sampling leaf positions')
tic
iLeaf = 0;
leafArea = 0;
iVarExt = 0;
while leafArea < candidateArea
    % Increase leaf index
    iLeaf = iLeaf + 1;
    % Sampling leaf position   
    [leafSP,PCSampling] = sample_leaf_position(aShape, ...
                              fDist_h,fDist_d,fDist_c, ...
                              maxfDist_h,maxfDist_d,maxfDist_c, ...
                              xMin,xMax,yMin,yMax,maxHeight, ...
                              maxHorzDist,stemCoordinates,PCSampling);
    leafStartPoints(iLeaf,:) = leafSP;

    % Leaf height value
    hLeaf = leafStartPoints(iLeaf,3)/maxHeight;
    % Leaf compass direction value
    cLeaf = xyz_to_compass_dir(leafStartPoints(iLeaf,:),stemCoordinates);
    % Leaf distance from stem value
    dLeaf = stem_to_alphashape_edge(aShape,hLeaf,cLeaf,maxHorzDist, ...
                                    maxHeight,stemCoordinates);

    % Sampling leaf normal orientation with LOD
    % Leaf inclination angle
    switch lod_inc_sampling
        case 'rejection sampling'
            % Sample inclination angle value with acceptance-rejection
            % sampling
            accepted = 0;
            while accepted == 0
                incProposal = rand(1)*pi/2;
                pars = fun_inc_params(hLeaf,dLeaf,cLeaf);
                funValue = f_inc(incProposal,pars);
                vertValue = rand(1)*max_f_inc(pars);
                if vertValue < funValue
                    incAngles(iLeaf) = incProposal;
                    accepted = 1;
                end
            end
        case 'inverse sampling'
            % Sample inclination angle by inverse transform sampling
            u = rand(1);
            pars = fun_inc_params(hLeaf,dLeaf,cLeaf);
            incAngles(iLeaf) = F_inc_inv(u,pars);
        case 'constant'
            pars = fun_inc_params(hLeaf,dLeaf,cLeaf);
            incAngles(iLeaf) = pars;
    end
    % Leaf azimuth angle
    switch lod_az_sampling
        case 'rejection sampling'
            % Sample azimuth angle value with acceptance-rejection
            % sampling
            accepted = 0;
            while accepted == 0
                azProposal = rand(1)*2*pi;
                pars = fun_az_params(hLeaf,dLeaf,cLeaf);
                funValue = f_az(azProposal,pars);
                vertValue = rand(1)*max_f_az(pars);
                if vertValue < funValue
                    azAngles(iLeaf) = azProposal;
                    accepted = 1;
                end
            end
        case 'constant'
            pars = fun_az_params(hLeaf,dLeaf,cLeaf);
            azAngles(iLeaf) = pars;
    end

    % Sampling leaf surface area with LSD
    switch dTypeLSD
        case 'uniform'
            pars = fun_size_params(hLeaf,dLeaf,cLeaf);
            sampledArea = (pars(2)-pars(1))*rand(1) + pars(1);
        case 'normal'
            pars = fun_size_params(hLeaf,dLeaf,cLeaf);
            sampledArea = sqrt(pars(2))*randn(1) + pars(1);
        case 'constant'
            pars = fun_size_params(hLeaf,dLeaf,cLeaf);
            sampledArea = pars;
    end
    leafArea = leafArea + sampledArea;
    % Store leaf scaling factor (same for all dimensions)
    leafScaleFactors(iLeaf,:) = sqrt(sampledArea/baseArea)*ones(1,3);

    % Unit vector of the leaf normal (y-axis assumed as north direction)
    leafNormal(iLeaf,:) = [sin(incAngles(iLeaf))*cos(azAngles(iLeaf)+pi/2),...
                           sin(incAngles(iLeaf))*sin(azAngles(iLeaf)+pi/2),...
                           cos(incAngles(iLeaf))];


    % Sampling leaf direction uniformly
    if abs(dot(leafNormal(iLeaf,:),[0 0 1])) > 0.999
        dirVec = [0 1 0];
    else
        dirVec = cross(leafNormal(iLeaf,:),[0 0 1]);
        dirVec = dirVec/norm(dirVec);
    end
    leafDir(iLeaf,:) = (rotation_matrix(leafNormal(iLeaf,:),2*pi*rand(1)) ...
                       *dirVec')';

    % Extend preallocated variables if needed
    if (iLeaf+1) > (iVarExt+1)*nInit
        leafStartPoints  = [leafStartPoints;  zeros(nInit,3)]; %#ok<AGROW>
        incAngles        = [incAngles;        zeros(nInit,1)]; %#ok<AGROW>
        azAngles         = [azAngles;         zeros(nInit,1)]; %#ok<AGROW>
        leafNormal       = [leafNormal;       zeros(nInit,3)]; %#ok<AGROW>
        leafDir          = [leafDir;          zeros(nInit,3)]; %#ok<AGROW>
        leafScaleFactors = [leafScaleFactors; zeros(nInit,3)]; %#ok<AGROW>
        iVarExt = iVarExt + 1;
    end
end

% Trim unnecessary rows
leafStartPoints = leafStartPoints(1:iLeaf,:);
leafNormal = leafNormal(1:iLeaf,:);
leafDir = leafDir(1:iLeaf,:);
leafScaleFactors = leafScaleFactors(1:iLeaf,:);
toc

%% Add leaves to the shape without intersections
if intersectionPrevention == true
    disp('Adding leaves to the model without intersections')
    tic
    Leaves = add_leaves_to_alphashape(aShape, ...
                                      Leaves, ...
                                      totalLeafArea, ....
                                      leafStartPoints, ...
                                      leafNormal, ...
                                      leafDir, ...
                                      leafScaleFactors);
    toc
else
    disp('Adding leaves to the model')
    tic
    iLeaf = 1;
    areaAdded = 0;
    while areaAdded < totalLeafArea && iLeaf < size(leafStartPoints,1)
        leafParent = NaN;
        petioleStart  = NaN;
        Leaves.add_leaf(leafStartPoints(iLeaf,:), ...
                        leafDir(iLeaf,:), ...
                        leafNormal(iLeaf,:), ...
                        leafScaleFactors(iLeaf,:), ...
                        leafParent, ...
                        petioleStart);
        areaAdded = areaAdded ...
                    + (leafScaleFactors(iLeaf,1)^2)*Leaves.base_area;
        iLeaf = iLeaf + 1;
    end
    toc
end
disp('Foliage generation finished')
