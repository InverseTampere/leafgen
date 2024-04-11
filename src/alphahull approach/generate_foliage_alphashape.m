function [Leaves,shp] = generate_foliage_alphashape(treePointCloud, ...
                                                    TargetDistributions, ...
                                                    totalLeafArea, ...
                                                    vertices, ...
                                                    tris, ...
                                                    varargin ...
                                                    )
%% Initialize values
alpha = 1;
flagStemCoordinates = false;
nBinsPC_h = 100;
nBinsPC_d = 100;
nBinsPC_c = 100;

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
                flagStemCoordinates = true;
            case 'pcpositionsampling'
                pcSamplingWeight = varargin{i+1};
                flagPCPositionSampling = true;
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
dType_h = TargetDistributions.dType_h;
p_h     = TargetDistributions.p_h;
% nBins_h = TargetDistributions.nBins_h;
switch dType_h
    case 'uniform'
        fDist_h = @(h) 1;
        maxfDist_h = 1;
    case 'polynomial'
        fDist_h = @(h) polyval(p_h,h);
        maxfDist_h = max([0,polynomial_upper_limit(p_h)]);
    case 'polynomialmixturemodel'
        nP = (length(p_h)-1)/2; % number of polynomial coefficients
        p1 = p_h(1:nP); % coefficients of the first polynomial
        p2 = p_h((nP+1):(2*nP)); % coefficients of the second polynomial
        w = p_h(end); % mixture model weight
        fDist_h = @(h) w*polyval(p1,h) + (1-w)*polyval(p2,h);
        maxfDist_h = max([0,polynomial_upper_limit(w*p1+(1-w)*p2)]);
    case 'weibull'
        l = p_h(1); % scale parameter
        k = p_h(2); % shape parameter
        fDist_h = @(h) fun_weibull(h,l,k);
        maxfDist_h = weibull_upper_limit(l,k);
    case 'weibullmixturemodel'
        l1 = p_h(1); k1 = p_h(2); % parameters of the first distribution
        l2 = p_h(3); k2 = p_h(4); % parameters of the second distribution
        w = p_h(5); % mixture model weight
        fDist_h = @(h) w*fun_weibull(h,l1,k1) + (1-w)*fun_weibull(h,l2,k2);
        maxfDist_h = weibull_mm_upper_limit(p_h);
    case 'beta'
        a = p_h(1);
        b = p_h(2);
        fDist_h = @(h) fun_beta(h,a,b);
        maxfDist_h = beta_upper_limit(a,b);
    case 'betamixturemodel'
        a1 = p_h(1); b1 = p_h(2); % parameters of the first distribution
        a2 = p_h(3); b2 = p_h(4); % parameters of the second distribution
        w = p_h(5); % mixture model weight
        fDist_h = @(h) w*fun_beta(h,a1,b1) + (1-w)*fun_beta(h,a2,b2);
        maxfDist_h = beta_mm_upper_limit(p_h);
end
% Distribution function and parameters for relative distance along
% sub-branch
dType_d = TargetDistributions.dType_d;
p_d     = TargetDistributions.p_d;
% nBins_d = TargetDistributions.nBins_d;
switch dType_d
    case 'uniform'
        fDist_d = @(d) 1;
        maxfDist_d = 1;
    case 'polynomial'
        fDist_d = @(d) polyval(p_d,d);
        maxfDist_d = max([0,polynomial_upper_limit(p_d)]);
    case 'polynomialmixturemodel'
        nP = (length(p_d)-1)/2; % order of polynomial
        p1 = p_d(1:nP); % coefficients of the first polynomial
        p2 = p_d((nP+1):(2*nP)); % coefficients of the second polynomial
        w = p_d(end); % mixture model weight
        fDist_d = @(d) w*polyval(p1,d) + (1-w)*polyval(p2,d);
        maxfDist_d = max([0,polynomial_upper_limit(w*p1+(1-w)*p2)]);
    case 'weibull'
        l = p_d(1); % scale parameter
        k = p_d(2); % shape parameter
        fDist_d = @(d) fun_weibull(d,l,k);
        maxfDist_d = weibull_upper_limit(l,k);
    case 'weibullmixturemodel'
        l1 = p_d(1); k1 = p_d(2); % parameters of the first distribution
        l2 = p_d(3); k2 = p_d(4); % parameters of the second distribution
        w = p_d(5); % mixture model weight
        fDist_d = @(d) w*fun_weibull(d,l1,k1) + (1-w)*fun_weibull(d,l2,k2);
        maxfDist_d = weibull_mm_upper_limit(p_d);
    case 'beta'
        a = p_d(1);
        b = p_d(2);
        fDist_d = @(d) fun_beta(d,a,b);
        maxfDist_d = beta_upper_limit(a,b);
    case 'betamixturemodel'
        a1 = p_d(1); b1 = p_d(2); % parameters of the first distribution
        a2 = p_d(3); b2 = p_d(4); % parameters of the second distribution
        w = p_d(5); % mixture model weight
        fDist_d = @(d) w*fun_beta(d,a1,b1) + (1-w)*fun_beta(d,a2,b2);
        maxfDist_d = beta_mm_upper_limit(p_d);
end
% Distribution function and parameters for compass direction
dType_c = TargetDistributions.dType_c;
p_c     = TargetDistributions.p_c;
% nBins_c = TargetDistributions.nBins_c;
switch dType_c
    case 'uniform'
        fDist_c = @(c) 1;
        maxfDist_c = 1;
    case 'vonmises'
        m = p_c(1); % mean
        k = p_c(2); % measure of concentration
        fDist_c = @(c) fun_vonmises(c,m,k);
        maxfDist_c = fDist_c(m);
    case 'vonmisesmixturemodel'
        m1 = p_c(1); k1 = p_c(2); % parameters of the first distribution
        m2 = p_c(3); k2 = p_c(4); % parameters of the second distribution
        w = p_c(5); % mixture model weight
        fDist_c = @(c) w*fun_vonmises(c,m1,k1) ...
                       + (1-w)*fun_vonmises(c,m2,k2);
        maxfDist_c = vonmises_mm_upper_limit(p_c);
end

%% Leaf orientation distribution functions
% Distribution funtion and parameter value function for leaf inclination
% angle distribution
dType_inc = TargetDistributions.dType_inc;
fun_inc_params = TargetDistributions.fun_inc_params;
switch dType_inc
    case 'dewit'
        % de Wit distribution function
        f_inc = @(x,a,b) (1 + a*cos(b*x))/(pi/2+(a/b)*sin(b*pi/2));
    case 'beta'
        % Beta distribution density function
        f_inc = @(x,a,b) (1/beta(a,b))*(x/(pi/2)).^(a-1) ...
                         .*(1-(x/(pi/2))).^(b-1);
        % Inverse of cumulative density function
        F_inc_inv = @(y,a,b) (pi/2)*betaincinv(y,a,b);
end
% Distribution function and parameter value function for leaf azimuth angle
% distribution
dType_az = TargetDistributions.dType_az;
fun_az_params = TargetDistributions.fun_az_params;
switch dType_az
    case 'vonmises'
        % Von Mises distribution density function
        f_az = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));
end
%% Leaf size distriubtion functions
dType_size = TargetDistributions.dType_size;
fun_size_params = TargetDistributions.fun_size_params;
%% Generate alpha shape on the point cloud
disp('---------------------------------------')
disp('Generating alphashape on point cloud')
tic
shp = alphaShape(treePointCloud,alpha);
toc
%% Extreme points of point cloud
% Maximum horizontal distance from origin to point cloud member
horzDistances = sqrt(treePointCloud(:,1).^2 + treePointCloud(:,2).^2);
maxHorzDist = max(horzDistances);

% Maximum height of the point cloud members
maxHeight = max(treePointCloud(:,3));

%% Point cloud probability voxelization
disp('Calculating point cloud voxel densities')
tic
if flagPCPositionSampling == true
    
    % Bin edges for each cylindrical variable
    binEdgesH = linspace(0,maxHeight,nBinsPC_h+1);
    binEdgesD = log(linspace(exp(0),exp(maxHorzDist),nBinsPC_d+1));
    binEdgesC = linspace(0,2*pi,nBinsPC_c+1);

    pcDensityInVolume = zeros(nBinsPC_h*nBinsPC_d*nBinsPC_c,1); %zeros(nBinsPC_h,nBinsPC_d,nBinsPC_c);
    volumeInd = zeros(nBinsPC_h*nBinsPC_d*nBinsPC_c,3);
    k = 0;
    for i_h = 1:nBinsPC_h
        % Part of point cloud within the height bin
        inHeightBin = all([(treePointCloud(:,3)<binEdgesH(i_h+1)) ...
                           (treePointCloud(:,3)>=binEdgesH(i_h))],2);
        pcHBin = treePointCloud(inHeightBin,:);

        for i_d = 1:nBinsPC_d
            % Part of point cloud within the height and distance bin
            pointDist = sqrt(pcHBin(:,1).^2+pcHBin(:,2).^2);
            inDistanceBin = all([(pointDist<binEdgesD(i_d+1)) ...
                                 (pointDist>=binEdgesD(i_d))],2);
            pcHDBin = pcHBin(inDistanceBin,:);

            for i_c = 1:nBinsPC_c
                % Increment loop counter
                k = k + 1;
                % Volume of cylinder slice
                hSlice = binEdgesH(i_h+1) - binEdgesH(i_h);
                cRatio = (binEdgesC(i_c+1)-binEdgesC(i_c))/(2*pi);
                if i_d == 1
                    vol = cRatio*hSlice*(pi*binEdgesD(i_d+1).^2);
                else
                    vol = cRatio*hSlice*(pi*binEdgesD(i_d+1).^2 ...
                          - pi*binEdgesD(i_d).^2);
                end
                % Volume indices
                volumeInd(k,:) = [i_h i_d i_c];
                % Find all point cloud elements within the slice
                hPoints = pcHDBin(:,3);
                dPoints = sqrt(pcHDBin(:,1).^2+pcHDBin(:,2).^2);
                cPoints = zeros(size(pcHDBin,1),1);
                xNeg = pcHDBin(:,1) <= 0;
                xPos = pcHDBin(:,1) >  0;
                xNegDir = [pcHDBin(xNeg,1:2) zeros(sum(xNeg),1)];
                xNegDir = xNegDir./sqrt(sum(xNegDir.^2,2));
                xPosDir = [pcHDBin(xPos,1:2) zeros(sum(xPos),1)];
                xPosDir = xPosDir./sqrt(sum(xPosDir.^2,2));
                cPoints(xNeg) = acos(xNegDir*[0 1 0]');
                cPoints(xPos) = 2*pi - acos(xPosDir*[0 1 0]');
                pointsInSlice = all([hPoints<binEdgesH(i_h+1) ...
                                     hPoints>=binEdgesH(i_h) ...
                                     dPoints<binEdgesD(i_d+1) ...
                                     dPoints>=binEdgesD(i_d) ...
                                     cPoints<binEdgesC(i_c+1) ...
                                     cPoints>=binEdgesC(i_c)] ...
                                     ,2);
                % Calculate the point density in slice
                pcDensityInVolume(k) = sum(pointsInSlice)/vol;
            end
        end
    end
    cumulativePCD = cumsum(pcDensityInVolume);
    cumulativeDistPCD = cumulativePCD./cumulativePCD(end);
    
end
toc

%% Initialize leaf object and candidate leaf area
Leaves = LeafModelTriangle(vertices,tris);
baseArea = Leaves.base_area;
candidateArea = 2*totalLeafArea;

%% Initializing leaf attribute variables
nInit = 10000;
leafStartPoints  = zeros(nInit,3);
incAngles        = zeros(nInit,1);
azAngles         = zeros(nInit,1);
leafNormal       = zeros(nInit,3);
leafDir          = zeros(nInit,3);
leafScaleFactors = zeros(nInit,3);

%% Sample leaf positions, orientations and sizes
disp('Sampling leaves from leaf distributions')
tic
iLeaf = 0;
leafArea = 0;
iVarExt = 0;
while leafArea < candidateArea
    % Increase leaf index
    iLeaf = iLeaf + 1;
    % Sampling leaf position
    accepted = 0;
    while accepted == 0
        if flagPCPositionSampling == true
            % Based on the point cloud sampling weight, draw wether the
            % position is sampled with point cloud density or LADD
            if rand(1) <= pcSamplingWeight
                % acceptedPCS = 0;
                % while acceptedPCS == 0
                %     % Proposal values for cylinder coordinates
                %     proposals = rand(1,3).*[maxHeight maxHorzDist 2*pi];
                %     % Bin indexes for the proposal coordinates
                %     ih = find((binEdgesH>=proposals(1)),1,'first') - 1;
                %     id = find((binEdgesD>=proposals(2)),1,'first') - 1;
                %     ic = find((binEdgesC>=proposals(3)),1,'first') - 1;
                %     % xyz-coordinates of the proposal point
                %     xyzProposal = rotation_matrix([0 0 1],proposals(3)) ...
                %                   *(proposals(2)*[0 1 0]') ... 
                %                   + [0 0 proposals(1)]';
                %     % Rejection sampling
                %     if inShape(shp,xyzProposal) && ...
                %        rand(1)*maxPCDensity < pcDensityInVolume(ih,id,ic)
                %         leafStartPoints(iLeaf,:) = proposals;
                %         acceptedPCS = 1;
                %     end
                % end
                % Sampling a volume block based on point cloud density
                u = rand(1);
                iPCD = find(cumulativeDistPCD>u,1,'first');
                inds = volumeInd(iPCD,:);
                hBin = binEdgesH(inds(1):(inds(1)+1));
                dBin = binEdgesD(inds(2):(inds(2)+1));
                cBin = binEdgesC(inds(3):(inds(3)+1));
                pointInShape = 0;
                while pointInShape == 0
                    proposals = rand(1,3) ...
                                .*[hBin(2)-hBin(1),dBin(2)-dBin(1), ...
                                   cBin(2)-cBin(1)] ...
                                + [hBin(1) dBin(1) cBin(1)];
                    xyzProposal = rotation_matrix([0 0 1],proposals(3)) ...
                                  *(proposals(2)*[0 1 0]') ... 
                                  + [0 0 proposals(1)]';
                    if inShape(shp,xyzProposal')
                        leafStartPoints(iLeaf,:) = xyzProposal';
                        pointInShape = 1;
                    end
                end
                % Calculate the structural variables of the sampled point
                hProposal = proposals(1)/maxHeight;
                cProposal = proposals(3);
                edgeDistance = stem_to_alphashape_edge(shp,hProposal, ...
                                                       cProposal, ...
                                                       maxHorzDist, ...
                                                       maxHeight, ...
                                                       stemCoordinates);
                dProposal = proposals(2)/edgeDistance;
                % Set the leaf acceptance flag to positive
                accepted = 1;
                % Skip the sampling of LADD
                continue
            end
        end
        % Proposal values for the variables
        hProposal = rand(1);
        dProposal = rand(1);
        cProposal = 2*pi*rand(1);
        % LADD value on propsal point
        funValues = [fDist_h(hProposal),fDist_d(dProposal), ...
                     fDist_c(cProposal)];
        vertValues = rand(1,3).*[maxfDist_h,maxfDist_d,maxfDist_c];
        if all(vertValues < funValues)
            % Probing the edge of alpha shape
            if flagStemCoordinates == true
                edgeValue = stem_to_alphashape_edge(shp,hProposal, ...
                                                    cProposal, ...
                                                    maxHorzDist, ...
                                                    maxHeight, ...
                                                    stemCoordinates);
            else
                edgeValue = stem_to_alphashape_edge(shp,hProposal, ...
                                                    cProposal, ...
                                                    maxHorzDist, ...
                                                    maxHeight);
            end
            % Check the if the distance to edge is nonzero
            if edgeValue == 0
                % The first two probe points are not inside the shape
                continue
            end
            edgeValue = yCoord(edgeIndex);
            % Position coordinates of the leaf start point
            hLeaf = maxHeight*hProposal;
            dLeaf = edgeValue*dProposal;
            cLeaf = cProposal;
            if flagStemCoordinates == true
                leafStartPoints(iLeaf,:) = (rotation_matrix([0 0 1],cLeaf) ...
                                           *[0 1 0]')'*dLeaf ...
                                           +[stemCen(1:2) hLeaf];
            else
                leafStartPoints(iLeaf,:) = (rotation_matrix([0 0 1],cLeaf) ...
                                           *[0 1 0]')'*dLeaf + [0 0 hLeaf];
            end
            accepted = 1;
        end
    end

    % Sampling leaf normal orientation with LOD
    % Leaf inclination angle
    switch dType_inc
        case 'dewit'
            % Sample inclination angle value with acceptance-rejection
            % sampling
            accepted = 0;
            while accepted == 0
                incProposal = rand(1)*pi/2;
                dParams = fun_inc_params(hProposal,dProposal,cProposal);
                funValue = f_inc(incProposal,dParams(1),dParams(2));
                vertValue = 1.5*rand(1); % all de Wit values are below 1.5
                if vertValue < funValue
                    incAngles(iLeaf) = incProposal;
                    accepted = 1;
                end
            end
        case 'beta'
            % Sample inclination angle by inverse transform sampling
            u = rand(1);
            dParams = fun_inc_params(hProposal,dProposal,cProposal);
            incAngles(iLeaf) = F_inc_inv(u,dParams(1),dParams(2));
    end
    % Leaf azimuth angle
    switch dType_az
        case 'vonmises'
            % Sample azimuth angle value with acceptance-rejection
            % sampling
            accepted = 0;
            while accepted == 0
                azProposal = rand(1)*2*pi;
                dParams = fun_az_params(hProposal,dProposal,cProposal);
                funValue = f_az(azProposal,dParams(1),dParams(2));
                vertValue = f_az(dParams(1),dParams(1),dParams(2))*rand(1);
                if vertValue < funValue
                    azAngles(iLeaf) = azProposal;
                    accepted = 1;
                end
            end
    end

    % Sampling leaf surface area with LSD
    switch dType_size
        case 'uniform'
            dParams = fun_size_params(hProposal,dProposal,cProposal);
            sampledArea = (dParams(2)-dParams(1))*rand(1) + dParams(1);
        case'normal'
            dParams = fun_size_params(hProposal,dProposal,cProposal);
            sampledArea = sqrt(dParams(2))*randn(1) + dParams(1);
    end
    leafArea = leafArea + sampledArea;
    % Store leaf scaling factor (same for all dimensions)
    leafScaleFactors(iLeaf,:) = (sampledArea/baseArea)*ones(1,3);

    % Unit vector of the leaf normal (y-axis assumed as north direction)
    leafNormal(iLeaf,:) = [sin(incAngles(iLeaf))*cos(azAngles(iLeaf)+pi/2),...
                           sin(incAngles(iLeaf))*sin(azAngles(iLeaf)+pi/2),...
                           cos(incAngles(iLeaf))];


    % Sampling leaf direction uniformly
    if abs(dot(leafNormal(iLeaf,:),[0 0 1])) > 0.999
        dirVec = [0 0 1];
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
disp('Adding leaves to the model without intersections')
tic
Leaves = add_leaves_to_alphashape(shp, ...
                                  Leaves, ...
                                  totalLeafArea, ....
                                  leafStartPoints, ...
                                  leafNormal, ...
                                  leafDir, ...
                                  leafScaleFactors);
toc
disp('Foliage generation finished')
