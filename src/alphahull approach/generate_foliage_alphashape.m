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
        maxfDist_h = polynomial_upper_limit(p_h);
    case 'polynomialmixturemodel'
        nP = (length(p_h)-1)/2; % order of polynomial
        p1 = p_h(1:nP); % coefficients of the first polynomial
        p2 = p_h((nP+1):(2*nP)); % coefficients of the second polynomial
        w = p_h(end); % mixture model weight
        fDist_h = @(h) w*polyval(p1,h) + (1-w)*polyval(p2,h);
        maxfDist_h = polynomial_upper_limit(w*p1+(1-w)*p2);
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
        maxfDist_d = polynomial_upper_limit(p_d);
    case 'polynomialmixturemodel'
        nP = (length(p_d)-1)/2; % order of polynomial
        p1 = p_d(1:nP); % coefficients of the first polynomial
        p2 = p_d((nP+1):(2*nP)); % coefficients of the second polynomial
        w = p_d(end); % mixture model weight
        fDist_d = @(d) w*polyval(p1,d) + (1-w)*polyval(p2,d);
        maxfDist_d = polynomial_upper_limit(w*p1+(1-w)*p2);
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
tic
shp = alphaShape(treePointCloud,alpha);
toc
%% Extreme points of point cloud
% Maximum horizontal distance from origin to point cloud member
horzDistances = sqrt(treePointCloud(:,1).^2 + treePointCloud(:,2).^2);
maxHorzDist = max(horzDistances);

% Maximum height of the point cloud members
maxHeight = max(treePointCloud(:,3));

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

%% Sample leaf positions and orientations
iLeaf = 0;
leafArea = 0;
iVarExt = 0;
while leafArea < candidateArea
    % Increase leaf index
    iLeaf = iLeaf + 1;
    % Sampling leaf position with LADD
    accepted = 0;
    while accepted == 0
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
            if flagStemCoordinates
                nPP = 2*100;
                yCoord = 2*maxHorzDist*linspace(0,1,nPP)';
                initPP = [zeros(nPP,1) yCoord zeros(nPP,1)];
                probePoints = (rotation_matrix([0 0 1],cProposal)*initPP')' ...
                    + [zeros(nPP,2) maxHeight*hProposal*ones(nPP,1)];
                iSC = find(stemCoordinates(:,3) > maxHeight*hProposal,1);
                relPos = (maxHeight*hProposal-stemCoordinates(iSC-1,3)) ...
                        /(stemCoordinates(iSC,3)-stemCoordinates(iSC-1,3));
                stemCen = relPos*(stemCoordinates(iSC,:) ...
                                  -stemCoordinates(iSC-1,:)) ...
                          + stemCoordinates(iSC-1,:);
                probePoints = probePoints + [stemCen(1:2) 0];
            else
                nPP = 100;
                yCoord = maxHorzDist*linspace(0,1,nPP)';
                initPP = [zeros(nPP,1) yCoord zeros(nPP,1)];
                probePoints = (rotation_matrix([0 0 1],cProposal)*initPP')' ...
                    + [zeros(nPP,2) maxHeight*hProposal*ones(nPP,1)];
            end
            tf = inShape(shp,probePoints);
            edgeIndex = find(~tf,1,'first') - 1;
            if edgeIndex < 2
                % The first two probe points are not inside the shape
                continue
            end
            edgeValue = yCoord(edgeIndex);
            % Position coordinates of the leaf start point
            hLeaf = maxHeight*hProposal;
            dLeaf = edgeValue*dProposal;
            cLeaf = cProposal;
            if flagStemCoordinates
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

%% Add leaves to the shape without intersections
Leaves = add_leaves_to_alphashape(shp, ...
                                  Leaves, ...
                                  totalLeafArea, ....
                                  leafStartPoints, ...
                                  leafNormal, ...
                                  leafDir, ...
                                  leafScaleFactors);
