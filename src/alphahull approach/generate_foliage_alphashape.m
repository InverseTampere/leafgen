function [Leaves,shp] = generate_foliage_alphashape(treePointCloud, ...
                                                    TargetDistributions, ...
                                                    totalLeafArea, ...
                                                    vertices, ...
                                                    tris, ...
                                                    varargin ...
                                                    )
%% Initialize values
alpha = 1;

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
%% Leaf area density distribution functions
% Distribution function and parameters for relative height
dType_h = TargetDistributions.dType_h;
p_h     = TargetDistributions.p_h;
% nBins_h = TargetDistributions.nBins_h;
switch dType_h
    case 'poly3'
        fDist_h = @(h) p_h(1)*h.^3 + p_h(2)*h.^2 + p_h(3)*h + p_h(4);
        [~,maxfDist_h] = fminbnd(-fDist_h,0,1);
    case 'weibull'
        % p = [scale, shape]
        sc = p_h(1);
        sh = p_h(2);
        fDist_h = @(h) (sh/sc)*(h/sc).^(sh-1).*exp(-(h/sc).^sh);
        if sh <= 1
            maxfDist_h = 5; % infinite mode at 0
        elseif sc*((sh-1)/sh).^(1/sh) < 1
            maxfDist_h = fDist_h(sc*((sh-1)/sh).^(1/sh)); % mode between 0 and 1
        else
            maxfDist_h = fDist_h(1); % mode at 1
        end
    case 'beta'
        a = p_h(1);
        b = p_h(2);
        fDist_h = @(h) (1/beta(a,b))*h.^(a-1).*(1-h).^(b-1);
        if a == b
            maxfDist_h = fDist_h(0.5);
        elseif a > 1 && b > 1
            maxfDist_h = fDist_h((a-1)/(a+b-2));
        elseif a < 1 && b < 1
            maxfDist_h = 10; % bimodal with infinite modes
        elseif a == 1 && b > 1
            maxfDist_h = fDist_h(0);
        elseif a > 1 && b == 1
            maxfDist_h = fDist_h(1);
        else
            maxfDist_h = 10; % infinite mode at 0 or 1
        end
        maxfDist_h = min(maxfDist_h,10);
end
% Distribution function and parameters for relative distance along
% sub-branch
dType_d = TargetDistributions.dType_d;
p_d     = TargetDistributions.p_d;
% nBins_d = TargetDistributions.nBins_d;
switch dType_d
    case 'poly4'
        fDist_d = @(d) p_d(1)*d.^4 + p_d(2)*d.^3 + p_d(3)*d.^2 ...
                       + p_d(4)*d + p_d(5);
        [~,maxfDist_d] = fminbnd(-fDist_d,0,1);
    case 'beta'
        a = p_d(1);
        b = p_d(2);
        fDist_d = @(d) (1/beta(a,b))*d.^(a-1).*(1-d).^(b-1);
        if a == b
            maxfDist_d = fDist_d(0.5);
        elseif a > 1 && b > 1
            maxfDist_d = fDist_d((a-1)/(a+b-2));
        elseif a < 1 && b < 1
            maxfDist_d = 10; % bimodal with infinite modes
        elseif a == 1 && b > 1
            maxfDist_d = fDist_d(0);
        elseif a > 1 && b == 1
            maxfDist_d = fDist_d(1);
        else
            maxfDist_d = 10; % infinite mode at 0 or 1
        end
        maxfDist_d = min(maxfDist_d,10);
end
% Distribution function and parameters for compass direction
dType_c = TargetDistributions.dType_c;
p_c     = TargetDistributions.p_c;
% nBins_c = TargetDistributions.nBins_c;
switch dType_c
    case 'vonmises'
        m = p_c(1); % mean
        k = p_c(2); % measure of concentration
        fDist_c = @(c) exp(k*cos(c-m))./(2*pi*besseli(0,k));
        maxfDist_c = fDist_c(m);
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
        vertValues = rand(1)*[maxfDist_h,maxfDist_d,maxfDist_c];
        if all(vertValues < funValues)
            % Probing the edge of alpha shape
            nPP = 100;
            yCoord = maxHorzDist*linspace(0,1,nPP)';
            initPP = [zeros(nPP,1) yCoord zeros(nPP,1)];
            probePoints = (rotation_matrix([0 0 1],cProposal)*initPP')' ...
                + [zeros(nPP,2) maxHeight*hProposal*ones(nPP,1)];
            if flagStemCoordinates
                iSC = find(stemCoordinates(:,3) > maxHeight*hProposal,1);
                relPos = (maxHeight*hProposal-stemCoordinates(iSC-1,3)) ...
                        /(stemCoordinates(iSC,3)-stemCoordinates(iSC-1,3));
                stemCen = relPos*(stemCoordinates(iSC,:) ...
                                  -stemCoordinates(iSC-1,:)) ...
                          + stemCoordinates(iSC-1,:);
                probePoints = probePoints + [stemCen(1:2) 0];
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
            leafStartPoints(iLeaf,:) = (rotation_matrix([0 0 1],cLeaf) ...
                                        *[0 1 0]')'*dLeaf + [0 0 hLeaf];
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
    if iLeaf > length(leafStartPoints(:,1))
        leafStartPoints  = [leafStartPoints;  zeros(nInit,3)]; %#ok<AGROW>
        incAngles        = [incAngles;        zeros(nInit,1)]; %#ok<AGROW>
        azAngles         = [azAngles;         zeros(nInit,1)]; %#ok<AGROW>
        leafNormal       = [leafNormal;       zeros(nInit,3)]; %#ok<AGROW>
        leafDir          = [leafDir;          zeros(nInit,3)]; %#ok<AGROW>
        leafScaleFactors = [leafScaleFactors; zeros(nInit,3)]; %#ok<AGROW>
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
