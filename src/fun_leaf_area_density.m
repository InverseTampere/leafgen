function relAreas = fun_leaf_area_density(CylinderParameters, ...
                                          TargetDistributions)

% Distribution function and parameters for relative height
dType_h = TargetDistributions.dType_h;
p_h     = TargetDistributions.p_h;
nBins_h = TargetDistributions.nBins_h;
switch dType_h
    case 'poly3'
        fDist_h = @(h) p_h(1)*h.^3 + p_h(2)*h.^2 + p_h(3)*h + p_h(4);
    case 'weibull'
        % p = [scale, shape]
        sc = p_h(1);
        sh = p_h(2);
        fDist_h = @(h) (sh/sc)*(h/sc).^(sh-1).*exp(-(h/sc).^sh);
    case 'beta'
        fDist_h = @(h) (1/beta(p_h(1),p_h(2)))*h.^(p_h(1)-1) ...
                       .*(1-h).^(p_h(2)-1);
end

% Distribution function and parameters for relative distance along
% sub-branch
dType_d = TargetDistributions.dType_d;
p_d     = TargetDistributions.p_d;
nBins_d = TargetDistributions.nBins_d;
switch dType_d
    case 'poly4'
        fDist_d = @(d) p_d(1)*d.^4 + p_d(2)*d.^3 + p_d(3)*d.^2 ...
                       + p_d(4)*d + p_d(5);
    case 'beta'
        fDist_d = @(d) (1/beta(p_d(1),p_d(2)))*d.^(p_d(1)-1) ...
                       .*(1-d).^(p_d(2)-1);
end

% Distribution function and parameters for compass direction
dType_c = TargetDistributions.dType_c;
p_c     = TargetDistributions.p_c;
nBins_c = TargetDistributions.nBins_c;
switch dType_c
    case 'vonmises'
        fDist_c = @(c) exp(p_c(2)*cos(c-p_c(1)))./(2*pi*besseli(0,p_c(2)));
end

%% QSM cylinder properties

relativeHeight = CylinderParameters.relative_height;
startPoint     = CylinderParameters.start_point;
branchIndex    = CylinderParameters.branch_index;
cylinderLength = CylinderParameters.length;

nCylinders = length(relativeHeight);

%% Partition of leaf surface area with respect to height

% Heightwise bins
binEdges_h = generate_bin_edges([0 1],nBins_h,fDist_h);

% Calculate normalized distribution function values for each height bin
binAreas_h = zeros(nBins_h,1);
for iBin = 1:nBins_h
    xTemp_h = linspace(binEdges_h(iBin),binEdges_h(iBin+1),1000);
    yTemp_h = fDist_h(xTemp_h);
    % Prevent infinite values at the edges of the interval
    if yTemp_h(1) == Inf
        yTemp_h(1) = yTemp_h(2) + (yTemp_h(2)-yTemp_h(3));
    end
    if yTemp_h(end) == Inf
        yTemp_h(end) = yTemp_h(end-1) + (yTemp_h(end-1)-yTemp_h(end-2));
    end
    % Integrate function values over bin width
    binAreas_h(iBin) = trapz(xTemp_h,yTemp_h);
end
% Normalize bin areas
binAreas_h = binAreas_h/sum(binAreas_h);

% Sort tree cylinders to height bins
binIndex_h = zeros(nCylinders,1);
for j = 1:nCylinders
    k = 1;
    while relativeHeight(j) > binEdges_h(k+1)
        k = k + 1;
    end
    binIndex_h(j) = k;
end

%% Partition of leaf surface area with respect to along-branch distance 
%  from the base of subbranch

% Relative along-branch distances from base
relativeDistanceFromBase = zeros(nCylinders,1);
indexVector = (1:1:nCylinders)';
for iBranch = 1:max(branchIndex) % stem index 0 is skipped automatically
    % Indexes of QSM cylinders belonging in branch
    bcIndexes = indexVector(branchIndex == iBranch);
    if sum(bcIndexes) == 0
        % the branch has no cylinders
        continue
    end
    bcIndexes = nonzeros(bcIndexes);
    % Lengths of corresponding cylinders
    bcLengths = zeros(length(bcIndexes),1);
    for j = 1:length(bcLengths)
        bcLengths(j) = cylinderLength(bcIndexes(j));
    end
    % Cumulative sum of the lengths
    bcLengthsCumulative = cumsum(bcLengths);
    % Cumulative distance of cylinder midpoints in sub-branch
    midPointCumulative = bcLengthsCumulative - 0.5*bcLengths;
    % Relative position in subbranch
    relDistInSubBranch = midPointCumulative/bcLengthsCumulative(end);
    k = 1;
    for j = 1:length(bcIndexes)
        bcInd = bcIndexes(j);
        relativeDistanceFromBase(bcInd) = relDistInSubBranch(k);
        k = k + 1;
    end
end

% Distancewise bins
binEdges_d = generate_bin_edges([0 1],nBins_d,fDist_d);

% Calculate normalized distribution function values for each distance bin
binAreas_d = zeros(nBins_d,1);
for iBin = 1:nBins_d
    xTemp_d = linspace(binEdges_d(iBin),binEdges_d(iBin+1),1000);
    yTemp_d = fDist_d(xTemp_d);
    % Prevent infinite values at the edges of the interval
    if yTemp_d(1) == Inf
        yTemp_d(1) = yTemp_d(2) + (yTemp_d(2)-yTemp_d(3));
    end
    if yTemp_d(end) == Inf
        yTemp_d(end) = yTemp_d(end-1) + (yTemp_d(end-1)-yTemp_d(end-2));
    end
    % Integrate function values over bin width
    binAreas_d(iBin) = trapz(xTemp_d,yTemp_d);
end
% Normalize bin areas
binAreas_h = binAreas_h/sum(binAreas_h);

% Sort tree cylinders to distance bins
binIndex_d = zeros(nCylinders,1);
for j = 1:nCylinders
    k = 1;
    while relativeDistanceFromBase(j) > binEdges_d(k+1)
        k = k + 1;
    end
    binIndex_d(j) = k;
end

%% Partition of leaf surface area with respect to angular direction

% % Coordinate of stem base as origin on the xy-plane
% xyOrigin = startPoint(1,1:2);
% Set the mean value of cylinder locations on xy-plane as the origin
xyOrigin = mean(startPoint(:,1:2));
% Coordinates of the cylinder start points on xy-plane
cylCoord = startPoint(:,1:2);
% Unit vectors pointing the direction of cylinder start point coordinates
cylDir = zeros(nCylinders,2);
for j = 1:nCylinders
    temp = cylCoord(j,:) - xyOrigin;
    cylDir(j,:) = temp/norm(temp);
end
% North set as the zero angle
northDir = [0 1];
% Angles of the cylinder direction unit vectors w.r.t north
phiCyl = zeros(nCylinders,1);
for j = 1:nCylinders
    if cylDir(j,1) <= 0
        phiCyl(j) = acos(dot(northDir,cylDir(j,:)));
    else
        phiCyl(j) = 2*pi - acos(dot(northDir,cylDir(j,:)));
    end
end

% Directionwise bins
binEdges_a = generate_bin_edges([0 2*pi],nBins_c,fDist_c);

% Calculate normalized distribution function values for each direction bin
binAreas_a = zeros(nBins_c,1);
for iBin = 1:nBins_c
    xTemp_a = linspace(binEdges_a(iBin),binEdges_a(iBin+1),1000);
    yTemp_a = fDist_c(xTemp_a);
    % Integrate function values over bin width (polar integration)
    binAreas_a(iBin) = 0.5*trapz(xTemp_a,yTemp_a.^2);
end
% Normalize bin areas
binAreas_a = binAreas_a/sum(binAreas_a);

% Sort tree cylinders to directional bins
binIndex_a = zeros(nCylinders,1);
for j = 1:nCylinders
    k = 1;
    while phiCyl(j) > binEdges_a(k+1)
        k = k + 1;
    end
    binIndex_a(j) = k;
end

%% Assigned relative leaf areas for the cylinders

% Block area factors for each cylinder
blockArea = zeros(nCylinders,1);
for iBin_h = 1:nBins_h
    for iBin_d = 1:nBins_d
        for iBin_a = 1:nBins_c
            % Indexes of cylinders with desired bin indexes
            bool_h = (binIndex_h == iBin_h);
            bool_d = (binIndex_d == iBin_d);
            bool_a = (binIndex_a == iBin_a);
            blockInds = logical(bool_h.*bool_d.*bool_a);
            % Accounting for cylinder length
            binCylLengths = cylinderLength(blockInds);
            binCylLengthsNorm = binCylLengths/sum(binCylLengths);
            % Distribute the leaf area of the bin
            blockArea(blockInds) = (binAreas_h(iBin_h) ...
                                    *binAreas_d(iBin_d) ...
                                    *binAreas_a(iBin_a) ...
                                    *binCylLengthsNorm);
        end
    end
end

% Normalized leaf area factors
relAreas = blockArea/sum(blockArea);


end