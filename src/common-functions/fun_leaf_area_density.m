function relAreas = fun_leaf_area_density(CylinderParameters, ...
                                          TargetDist)

%% Function definitions

fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_weibull = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));

%% Distribution functions and parameters 

% Distribution function and parameters for relative height
dTypeLADDh = TargetDist.dTypeLADDh;
hParams    = TargetDist.hParams;
switch dTypeLADDh
    case 'uniform'
        fDist_h = @(h) 1*ones(size(h));
    case 'polynomial'
        fDist_h = @(h) polyval(hParams,h);
    case 'polynomialmixture'
        nP = (length(hParams)-1)/2; % number of polynomial coefficients
        p1 = hParams(1:nP); % coefficients of the first polynomial
        p2 = hParams((nP+1):(2*nP)); % coefficients of the second polynom.
        w = hParams(end); % mixture model weight
        fDist_h = @(h) w*polyval(p1,h) + (1-w)*polyval(p2,h);
    case 'weibull'
        l = hParams(1); % scale parameter
        k = hParams(2); % shape parameter
        fDist_h = @(h) fun_weibull(h,l,k);
    case 'weibullmixture'
        l1 = hParams(1); k1 = hParams(2); % parameters of the first dist.
        l2 = hParams(3); k2 = hParams(4); % parameters of the second dist.
        w = hParams(5); % mixture model weight
        fDist_h = @(h) w*fun_weibull(h,l1,k1) + (1-w)*fun_weibull(h,l2,k2);
    case 'beta'
        a = hParams(1);
        b = hParams(2);
        fDist_h = @(h) fun_beta(h,a,b);
    case 'betamixture'
        a1 = hParams(1); b1 = hParams(2); % parameters of the first dist.
        a2 = hParams(3); b2 = hParams(4); % parameters of the second dist.
        w = hParams(5); % mixture model weight
        fDist_h = @(h) w*fun_beta(h,a1,b1) + (1-w)*fun_beta(h,a2,b2);
end

% Distribution function and parameters for relative distance along
% sub-branch
dTypeLADDd = TargetDist.dTypeLADDd;
dParams    = TargetDist.dParams;
switch dTypeLADDd
    case 'uniform'
        fDist_d = @(d) 1*ones(size(d));
    case 'polynomial'
        fDist_d = @(d) polyval(dParams,d);
    case 'polynomialmixture'
        nP = (length(dParams)-1)/2; % order of polynomial
        p1 = dParams(1:nP); % coefficients of the first polynomial
        p2 = dParams((nP+1):(2*nP)); % coefficients of the second polynom.
        w = dParams(end); % mixture model weight
        fDist_d = @(d) w*polyval(p1,d) + (1-w)*polyval(p2,d);
    case 'weibull'
        l = dParams(1); % scale parameter
        k = dParams(2); % shape parameter
        fDist_d = @(d) fun_weibull(d,l,k);
    case 'weibullmixture'
        l1 = dParams(1); k1 = dParams(2); % parameters of the first dist.
        l2 = dParams(3); k2 = dParams(4); % parameters of the second dist.
        w = dParams(5); % mixture model weight
        fDist_d = @(d) w*fun_weibull(d,l1,k1) + (1-w)*fun_weibull(d,l2,k2);
    case 'beta'
        a = dParams(1);
        b = dParams(2);
        fDist_d = @(d) fun_beta(d,a,b);
    case 'betamixture'
        a1 = dParams(1); b1 = dParams(2); % parameters of the first dist.
        a2 = dParams(3); b2 = dParams(4); % parameters of the second dist.
        w = dParams(5); % mixture model weight
        fDist_d = @(d) w*fun_beta(d,a1,b1) + (1-w)*fun_beta(d,a2,b2);
end

% Distribution function and parameters for compass direction
dTypeLADDc = TargetDist.dTypeLADDc;
cParams    = TargetDist.cParams;
switch dTypeLADDc
    case 'uniform'
        fDist_c = @(c) 1/(2*pi)*ones(size(c));
    case 'vonmises'
        m = cParams(1); % mean
        k = cParams(2); % measure of concentration
        fDist_c = @(c) fun_vonmises(c,m,k);
    case 'vonmisesmixture'
        m1 = cParams(1); k1 = cParams(2); % parameters of the first dist.
        m2 = cParams(3); k2 = cParams(4); % parameters of the second dist.
        w = cParams(5); % mixture model weight
        fDist_c = @(c) w*fun_vonmises(c,m1,k1) ...
                       + (1-w)*fun_vonmises(c,m2,k2);
end

%% QSM cylinder properties

relativeHeight = CylinderParameters.relative_height;
startPoint     = CylinderParameters.start_point;
branchIndex    = CylinderParameters.branch_index;
cylinderLength = CylinderParameters.length;
cylinderAxis   = CylinderParameters.axis;

nCylinders = length(relativeHeight);

%% Bins for distribution of leaf area

% Heightwise bins for distributing leaf area
if isfield(TargetDist,'nBins_h')
    nBinsH = TargetDist.nBinsLADDh;
else
    nBinsH = 10;
end
% Distancewise bins for distributing leaf area
if isfield(TargetDist,'nBins_d')
    nBinsD = TargetDist.nBinsLADDd;
else
    nBinsD = 10;
end
% Directionwise bins for distributing leaf area
if isfield(TargetDist,'nBins_c')
    nBinsC = TargetDist.nBinsLADDc;
else
    nBinsC = 10;
end


%% Partition of leaf surface area with respect to height

% Heightwise bins
binEdges_h = generate_bin_edges([0 1],nBinsH,fDist_h);

% Calculate normalized distribution function values for each height bin
binAreas_h = zeros(nBinsH,1);
for iBin = 1:nBinsH
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

% Set cylinder endpoint as the reference point for position in branch
referencePoint = "end";

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
    if referencePoint == "end"
        % Relative position in subbranch
        relDistInSubBranch = bcLengthsCumulative/bcLengthsCumulative(end);
    elseif referencePoint == "middle"
        % Cumulative distance of cylinder midpoints in sub-branch
        midPointCumulative = bcLengthsCumulative - 0.5*bcLengths;
        % Relative position in subbranch
        relDistInSubBranch = midPointCumulative/bcLengthsCumulative(end);
    end
    k = 1;
    for j = 1:length(bcIndexes)
        bcInd = bcIndexes(j);
        relativeDistanceFromBase(bcInd) = relDistInSubBranch(k);
        k = k + 1;
    end
end

% Distancewise bins
binEdges_d = generate_bin_edges([0 1],nBinsD,fDist_d);

% Calculate normalized distribution function values for each distance bin
binAreas_d = zeros(nBinsD,1);
for iBin = 1:nBinsD
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
binAreas_d = binAreas_d/sum(binAreas_d);

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

% Cylinder middle points
midPoint = startPoint + 0.5*cylinderLength.*cylinderAxis;
% % Coordinate of stem base as origin on the xy-plane
% xyOrigin = startPoint(1,1:2);
% Set the mean value of cylinder locations on xy-plane as the origin
xyOrigin = mean(midPoint(:,1:2));
% Coordinates of the cylinder middle points on xy-plane
cylCoord = midPoint(:,1:2);
% Unit vectors pointing the direction of cylinder middle point coordinates
cylDir = zeros(nCylinders,2);
for j = 1:nCylinders
    if norm(cylCoord(j,:)-xyOrigin) < 1e-6 % midpoint right above set origin
        % Allocate cylinder direction randomly
        temp = (rotation_matrix([0 0 1],rand(1)*2*pi)*[0 1 0]')';
        cylDir(j,:) = temp(1:2)/norm(temp(1:2));
    else 
        % Find cylinder midpoint direction in xy plane
        temp = cylCoord(j,:) - xyOrigin;
        cylDir(j,:) = temp/norm(temp);
    end
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
binEdges_a = generate_bin_edges([0 2*pi],nBinsC,fDist_c);

% Calculate normalized distribution function values for each direction bin
binAreas_a = zeros(nBinsC,1);
for iBin = 1:nBinsC
    xTemp_a = linspace(binEdges_a(iBin),binEdges_a(iBin+1),1000);
    yTemp_a = fDist_c(xTemp_a);
    % Integrate function values over bin width
    binAreas_a(iBin) = trapz(xTemp_a,yTemp_a);
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
for iBin_h = 1:nBinsH
    for iBin_d = 1:nBinsD
        for iBin_a = 1:nBinsC
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