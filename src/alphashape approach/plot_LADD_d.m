function f = plot_LADD_d(aShape,Leaves,TargetDistributions,varargin)

% Initialize values
nBins = 10;
pCloud = aShape.Points;
flagStemCoordinates = false;

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'nbins'
                nBins = varargin{i+1};
            case 'stemcoordinates'
                stemCoordinates = varargin{i+1};
                flagStemCoordinates = true;
        end
    end
    i = i + 1;
end

% Issue warning if there is no stem coordinate information
if flagStemCoordinates == false
    warning("No stem coordinate information supplied, assuming z-axis as the stem.")
end

% Initialize figure object
f = figure; clf, hold on

% Check validity of distribution function type
dType = TargetDistributions.dType_d;
if ~any(strcmp(dType,{'','uniform','polynomial', ...
        'polynomialmixturemodel','weibull','weibullmixturemodel', ...
        'beta','betamixturemodel'}))
    error("LADD distance from stem distribution type not recognized.")
end

% Read target distribution parameters
p = TargetDistributions.p_d;

% Functions
fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_weibull = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);

% Relative distance discretization
xx = 0:0.001:1;

% Bins for the histogram of accepted leaves
binEdges = linspace(0,1,nBins+1);

%% Plot the target distribution function

if dType ~= ""
    switch dType
        case 'uniform'
            yy = ones(size(xx));
        case  'polynomial'
            yy = polyval(p,xx);
        case 'polynomialmixturemodel'
            nP = (length(p)-1)/2; % number of polynomial coefficients
            p1 = p(1:nP); % coefficients of the first polynomial
            p2 = p((nP+1):(2*nP)); % coefficients of the second polynomial
            w = p(end); % mixture model weight
            yy = w*polyval(p1,xx) + (1-w)*polyval(p2,xx);
        case 'weibull'
            l = p(1); % scale parameter
            k = p(2); % shape parameter
            yy = fun_weibull(xx,l,k);
        case 'weibullmixturemodel'
            l1 = p(1); k1 = p(2); % parameters of the first distribution
            l2 = p(3); k2 = p(4); % parameters of the second distribution
            w = p(5); % mixture model weight
            yy = w*fun_weibull(xx,l1,k1) + (1-w)*fun_weibull(xx,l2,k2);
        case 'beta'
            a = p(1);
            b = p(2);
            yy = fun_beta(xx,a,b);
        case 'betamixturemodel'
            a1 = p(1); b1 = p(2); % parameters of the first distribution
            a2 = p(3); b2 = p(4); % parameters of the second distribution
            w = p(5); % mixture model weight
            yy = w*fun_beta(xx,a1,b1) + (1-w)*fun_beta(xx,a2,b2);
    end
    % Prevent infinite values at the edges of the interval
    if yy(1) == Inf
        yy(1) = yy(2) + (yy(2)-yy(3));
    end
    if yy(end) == Inf
        yy(end) = yy(end-1) + (yy(end-1)-yy(end-2));
    end
    % Normalization
    yy = yy/trapz(xx,yy);
    plot(xx,yy,'r-','LineWidth',2,'DisplayName',"Target distribution")
end

%% Plot the histogram of accepted leaves

% Extracting leaf information
leafCount = Leaves.leaf_count;
leafScale = Leaves.leaf_scale;
leafBaseArea = Leaves.base_area;
leafStartPoints = Leaves.leaf_start_point;

% Area of each leaf (leaf scaling identical in every dimension)
leafAreas = (leafScale(:,1).^2)*leafBaseArea;

% Maximum horizontal distance of point cloud from origin
maxHorzDist = max(sqrt(sum(pCloud(:,1:2).^2,2)));

% Initialize variable for relative distance from stem
relDistFromStem = zeros(leafCount,1);

% Variables for missed leaves due to problems in defining distance to the
% alphashape edge
nMissed = 0;
indMissed = [];
areaMissed = 0;

if flagStemCoordinates == true % Stem coordinates supplied as input
    for iLeaf = 1:leafCount
        % Index of stem coordinate with z-value below the leaf height
        iSC = find(stemCoordinates(:,3) > leafStartPoints(iLeaf,3),1);
        % Relative z-position of leaf between the stem coordinate nodes
        relPos = (leafStartPoints(iLeaf,3)- stemCoordinates(iSC-1,3)) ...
                 /(stemCoordinates(iSC,3) - stemCoordinates(iSC-1,3));
        % Location of the stem center on the height of the leaf
        stemCen = relPos*(stemCoordinates(iSC,:) ...
                          -stemCoordinates(iSC-1,:)) ...
                  + stemCoordinates(iSC-1,:);
        % Horizontal distance from stem to leaf
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0] - [stemCen(1:2) 0];
        distStemToLeaf = norm(stemToLeaf); 
        % Counter-clockwise angle between north and leaf direction wrt. the
        % stem
        unitStemToLeaf = stemToLeaf/distStemToLeaf;
        if unitStemToLeaf(1) < 0
            compassDir = acos(dot(unitStemToLeaf,[0 1 0]));
        else 
            compassDir = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
        end
        % Probing the edge of alpha shape
        nPP = 2*100;
        yCoord = 2*maxHorzDist*linspace(0,1,nPP)';
        initPP = [zeros(nPP,1) yCoord zeros(nPP,1)];
        probePoints = (rotation_matrix([0 0 1],compassDir)*initPP')' ...
                      + stemCen;
        tf = inShape(aShape,probePoints);
        edgeIndex = find(~tf,1,'first') - 1;
        if edgeIndex < 1
            % The stem is outside alphashape
            nMissed = nMissed + 1;
            indMissed = [indMissed, iLeaf];
            areaMissed = areaMissed ...
                         + Leaves.base_area*Leaves.leaf_scale(iLeaf,2).^2;
            continue
        end
        edgeValue = yCoord(edgeIndex);
        % The relative horizontal distance from the stem
        relDistFromStem(iLeaf) = distStemToLeaf/edgeValue;
    end
else % Assuming stem to be the z-axis
    for iLeaf = 1:leafCount
        % Unit direction vector on xy-plane
        xyDir = leafStartPoints(iLeaf,1:2) ...
                /sqrt(sum(leafStartPoints(iLeaf,1:2).^2));
        % Probing the edge of alpha shpe
        nPP = 100;
        probeDist = maxHorzDist*linspace(0,1,nPP)';
        probePoints = [probeDist*xyDir, ...
                       leafStartPoints(iLeaf,3)*ones(nPP,1)];
        tf = inShape(aShape,probePoints);
        edgeIndex = find(~tf,1,'first') - 1;
        if edgeIndex < 1
            % The stem is outside alphashape
            nMissed = nMissed + 1;
            indMissed = [indMissed, iLeaf];
            areaMissed = areaMissed ...
                         + Leaves.base_area*Leaves.leaf_scale(iLeaf,2).^2;
            continue
        end
        edgeValue = norm(probePoints(edgeIndex,:));
        % The relative horizontal distance from the stem
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0];
        relDistFromStem(iLeaf) = norm(stemToLeaf)/edgeValue;
    end
end

if nMissed > 0
    warning('%d leaves not included in the \"relative distance from stem\" histogram due to the defined stem being outside alphashape. (Total missed area %.2f m^2)',nMissed,areaMissed)
    leafCount = leafCount - nMissed;
    relDistFromStem(indMissed) = [];
end

% Calculate weighted histogram for the leaf area wrt. height
leafHist = zeros(nBins,1);
for iLeaf = 1:leafCount
    for iBin = 1:nBins
        if relDistFromStem(iLeaf) > binEdges(iBin) && ...
            relDistFromStem(iLeaf) <= binEdges(iBin+1)
               leafHist(iBin) = leafHist(iBin) + leafAreas(iLeaf);
        end
    end
end

% Calculate accepted leaf area frequency density in bins
leafHistFD = zeros(nBins,1);
for iBin = 1:nBins
    leafHistFD(iBin) = leafHist(iBin)/(binEdges(iBin+1)-binEdges(iBin));
end
% Normalization
leafHistFD = leafHistFD/sum(leafHistFD);

% Divide with bin widths so that total bar area equals to 1
leafHistFD = leafHistFD./diff(binEdges);

% Plotting the histogram 
custom_bar_plot(binEdges,leafHistFD,'FaceColor','b','FaceAlpha',0.5,...
                'DisplayName','Accepted leaf area')
xlabel("relative distance from stem")
ylabel("leaf area frequency density [m^2]")
axis tight
legend('Location','southeast')

end