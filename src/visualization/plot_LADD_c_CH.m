function f = plot_LADD_c_CH(aShape,Leaves,TargetDistributions,varargin)

% Initialize values
nBins = 10;
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
dType = TargetDistributions.dTypeLADDc;
if ~any(strcmp(dType,{'','uniform','vonmises','vonmisesmixture'}))
    error("LADD compass direction distribution type not recognized.")
end

% Read target distribution parameters
p = TargetDistributions.pLADDc;

% Functions
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));

% Compass direction discretization
xx = 0:0.001:2*pi;

% Bins for the histogram of accepted leaves
binEdges = linspace(0,2*pi,nBins+1);

%% Plot the target distribution function
if dType ~= ""
    switch dType
        case 'uniform'
            yy = 1/(2*pi)*ones(size(xx));
        case 'vonmises'
            m = p(1); % mean
            k = p(2); % measure of concentration
            yy = fun_vonmises(xx,m,k);
        case 'vonmisesmixture'
            m1 = p(1); k1 = p(2); % parameters of the first distribution
            m2 = p(3); k2 = p(4); % parameters of the second distribution
            w = p(5); % mixture model weight
            yy = w*fun_vonmises(xx,m1,k1) + (1-w)*fun_vonmises(xx,m2,k2);
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
    plot(xx,yy,'r:','LineWidth',2,'DisplayName',"Target distribution")
end

%% Plot the histogram of accepted leaves

% Extracting leaf information
leafCount = Leaves.leaf_count;
leafScale = Leaves.leaf_scale;
leafBaseArea = Leaves.base_area;
leafStartPoints = Leaves.leaf_start_point;

% Area of each leaf (leaf scaling identical in every dimension)
leafAreas = (leafScale(:,1).^2)*leafBaseArea;

% initialize variable for compass direction wrt. the stem
compassDir = zeros(leafCount,1);

if flagStemCoordinates == true
    for iLeaf = 1:leafCount
        % Index of stem coordinate with z-value below the leaf height
        if stemCoordinates(end,3) < leafStartPoints(iLeaf,3)
            iSC = size(stemCoordinates,1);
        else
            iSC = find(stemCoordinates(:,3) > leafStartPoints(iLeaf,3),1);
        end
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
            compassDir(iLeaf) = acos(dot(unitStemToLeaf,[0 1 0]));
        else 
            compassDir(iLeaf) = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
        end
    end
else
    for iLeaf = 1:leafCount
        % Unit horizontal direction vector from the stem to leaf
        stemToLeaf = [leafStartPoints(iLeaf,1:2) 0];
        unitStemToLeaf = stemToLeaf/norm(stemToLeaf);
        % Counter-clockwise compass direction with respect to north 
        % direction
        if unitStemToLeaf(1) < 0
            compassDir(iLeaf) = acos(dot(unitStemToLeaf,[0 1 0]));
        else 
            compassDir(iLeaf) = 2*pi - acos(dot(unitStemToLeaf,[0 1 0]));
        end
    end
end

% Calculate weighted histogram for the leaf area wrt. height
leafHist = zeros(nBins,1);
for iLeaf = 1:leafCount
    for iBin = 1:nBins
        if compassDir(iLeaf) > binEdges(iBin) && ...
            compassDir(iLeaf) <= binEdges(iBin+1)
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
custom_bar_plot(binEdges,leafHistFD,'FaceColor','b','FaceAlpha',0.3,...
                'DisplayName','Accepted leaf area')
xlabel("compass direction")
ylabel("leaf area frequency density [m^2]")
axis tight
legend('Location','southeast')


end