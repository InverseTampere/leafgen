function f = plot_LADD_h(aShape,Leaves,TargetDistributions,varargin)

% Initialize values
nBins = 10;
pCloud = aShape.Points;

% Check additional parameters
i = 1;
NArg = numel(varargin);
while i <= NArg
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'nbins'
                nBins = varargin{i+1};
        end
    end
    i = i + 1;
end

% Initialize figure object
f = figure; clf, hold on

% Check validity of distribution function type
dType = TargetDistributions.dType_h;
if dType ~= "poly3" && dType ~= "weibull" && dType ~= "beta" && ...
   dType ~= "none"
    error("LADD height distribution type not recognized.")
end

% Read target distribution parameters
p = TargetDistributions.p_h;

% Functions
fun_weibull = @(p,x) (p(2)/p(1))*(x/p(1)).^(p(2)-1).*exp(-(x/p(1)).^p(2));
fun_beta = @(p,x) (1/beta(p(1),p(2)))*x.^(p(1)-1).*(1-x).^(p(2)-1);

% Relative height discretization
xx = 0:0.001:1;

% Bins for the histogram of accepted leaves
binEdges = linspace(0,1,nBins+1);

%% Plot the target distribution function
if dType ~= "none"
    if dType == "poly3"
        % Normalization
        q = polyint(p);
        yy = polyval(p,xx)/diff(polyval(q,[0,1]));
        % Plotting the curve
        plot(xx,yy,'r-','LineWidth',2,'DisplayName',"Target distribution")
    elseif dType == "weibull"
        yy = fun_weibull(p,xx);
        % Normalization
        yy = yy/trapz(xx,yy);
        % Plotting the curve
        plot(xx,yy,'r-','LineWidth',2,'DisplayName',"Target distribution")
    elseif dType == "beta"
        yy = fun_beta(p,xx);
        % Prevent infinite values at the edges of the interval
        if yy(1) == Inf
            yy(1) = yy(2) + (yy(2)-yy(3));
        end
        if yy(end) == Inf
            yy(end) = yy(end-1) + (yy(end-1)-yy(end-2));
        end
        % Normalization
        yy = yy/trapz(xx,yy);
        % Plotting the curve
        plot(xx,yy,'r-','LineWidth',2,'DisplayName',"Target distribution")
    end
end

%% Plot the histogram of accepted leaves

% Extracting leaf information
leafCount = Leaves.leaf_count;
leafScale = Leaves.leaf_scale;
leafBaseArea = Leaves.base_area;
leafStartPoints = Leaves.leaf_start_point;

% Area of each leaf (leaf scaling identical in every dimension)
leafAreas = (leafScale(:,1).^2)*leafBaseArea;

% Heights of the leaves
hLeaves = leafStartPoints(:,3);

% Set the maximum and minimum height based on the point cloud
h_min = min(pCloud(:,3)); 
h_max = max(pCloud(:,3)); 
hLeaves(hLeaves<h_min) = h_min;
hLeaves(hLeaves>h_max) = h_max;

% Normalize leaf heights
relativeLeafHeights = (hLeaves-h_min)/(h_max-h_min);

% Calculate weighted histogram for the leaf area wrt. height
leafHist = zeros(nBins,1);
for iLeaf = 1:leafCount
    for iBin = 1:nBins
        if relativeLeafHeights(iLeaf) > binEdges(iBin) && ...
            relativeLeafHeights(iLeaf) <= binEdges(iBin+1)
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
                'DisplayName','Accepted leaf area','flipxy',0)
xlabel("relative height")
ylabel("leaf area frequency density [m^2]")
axis tight
legend('Location','southeast')