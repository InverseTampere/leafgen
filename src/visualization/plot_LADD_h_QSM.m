function f = plot_LADD_h_QSM(QSM,Leaves,TargetLADD,varargin)

% Initialize values
nBins = 10;

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
dType = TargetLADD.dTypeLADDh;
if ~any(strcmp(dType,{'none','uniform','polynomial','polynomialmixture',...
        'weibull','weibullmixture','beta','betamixture'}))
    error("LADD height distribution type not recognized.")
end

% Read target distribution parameters
p = TargetLADD.hParams;

% Functions
fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_weibull = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);

% Relative height discretization
xx = 0:0.001:1;

% Bins for the histogram of accepted leaves
binEdges = linspace(0,1,nBins+1);

%% Plot the target distribution function

if dType ~= "none"
    switch dType
        case 'uniform'
            yy = ones(size(xx));
        case  'polynomial'
            yy = polyval(p,xx);
        case 'polynomialmixture'
            nP = (length(p)-1)/2; % number of polynomial coefficients
            p1 = p(1:nP); % coefficients of the first polynomial
            p2 = p((nP+1):(2*nP)); % coefficients of the second polynomial
            w = p(end); % mixture model weight
            yy = w*polyval(p1,xx) + (1-w)*polyval(p2,xx);
        case 'weibull'
            l = p(1); % scale parameter
            k = p(2); % shape parameter
            yy = fun_weibull(xx,l,k);
        case 'weibullmixture'
            l1 = p(1); k1 = p(2); % parameters of the first distribution
            l2 = p(3); k2 = p(4); % parameters of the second distribution
            w = p(5); % mixture model weight
            yy = w*fun_weibull(xx,l1,k1) + (1-w)*fun_weibull(xx,l2,k2);
        case 'beta'
            a = p(1);
            b = p(2);
            yy = fun_beta(xx,a,b);
        case 'betamixture'
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
    plot(xx,yy,'r:','LineWidth',2,'DisplayName',"Target distribution")
end

%% Plot the histogram of accepted leaves

% Extracting leaf information
leafCount = Leaves.leaf_count;
leafScale = Leaves.leaf_scale;
leafBaseArea = Leaves.base_area;
leafStartPoints = Leaves.leaf_start_point;
cylinderMidPoints = QSM.cylinder_mid_point;

% Area of each leaf (leaf scaling identical in every dimension)
leafAreas = (leafScale(:,1).^2)*leafBaseArea;

% Heights of the leaves and cylinders
hLeaves = leafStartPoints(:,3);
hCylinders = cylinderMidPoints(:,3);

% Set the maximum and minimum height based on the cylinder heights
h_min = min(hCylinders); 
h_max = max(hCylinders); 
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
custom_bar_plot(binEdges,leafHistFD,'FaceColor','b','FaceAlpha',0.3,...
                'DisplayName','Accepted leaf area','flipxy',0)
xlabel("relative height")
ylabel("leaf area frequency density [m^2]")
axis tight
legend('Location','northwest')