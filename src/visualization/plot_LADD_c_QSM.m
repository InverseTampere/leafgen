function f = plot_LADD_c_QSM(QSM,Leaves,TargetLADD,varargin)

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
dType = TargetLADD.dTypeLADDc;
if ~any(strcmp(dType,{'none','uniform','vonmises','vonmisesmixture'}))
    error("LADD compass direction distribution type not recognized.")
end

% Read target distribution parameters
p = TargetLADD.pLADDc;

% Functions
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));

% Compass direction discretization
xx = 0:0.001:2*pi;

% Bins for the histogram of accepted leaves
binEdges = linspace(0,2*pi,nBins+1);

%% Plot the target distribution function
if dType ~= "none"
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

%% Histogram based on accepted leaves

% % Extracting QSM and leaf information
leafScale      = Leaves.leaf_scale;
leafBaseArea   = Leaves.base_area;
leafStartPoint = Leaves.leaf_start_point;
leafCount      = Leaves.leaf_count;
cylinderStartPoint = QSM.cylinder_start_point;
cylinderAxis       = QSM.cylinder_axis;
cylinderLength     = QSM.cylinder_length;

% Maximum cylinder length
maxLen = max(cylinderLength);
% Normalized cylinder lengths
cylLenNorm = cylinderLength./maxLen;
% Cylinder midpoints
midPoints = cylinderStartPoint + 0.5*cylinderLength.*cylinderAxis;
% Set the mean value of cylinder locations on xy-plane weighted with 
% cylinder length as the origin
xyOrigin = sum(cylLenNorm.*midPoints(:,1:2),1)./sum(cylLenNorm);

% Area of each leaf (leaf scaling identical in every dimension)
leafAreas = (leafScale(:,1).^2)*leafBaseArea;

% Coordinates of the leaf start points on xy-plane
leafCoord = leafStartPoint(:,1:2);
% Unit vectors pointing the direction of cylinder start point coordinates
leafDirection = zeros(leafCount,2);
for j = 1:leafCount
    temp = leafCoord(j,:) - xyOrigin;
    leafDirection(j,:) = temp/norm(temp);
end
% North set as the zero angle
northDirection = [0 1];
% Angles of the cylinder direction unit vectors w.r.t north
phiLeaf = zeros(leafCount,1);
for j = 1:leafCount
    if leafDirection(j,1) <= 0
        phiLeaf(j) = acos(dot(northDirection,leafDirection(j,:)));
    else
        phiLeaf(j) = 2*pi - acos(dot(northDirection,leafDirection(j,:)));
    end
end

% Calculate weighted histogram for the leaf area wrt. direction
leafHist = zeros(nBins,1);
for iLeaf = 1:leafCount
    leafAngle = phiLeaf(iLeaf);
    for iBin = 1:nBins
        if leafAngle > binEdges(iBin) && leafAngle <= binEdges(iBin+1)
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
xlabel("direction")
ylabel("leaf area frequency density [m^2]")
axis tight
legend('Location','northwest')

end