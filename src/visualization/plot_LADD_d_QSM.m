function f = plot_LADD_d_QSM(QSM,Leaves,TargetLADD,varargin)

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
dType = TargetLADD.dTypeLADDd;
if ~any(strcmp(dType,{'none','uniform','polynomial','polynomialmixture',...
        'weibull','weibullmixture','beta','betamixture'}))
    error("LADD distance from stem distribution type not recognized.")
end

% Read target distribution parameters
p = TargetLADD.pLADDd;

% Functions
fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_weibull = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);

% Relative distance discretization
xx = 0:0.001:1;

% Bins for the histogram of accepted leaves
binEdges = linspace(0,1,nBins+1)';

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

%% Histogram based on accepted leaves

% Extracting QSM and leaf information
cylinderStartPoint  = QSM.cylinder_start_point;
cylinderAxis        = QSM.cylinder_axis;
cylinderBranchIndex = QSM.cylinder_branch_index;
cylinderLength      = QSM.cylinder_length;
cylinderCount       = QSM.block_count;
petioleStartPoint = Leaves.petiole_start_point;
leafParents       = Leaves.leaf_parent;
leafCount         = Leaves.leaf_count;
leafScale         = Leaves.leaf_scale;
leafBaseArea      = Leaves.base_area;

% Area of each leaf (leaf scaling identical in every dimension)
leafAreas = (leafScale(:,1).^2)*leafBaseArea;

% Cylinder index vector
cylIndices = (1:1:cylinderCount)';

% Leaf index vector
leafIndices = (1:1:leafCount)';

%Leaf area histogram variable
leafHist = zeros(nBins,1);

for iBranch = 0:max(cylinderBranchIndex)
    % Indices of QSM cylinders belonging in branch
    bcIndices = cylIndices(cylinderBranchIndex == iBranch);
    if sum(bcIndices) == 0
        % the branch has no cylinders
        continue
    end
    % Lengths of corresponding cylinders
    bcLengths = zeros(length(bcIndices),1);
    for j = 1:length(bcLengths)
        bcLengths(j) = cylinderLength(bcIndices(j));
    end
    % Cumulative sum of the lengths
    bcLengthsCumulative = cumsum(bcLengths);

    % Branch cylinder edges relative to subbranch distance
    bcRelEdges = [0; bcLengthsCumulative]/bcLengthsCumulative(end);

    % Loop over branch cylinders
    for iBC = 1:length(bcIndices)
        % Find the indices of the leaves attached to the cylinder
        childLeafInds = leafIndices(leafParents == bcIndices(iBC));
        % If cylinder has no leaves, advance to the next cylinder
        if isempty(childLeafInds) == true
            continue
        end
        % Vetor from cylinder start to petiole start
        cylStartToPetioleStart = petioleStartPoint(childLeafInds,:) ...
                              - cylinderStartPoint(bcIndices(iBC),:);
        % Projection to cylinder axis
        cylAxis = cylinderAxis(bcIndices(iBC),:);
        petioleStartProj = cylStartToPetioleStart*cylAxis' ...
                       .*cylAxis./(norm(cylAxis).^2);
        % Leaf start point values as relative branch distance
        leafRelPos = sqrt(sum(petioleStartProj.^2,2)) ...
                     ./cylinderLength(bcIndices(iBC)) ...
                     *(bcRelEdges(iBC+1)-bcRelEdges(iBC)) ...
                     + bcRelEdges(iBC);
        % Collect leaf area to correct bins
        for iBin = 1:nBins
            indsInBin = childLeafInds(and(leafRelPos>binEdges(iBin), ...
                                          leafRelPos<=binEdges(iBin+1)));
            if isempty(indsInBin) == false
                leafHist(iBin) = leafHist(iBin)+sum(leafAreas(indsInBin));
            end
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
xlabel("relative sub-branch distance")
ylabel("leaf area frequency density [m^2]")
axis tight
legend('Location','northwest')

end