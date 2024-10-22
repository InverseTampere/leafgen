function [] = custom_bar_plot(binEdges,binValues,varargin)

faceColor = 'b';
faceAlpha = 0.5;
displayName = '';
flipxy = 0;

% Check additional parameters.
i = 1;
NArg = numel(varargin);
while i <= NArg
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'facecolor'
                faceColor = varargin{i+1};
                i = i + 1;
            case 'facealpha'
                faceAlpha = varargin{i+1};
                i = i + 1;
            case 'displayname'
                displayName = varargin{i+1};
                i = i + 1;
            case 'flipxy'
                flipxy = varargin{i+1};
                i = i + 1;
        end
    end
    i = i + 1;
end

% Number of bins
nBins = length(binValues);

% Define variables for the plot
xx = zeros(4*nBins,1);
yy = zeros(4*nBins,1);
for iBin = 1:nBins
    xx(4*(iBin-1)+1 : 4*(iBin-1)+2) = binEdges(iBin);
    xx(4*(iBin-1)+3 : 4*(iBin-1)+4) = binEdges(iBin+1);
    yy(4*(iBin-1)+2 : 4*(iBin-1)+3) = binValues(iBin);
end

% Plotting
area(xx,yy,'FaceColor',faceColor,'FaceAlpha',faceAlpha,'DisplayName',...
     displayName)

if flipxy == 1
    view([90 -90])
end