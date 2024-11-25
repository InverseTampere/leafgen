function binEdges = generate_bin_edges(interval,nBins,fun,binPartition)

if nargin == 4
    switch lower(binPartition)
        case 'probabilitymass'
            k_max = 10000;
            xx = linspace(interval(1),interval(2),k_max);
            yy = fun(xx);
            if yy(1) == Inf
                yy(1) = yy(2) + (yy(2)-yy(3));
            end
            if yy(end) == Inf
                yy(end) = yy(end-1) + (yy(end-1)-yy(end-2));
            end
            totalArea = trapz(xx,yy);
            binArea = totalArea/nBins;
            binEdges = zeros(nBins+1,1);
            binEdges(1) = interval(1);
            binEdges(end) = interval(2);
            k = 2;
            prevEdgeIndex = 1;
            for i = 2:nBins
                tempArea = 0;
                while tempArea < binArea
                    tempArea = trapz(xx(prevEdgeIndex:k), ...
                                     yy(prevEdgeIndex:k));
                    k = k + 1;
                end
                binEdges(i) = xx(k);
                prevEdgeIndex = k;
                k = k + 1;
            end
    end
else
    binEdges = linspace(interval(1),interval(2),nBins+1);
    binEdges = binEdges(:);
end

end