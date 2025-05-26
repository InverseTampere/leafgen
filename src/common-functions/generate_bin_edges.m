% This file is part of LeafGen
% 
% LeafGen is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% LeafGen is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with LeafGen.  If not, see <https://www.gnu.org/licenses/>.

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