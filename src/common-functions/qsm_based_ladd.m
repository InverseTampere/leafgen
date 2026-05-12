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

function relAreas = qsm_based_ladd(CylinderParameters, ...
                                   TargetDistributions, ...
                                   leaflessBranchInds)

% Scaling function for the relative distance along subbranch (here the 
% default x^4 is arbitrarily chosen)
if isempty(TargetDistributions.pLADDd)
    f_dist_scaling = @(x) x.^4;
else
    p = TargetDistributions.pLADDd;
    f_dist_scaling = @(x) x.^p;
end

% Set cylinder leaf area budgets based on the scaling function
cylAreas = f_dist_scaling( ...
               CylinderParameters.relative_distance_along_subbranch);

% Scale the areas with cylinder length
cylAreas = cylAreas.*CylinderParameters.length;

% Prevent leaf generation on the tree stem (or other branches set leafless)
for i = leaflessBranchInds
    cylAreas(CylinderParameters.branch_index == i) = 0;
end

% Normalize to get the relative areas
totalArea = sum(cylAreas);
if totalArea > 0
    relAreas = cylAreas./totalArea;
else
    % All cylinders are leafless
    relAreas = 0;
end

end