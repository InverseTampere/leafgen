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

function relAreas = qsm_based_ladd(CylinderParameters)

% Scaling function for the relative distance along subbranch (here the 
% x^4 is arbitrarily chosen, feel free to modify this)
f_dist_scaling = @(x) x.^4;

% Set cylinder leaf area budgets based on the scaling function
cylAreas = f_dist_scaling( ...
               CylinderParameters.relative_distance_along_subbranch);

% Scale the areas with cylinder length
cylAreas = cylAreas.*CylinderParameters.length;

% Prevent leaf generation on the tree stem
cylAreas(CylinderParameters.branch_index == 1) = 0;

% Normalize to get the relative areas
relAreas = cylAreas./sum(cylAreas);

end