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

function fMax = polynomial_upper_limit(p)
% Gives maximum value of polynomial on closed interval [0,1]
derParams = polyder(p);
derRoots = roots(derParams);
derRootsOnInterval = derRoots(derRoots>=0);
derRootsOnInterval = derRootsOnInterval(derRootsOnInterval<=1);
criticalPoints = [0; 1; derRootsOnInterval];
fMax = max(polyval(p,criticalPoints));
end