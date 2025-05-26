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

function [V,I] = find_closest_node(nodes,value)

% Find the node closest to the given value

if isscalar(nodes)
    I = 1;
elseif value > max(nodes)
    I = length(nodes);
else
    % First node above or equal to value
    iNodeAbove = find(nodes >= value,1,'first');
    if iNodeAbove == 1
        I = iNodeAbove;
    else
        % First node below value
        iNodeBelow = iNodeAbove - 1;
        % Choose the node closer to value
        if abs(value-nodes(iNodeBelow)) < abs(value-nodes(iNodeAbove))
            I = iNodeBelow;
        else
            I = iNodeAbove;
        end
    end
end
V = nodes(I);

end