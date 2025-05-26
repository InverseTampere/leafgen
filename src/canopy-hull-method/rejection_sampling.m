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

function x = rejection_sampling(fun,maxValue)
accepted = 0;
while accepted == 0
    % Proposal value
    proposal = rand(1);
    % Function value on propsal point
    funValue = fun(proposal);
    vertValue = rand(1)*maxValue;
    if vertValue < funValue
        x = proposal;
        accepted = 1;
    end
end
end