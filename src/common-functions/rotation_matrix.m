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

function R = rotation_matrix(u,k)
% Matrix for rotation around u by k radians

    R = zeros(3,3);
    c = cos(k); 
    s = sin(k);
    R(1,:) = [u(1)^2+(1-u(1)^2)*c,    ...
              u(1)*u(2)*(1-c)-u(3)*s, ...
              u(1)*u(3)*(1-c)+u(2)*s];
    
    R(2,:) = [u(1)*u(2)*(1-c)+u(3)*s, ...
              u(2)^2+(1-u(2)^2)*c,    ...
              u(2)*u(3)*(1-c)-u(1)*s];
    
    R(3,:) = [u(1)*u(3)*(1-c)-u(2)*s, ...
              u(2)*u(3)*(1-c)+u(1)*s, ...
              u(3)^2+(1-u(3)^2)*c];

end