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

function fMax = vonmises_mm_upper_limit(p)
% Gives maximum value of mixture model of Von Mises distributions  on 
% closed  interval [0,2*pi]. Note that might sometimes find only a local
% maximum.
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));
m1 = p(1); k1 = p(2); % parameters of the first distribution
m2 = p(3); k2 = p(4); % parameters of the second distribution
w = p(5); % mixture model weight
fun_mm = @(x) w*fun_vonmises(x,m1,k1) + (1-w)*fun_vonmises(x,m2,k2);
fun_mm_neg = @(x) -(w*fun_vonmises(x,m1,k1) + (1-w)*fun_vonmises(x,m2,k2));
fMax = max([fun_mm(0), fun_mm(2*pi), fun_mm(fminbnd(fun_mm_neg,0,2*pi)),...
            fun_mm(m1), fun_mm(m2)]);
end