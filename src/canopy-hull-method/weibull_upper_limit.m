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

function fMax = weibull_upper_limit(l,k)
% Gives maximum value of Weibull distribution limited to value 5 on closed 
% interval [0,1]
fun_weibull = @(x) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);
if k <= 1
    fMax = min([fun_weibull(0.001),5]); % infinite mode at 0
elseif l*((k-1)/k).^(1/k) < 1
    fMax = fun_weibull(l*((k-1)/k).^(1/k)); % mode between 0 and 1
else
    fMax = fun_weibull(1); % mode at 1
end