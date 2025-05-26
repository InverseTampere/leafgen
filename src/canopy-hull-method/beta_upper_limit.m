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

function fMax = beta_upper_limit(a,b)
% Gives maximum value of beta distribution limited to value 10 on closed 
% interval [0,1]
fun_beta = @(x) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
if a == b
    fMax = fun_beta(0.5);
elseif a > 1 && b > 1
    fMax = fun_beta((a-1)/(a+b-2));
elseif a < 1 && b < 1
    fMax = 10; % bimodal with infinite modes
elseif a == 1 && b > 1
    fMax = fun_beta(0);
elseif a > 1 && b == 1
    fMax = fun_beta(1);
else
    fMax = 10; % infinite mode at 0 or 1
end
fMax = min(fMax,10); % limit the maximum value to 10
end
