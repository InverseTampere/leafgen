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

function fMax = beta_mm_upper_limit(p)
% Gives maximum value of mixture model of beta distributions limited to 
% value 10 on closed  interval [0,1]. Note that might sometimes find only a
% local maximum.
fun_beta = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
a1 = p(1); b1 = p(2); % parameters of the first distribution
a2 = p(3); b2 = p(4); % parameters of the second distribution
w = p(5); % mixture model weight
fun_mm = @(x) w*fun_beta(x,a1,b1) + (1-w)*fun_beta(x,a2,b2);
fun_mm_neg = @(x) -(w*fun_beta(x,a1,b1) + (1-w)*fun_beta(x,a2,b2));
% Find possible maximum values
if (a1<0 || a2<0) && (b1<0 || b2<0)
    fMax = min([max([fun_mm(0.001) fun_mm(0.999)]) 10]);
else
    [xMax,negSearchMax] = fminbnd(fun_mm_neg,0,1);
    if 0.001 < xMax && xMax < 0.999
        searchMax = -negSearchMax;
    else
        searchMax = 0;
    end
    leftMax = min([fun_mm(0.001) 10]);
    rightMax = min([fun_mm(0.999) 10]);
    maxProposals = [searchMax leftMax rightMax];
    if a1 > 1 && b1 > 1
        dist1Max = fun_mm((a1-1)/(a1+b1-2));
        maxProposals = [maxProposals dist1Max];
    end
    if a2 > 1 && b2 > 1
        dist2Max = fun_mm((a2-1)/(a2+b2-2));
        maxProposals = [maxProposals dist2Max];
    end
    fMax = max(maxProposals);
end