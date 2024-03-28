function fMax = weibull_mm_upper_limit(p)
% Gives maximum value of mixture model of Weibull distributions limited to 
% value 5 on closed  interval [0,1]. Note that might sometimes find only a
% local maximum.
fun_weibull = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k);
l1 = p(1); k1 = p(2); % parameters of the first distribution
l2 = p(3); k2 = p(4); % parameters of the second distribution
w = p(5); % mixture model weight
fun_mm = @(x) w*fun_weibull(x,l1,k1) + (1-w)*fun_weibull(x,l2,k2);
fun_mm_neg = @(x) -(w*fun_weibull(x,l1,k1)+(1-w)*fun_weibull(x,l2,k2));
if k1 <= 1 || k2 <= 1
    fMax = min([fun_mm(0.001),5]); % infinite node at 0
else
    % Find possible locations of maximum
    xx = [fminbnd(fun_mm_neg,0,1), l1*((k1-1)/k1).^(1/k1), ...
          l2*((k2-1)/k2).^(1/k2), 0, 1];
    xx = xx(xx<=1);
    fMax = max(fun_mm(xx));
end  