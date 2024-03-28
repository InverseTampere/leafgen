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