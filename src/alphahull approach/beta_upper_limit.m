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
