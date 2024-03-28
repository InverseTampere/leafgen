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