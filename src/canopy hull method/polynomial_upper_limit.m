function fMax = polynomial_upper_limit(p)
% Gives maximum value of polynomial on closed interval [0,1]
derParams = polyder(p);
derRoots = roots(derParams);
derRootsOnInterval = derRoots(derRoots>=0);
derRootsOnInterval = derRootsOnInterval(derRootsOnInterval<=1);
criticalPoints = [0; 1; derRootsOnInterval];
fMax = max(polyval(p,criticalPoints));
end