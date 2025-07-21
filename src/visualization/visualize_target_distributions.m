function f = visualize_target_distributions(TargetDistributions, ...
                                            structuralVariablePoint)

% Structural variable values
rh = structuralVariablePoint(1); % relative height
rd = structuralVariablePoint(2); % relative distance
cd = structuralVariablePoint(3); % compass direction

% Check which fields exist in the input struct
flagLADDh  = isfield(TargetDistributions,'dTypeLADDh');
flagLADDd  = isfield(TargetDistributions,'dTypeLADDd');
flagLADDc  = isfield(TargetDistributions,'dTypeLADDc');
flagLODinc = isfield(TargetDistributions,'dTypeLODinc');
flagLODaz  = isfield(TargetDistributions,'dTypeLODaz');
flagLSD    = isfield(TargetDistributions,'dTypeLSD');

% Function definitions
fun_beta     = @(x,a,b) (1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
fun_weibull  = @(x,l,k) (k/l)*(x/l).^(k-1).*exp(-(x/l).^k) ...
                        /(1-exp(-(1/l)^k));
fun_vonmises = @(x,m,k) exp(k*cos(x-m))./(2*pi*besseli(0,k));
fun_dewit    = @(x,a,b) (1 + a*cos(b*x))/(pi/2+(a/b)*sin(b*pi/2));
fun_normal   = @(x,m,v) exp(-(x-m).^2/(2*v))/sqrt(2*pi*v);

%% Assign function handles to be plotted

% LADD
if flagLADDh
    dTypeLADDh = TargetDistributions.dTypeLADDh;
    pLADDh     = TargetDistributions.pLADDh;
    switch dTypeLADDh
        case 'uniform'
            fDistLADDh = @(h) 1*ones(size(h));
        case 'polynomial'
            fDistLADDh = @(h) polyval(pLADDh,h);
        case 'polynomialmixture'
            nP = (length(pLADDh)-1)/2; % number of polynomial coefficients
            p1 = pLADDh(1:nP); % coefficients of the first polynomial
            p2 = pLADDh((nP+1):(2*nP)); % coefficients of the second polyn.
            w = pLADDh(end); % mixture model weight
            fDistLADDh = @(h) w*polyval(p1,h) + (1-w)*polyval(p2,h);
        case 'weibull'
            l = pLADDh(1); % scale parameter
            k = pLADDh(2); % shape parameter
            fDistLADDh = @(h) fun_weibull(h,l,k);
        case 'weibullmixture'
            l1 = pLADDh(1); k1 = pLADDh(2); % parameters of the first dist.
            l2 = pLADDh(3); k2 = pLADDh(4); % parameters of the second dist.
            w = pLADDh(5); % mixture model weight
            fDistLADDh = @(h) w*fun_weibull(h,l1,k1) ...
                              + (1-w)*fun_weibull(h,l2,k2);
        case 'beta'
            a = pLADDh(1);
            b = pLADDh(2);
            fDistLADDh = @(h) fun_beta(h,a,b);
        case 'betamixture'
            a1 = pLADDh(1); b1 = pLADDh(2); % parameters of the first dist.
            a2 = pLADDh(3); b2 = pLADDh(4); % parameters of the second dist.
            w = pLADDh(5); % mixture model weight
            fDistLADDh = @(h) w*fun_beta(h,a1,b1) + (1-w)*fun_beta(h,a2,b2);
    end
end

if flagLADDd
    dTypeLADDd = TargetDistributions.dTypeLADDd;
    pLADDd     = TargetDistributions.pLADDd;
    switch dTypeLADDd
        case 'uniform'
            fDistLADDd = @(d) 1*ones(size(d));
        case 'polynomial'
            fDistLADDd = @(d) polyval(pLADDd,d);
        case 'polynomialmixture'
            nP = (length(pLADDd)-1)/2; % order of polynomial
            p1 = pLADDd(1:nP); % coefficients of the first polynomial
            p2 = pLADDd((nP+1):(2*nP)); % coefficients of the second polyn.
            w = pLADDd(end); % mixture model weight
            fDistLADDd = @(d) w*polyval(p1,d) + (1-w)*polyval(p2,d);
        case 'weibull'
            l = pLADDd(1); % scale parameter
            k = pLADDd(2); % shape parameter
            fDistLADDd = @(d) fun_weibull(d,l,k);
        case 'weibullmixture'
            l1 = pLADDd(1); k1 = pLADDd(2); % parameters of the first dist.
            l2 = pLADDd(3); k2 = pLADDd(4); % parameters of the second dist.
            w = pLADDd(5); % mixture model weight
            fDistLADDd = @(d) w*fun_weibull(d,l1,k1) ...
                              + (1-w)*fun_weibull(d,l2,k2);
        case 'beta'
            a = pLADDd(1);
            b = pLADDd(2);
            fDistLADDd = @(d) fun_beta(d,a,b);
        case 'betamixture'
            a1 = pLADDd(1); b1 = pLADDd(2); % parameters of the first dist.
            a2 = pLADDd(3); b2 = pLADDd(4); % parameters of the second dist.
            w = pLADDd(5); % mixture model weight
            fDistLADDd = @(d) w*fun_beta(d,a1,b1) + (1-w)*fun_beta(d,a2,b2);
    end
end

if flagLADDc
    dTypeLADDc = TargetDistributions.dTypeLADDc;
    pLADDc     = TargetDistributions.pLADDc;
    switch dTypeLADDc
        case 'uniform'
            fDistLADDc = @(c) 1/(2*pi)*ones(size(c));
        case 'vonmises'
            m = pLADDc(1); % mean
            k = pLADDc(2); % measure of concentration
            fDistLADDc = @(c) fun_vonmises(c,m,k);
        case 'vonmisesmixture'
            m1 = pLADDc(1); k1 = pLADDc(2); % parameters of the first dist.
            m2 = pLADDc(3); k2 = pLADDc(4); % parameters of the second dist.
            w = pLADDc(5); % mixture model weight
            fDistLADDc = @(c) w*fun_vonmises(c,m1,k1) ...
                              + (1-w)*fun_vonmises(c,m2,k2);
    end
end

% LOD
if flagLODinc
    dTypeLODinc = TargetDistributions.dTypeLODinc;
    pLODinc     = TargetDistributions.fun_pLODinc(rh,rd,cd);
    switch dTypeLODinc
        case 'uniform'
            % Uniform distribution
            fDistLODinc = @(x) 2/pi*ones(size(x));
        case 'spherical'
            % Spherical distribution function
            fDistLODinc = @(x) sin(x);
        case 'dewit'
            % Generalized de Wit's distribution function
            a = pLODinc(1);
            b = pLODinc(2);
            fDistLODinc = @(x) fun_dewit(x,a,b);
        case 'beta'
            % Beta distribution density function
            a = pLODinc(1);
            b = pLODinc(2);
            fDistLODinc = @(x) fun_beta(2*x/pi,a,b);
        case 'constant'
            % SPECIAL CASE
    end
end

if flagLODaz
    dTypeLODaz = TargetDistributions.dTypeLODaz;
    pLODaz     = TargetDistributions.fun_pLODaz(rh,rd,cd);
    switch dTypeLODaz
        case 'uniform'
            % Uniform distribution
            fDistLODaz = @(x) 1/(2*pi)*ones(size(x));
        case 'vonmises'
            % Von Mises distribution density function
            m = pLODaz(1);
            k = pLODaz(2);
            fDistLODaz = @(x) fun_vonmises(x,m,k);
        case 'constant'
            % SPECIAL CASE
    end
end

% LSD
if flagLSD
    dTypeLSD = TargetDistributions.dTypeLSD;
    pLSD     = TargetDistributions.fun_pLSD(rh,rd,cd);
    switch dTypeLSD
        case 'uniform'
            % Uniform distribution on the interval [a,b]
            a = pLSD(1);
            b = pLSD(2);
            fDistLSD = @(x) ones(size(x))/(b-a);
        case 'normal'
            % Normal distribution with mean m and variance v
            m = pLSD(1);
            v = pLSD(2);
            fDistLSD = @(x) fun_normal(x,m,v);
        case 'constant'
            % SPECIAL CASE
    end
end

%% Plotting the functions

% Figure initialization
f = figure;
fPos = get(f,"Position");
newPos = [fPos(1)-0.5*fPos(3) fPos(2) 1.5*fPos(3) fPos(4)];
set(f,'Position',newPos);
nLADD = sum([flagLADDh flagLADDd flagLADDc]);
nLOD  = sum([flagLODinc flagLODaz]);
nLSD  = sum(flagLSD);
switch nLADD + nLOD + nLSD
    case 1
        rows = 1;
        cols = 1;
    case 2
        rows = 2;
        cols = 1;
    case 3
        rows = 2;
        cols = 2;
    case 4
        rows = 2;
        cols = 2;
    case 5
        rows = 2;
        cols = 3;
    case 6
        rows = 2;
        cols = 3;
end
t = tiledlayout(rows,cols,'TileSpacing','compact');
title(t,"Target leaf distributions at point [" + num2str(rh) + " " ...
      + num2str(rd) + " " + num2str(cd) + "]")

% LADD
if nLADD > 0
    plottedLADD = [false false false];
    for ii = 1:nLADD
        nexttile
        if flagLADDh == true && plottedLADD(1) == false
            xx = linspace(0,1,200);
            yy = fDistLADDh(xx);
            plot(xx,yy,'-','LineWidth',2,'Color',"blue")
            axis([0 1 0 1.1*max(yy)])
            title("LADD relative height")
            xlabel("relative height")
            plottedLADD(1) = true;
        elseif flagLADDd == true && plottedLADD(2) == false
            xx = linspace(0,1,200);
            yy = fDistLADDd(xx);
            plot(xx,yy,'-','LineWidth',2,'Color',"blue")
            axis([0 1 0 1.1*max(yy)])
            title("LADD relative distance")
            xlabel("relative distance")
            plottedLADD(2) = true;
        elseif flagLADDc == true && plottedLADD(3) == false
            xx = linspace(0,2*pi,200);
            yy = fDistLADDc(xx);
            plot(xx,yy,'-','LineWidth',2,'Color',"blue")
            axis([0 2*pi 0 1.1*max(yy)])
            title("LADD compass direction")
            xlabel("compass direction [rad]")
            xticks([0 pi/2 pi 3*pi/2 2*pi])
            xticklabels(["0" "\pi/2" "\pi" "3\pi/2" "2\pi"])
            plottedLADD(3) = true;
        end
    end
end

% LOD
if nLOD > 0
    plottedLOD = [false false];
    for jj = 1:nLOD
        nexttile
        if flagLODinc == true && plottedLOD(1) == false
            if dTypeLODinc == "constant"
                xx = [0 pLODinc pLODinc pLODinc pi/2];
                yy = [0 0 1 0 0];
            else
                xx = linspace(0,pi/2,200);
                yy = fDistLODinc(xx);
            end
            plot(xx,yy,'-','LineWidth',2,'Color',"#9400D3")
            axis([0 pi/2 0 1.1*max(yy)])
            title("LOD inclination angle")
            xlabel("inclination angle [rad]")
            xticks([0 pi/8 pi/4 3*pi/8 pi/2])
            xticklabels(["0" "\pi/8" "\pi/4" "3\pi/8" "\pi/2"])
            plottedLOD(1) = true;
        elseif flagLODaz == true && plottedLOD(2) == false
            if dTypeLODaz == "constant"
                xx = [0 pLODaz pLODaz pLODaz 2*pi];
                yy = [0 0 1 0 0];
            else
                xx = linspace(0,2*pi,200);
                yy = fDistLODaz(xx);
            end
            plot(xx,yy,'-','LineWidth',2,'Color',"#9400D3")
            axis([0 2*pi 0 1.1*max(yy)])
            title("LOD azimuth angle")
            xlabel("azimuth angle [rad]")
            xticks([0 pi/2 pi 3*pi/2 2*pi])
            xticklabels(["0" "\pi/2" "\pi" "3\pi/2" "2\pi"])
            plottedLOD(2) = true;
        end
    end
end

% LSD
if nLSD > 0
    nexttile
    switch dTypeLSD
        case 'uniform'
            lb = pLSD(1);
            ub = pLSD(2);
        case 'normal'
            m = pLSD(1);
            v = pLSD(2);
            lb = m - 4*sqrt(v);
            ub = m + 4*sqrt(v);
        case 'constant'
            lb = 0.9*pLSD;
            ub = 1.1*pLSD;
    end
    if dTypeLSD == "constant"
        xx = [lb pLSD pLSD pLSD ub];
        yy = [0 0 1 0 0];
    else
        xx = linspace(lb,ub,200);
        yy = fDistLSD(xx);
    end
    plot(xx,yy,'-','LineWidth',2,'Color',"#00CD94")
    axis([lb ub 0 1.1*max(yy)])
    title("LSD")
    xlabel("leaf area [m^2]")
end
end