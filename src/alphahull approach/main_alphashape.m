clear, close all
addpath("C:\Users\bppimo\Documents\GitHub\foliage-generation-repo\src")
addpath("C:\Users\bppimo\Documents\GitHub\foliage-generation-repo\src\classes")

%% Point cloud

load('examplePC.mat')
pCloud = double(Small);

% Translate the lowest point of the cloud to the plane z=0
pCloud = pCloud - min(pCloud(:,3));

% Translate the trunk of the tree to origin
tfBottom = pCloud(:,3) < 1;
pCloud = pCloud - [mean(pCloud(tfBottom,1:2)) 0];

%% Initialize leaf base geometry

% Vertices of the leaf basis geometry
vertices = [-0.04  0.0   0.0;
            0.0    0.08  0.0;
            0.04   0.0   0.0];

% Triangles of the leaf basis geometry
tris = [1, 2, 3];

%% Define target leaf distributions

% LADD relative height
TargetDistributions.dType_h = 'beta';
TargetDistributions.p_h = [4.5 1];
TargetDistributions.nBins_h = 10;

% LADD relative distance from stem
TargetDistributions.dType_d = 'beta';
TargetDistributions.p_d = [2 1];
TargetDistributions.nBins_d = 10;

% LADD compass direction
TargetDistributions.dType_c = 'vonmises';
TargetDistributions.p_c = [pi 0.1];
TargetDistributions.nBins_c = 10;

% LOD inclination angle
TargetDistributions.dType_inc = 'dewit';
TargetDistributions.fun_inc_params = @(h,d,c) [1,2];
% TargetDistributions.p_inc = [1 2];

% LOD azimuth angle
TargetDistributions.dType_az = 'vonmises';
TargetDistributions.fun_az_params = @(h,d,c) [3.3, 0.25];

% LSD
TargetDistributions.dType_size = 'uniform';
TargetDistributions.fun_size_params = @(h,d,c) [0.0021, 0.0038];

%% Generate foliage

totalLeafArea = 50;

[Leaves,aShape] = generate_foliage_alphashape(pCloud,TargetDistributions, ...
                                     totalLeafArea,vertices,tris, ...
                                     'alpha',1);

%% Visualize the foliage

figure(1), clf

if 1
    % Plot alphashape
    tiledlayout(1,2)
    ax1 = nexttile;
    pc = aShape.Points;
    plot3(pc(:,1),pc(:,2),pc(:,3),'g.','MarkerSize',3)
    hold on, grid on
    plot(aShape,'FaceColor','m','FaceAlpha',0.2)
    pcTop = max(pc(:,3));
    plot3([0 0],[0 0],[0 pcTop],'c-','LineWidth',3)
    ax2 = nexttile;
    Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
end

% Plot leaves
hLeaf = Leaves.plot_leaves();
% Set leaf color
set(hLeaf,'FaceColor',[0,150,0]./255,'EdgeColor','none');
grid on
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

%% Plot leaf distributions
return
plot_LADD_h(Leaves,TargetDistributions);
plot_LADD_d(Leaves,TargetDistributions);
plot_LADD_c(Leaves,TargetDistributions);

% plot_LOD_inc(Leaves,fType_inc,p_inc,6)%,'HeightBins',3);
% plot_LOD_az(Leaves,fType_az,p_az,5);

% plot_LSD(Leaves,fType_size,p_size,10);