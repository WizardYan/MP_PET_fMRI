%% make radar plots
clear all; clc
basedir = './';
cd(basedir);
addpath(genpath('from_github'));
%% set inputs
numClusters = 6;
split='main'

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

%% Dendrogram
Z = linkage(centroids', 'ward','correlation'); % Using Ward's method
figure;
dendrogram(Z);
title('Dendrogram of k-means Centroids, method: ward, metric: correlation');
xlabel('Cluster Index');
ylabel('Distance');
            
            
%%
overallNames = clusterNames;
[nparc,numClusters] = size(centroids);
[~,~,~,net10angle] = NAME_CLUSTERS_ANGLE(centroids);

%% plot

YeoColors = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0];
YeoColors = [YeoColors;YeoColors];

[~,~,net10angle_Up,net10angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids);
YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN','SUB'};

numNets = numel(YeoNetNames);


%% make radial plots

clusterColors = GET_CLUSTER_COLORS(numClusters);

clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; 

f=figure;
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[net10angle_Up(K,:) net10angle_Up(K,1)],'k');
    polarplot(netAngle,[net10angle_Down(K,:) net10angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rlim([0.0 0.7]);
    rticks([0.2 0.4 0.8]); rticklabels({'','0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',8);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',10);
end

f.PaperUnits = 'inches';
f.PaperSize = [8, 1.5];
f.PaperPosition = [0 0 8 1.5];

%% This part has not completed, still need to save the figure mannually.  - Weizheng
saveas(f, fullfile(savedir,['systems_plot_sps_k',num2str(numClusters),'.pdf']))
