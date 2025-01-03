%  make elbow plot (Fig S2a-b) of within cluster / total variance explained by clusters at different k values
% assumes correlation distance

%% init

%% SI Figure 13 (need to run repeatkmeans_sps.m first)

clear all; close all;clc
basedir = '/Users/yanw4/Documents/LNI_P1_MP_DFNC_202212';
cd(basedir);

split = 'main'

savedir = fullfile(basedir,'results','example');mkdir(savedir);
load(fullfile(basedir,['data/',split,'.mat'])) % make sure this file matches the split you want to run below

%% inputs
N = size(concTS,1); % number of observations
maxk = round(sqrt(TR))-1; %<- check, may not need "-1" maxk^2 should be < total TR in scan

k_rng = 2:maxk;

VarianceExplained = zeros(length(k_rng),1);

for numClusters = k_rng
	disp(['K = ',num2str(numClusters)])
	load(fullfile(basedir,['/repkmeans/kmeans',num2str(split),'k_',num2str(numClusters),'.mat']));
	kClusterCentroids = GET_CENTROIDS(concTS,partition,numClusters);
	VarianceExplained(numClusters - 1) = VAREXPLAINED(concTS,partition,kClusterCentroids,numClusters);
end

save(fullfile(savedir,['VarianceExplained',num2str(split),'.mat']),'VarianceExplained','k_rng');

% Supplimentary figure
f=figure;
plot(k_rng,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['ElbowPlot',num2str(split),'.pdf']));

f = figure;
plot(k_rng(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['GainInVarianceExplained',num2str(split),'.pdf']));
