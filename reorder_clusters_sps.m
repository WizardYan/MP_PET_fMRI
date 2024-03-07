%% This script may be used to combine cluster identities in a partition, 
% also (more commonly) can be used to reorder and rename clusters in meta-state hierarchy
% there is a function in the Cornblath repo (reorderClusters.m) to do this automatically, however i did not investigate using it

% use max correlation to a set of centroids or visualize with
% systems_plot.sps first to determin new order and manually input order
% below

clear all; clc
basedir = '/Users/yanw4/Documents/LNI_P1_MP_DFNC_202212';
cd(basedir);

%% set inputs
numClusters = 6;
split='main' %must match file you want to load

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

[X,~]=size(partition);

%% Reassign clusters in partition 

clusterNames={'FPN-';'FPN+';'SOM-';'SOM+';'VIS-';'VIS+'}; %new cluster names in the new order


c_partition = NaN(X,1); %initialize new partition 

c_partition(partition == 1) = 3;
c_partition(partition == 2) = 4;
c_partition(partition == 3) = 5;
c_partition(partition == 4) = 1;
c_partition(partition == 5) = 6; 
c_partition(partition == 6) = 2; 

partition = c_partition; %rename to original naming for streamlined analysis.

%% reorder centroids (skip this section if combining instead of re-ordering)
[nparc,~]=size(centroids);
r_centroids = NaN(nparc,numClusters);

% Check to make sure the new order here matches the one above! (LH / RH
% sides of assignment equations should swap)
r_centroids(:,1) = centroids(:,4);
r_centroids(:,2) = centroids(:,6);
r_centroids(:,3) = centroids(:,1);
r_centroids(:,4) = centroids(:,2);
r_centroids(:,5) = centroids(:,3);
r_centroids(:,6) = centroids(:,5);

centroids = r_centroids;

%% save partition in separate file if combining (add prefix 'c' to filename) or overwrite old file if just re-ordering

save(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids');
