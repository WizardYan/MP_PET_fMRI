%% script to get violin plot data - can alternatively use Cornblath functions
% this data goes with violinplots_sps.R to make Figure 3 b/c i-iii

clear all; close all;clc
basedir = './';
cd(basedir);

%% set inputs

split='main'
load(fullfile(['./data/',split,'.mat']))
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory


% MP_stop = nscans;
% PL_start = MP_stop+1;
% tot = nscans*2;
MP_stop=37;
PL_start=MP_stop+1;
tot = 74;

%% Because of scrubbing some time points out, Weizheng Modified the scripts for calcualting the 
% number of clusters percan
subjInd(18772:end) = subjInd(18772:end) + nsubjs;


%% big loop
for numClusters=[6]
    load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
    %% count number of each cluster per scan
%     A=reshape(partition,TR,[]); 
    count = zeros(numClusters,tot);
    for b=1:tot
        index_for_specific_subject = find(subjInd==b);
        [count(:,b),~] = hist(partition(index_for_specific_subject),1:numClusters);
    end
    
    %% Calculate Fractional Occupancy

    MP=count(:,1:MP_stop); 
    MPfo = MP./sum(MP,1); % data for all scans
    
 
    PL=count(:,PL_start:tot);
    PLfo = PL./sum(PL,1); 
    
    %% Calculate Dwell Time and Appearance Rate
    dwell = []; %zeros(numClusters,220,89);
    sumApp=0;
    for c=1:numClusters
        for b=1:tot

            index_for_specific_subject = find(subjInd==b);

            appear=find(partition(index_for_specific_subject)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    dwell(c,i,b)=s;
                    a=a+1;
                else
                    dwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(dwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    
    MP = dwell(:,:,1:MP_stop); 
    MP_count=zeros(numClusters,3,nscans); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    MPdt1=[];
    MPar1=[];
    [~,maxIndex,~]=size(dwell);
    for c=1:numClusters
        for b=1:nsubjs
            MP_count(c,1,b)=MP_count(c,1,b)+sum(MP(c,:,b));
            a=1;
            while a<=maxIndex
                if MP(c,a,b)~= 0
                    MP_count(c,2,b)=MP_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            MPdt(c,b)=MP_count(c,1,b)/MP_count(c,2,b)*rep_time; %total time/#appear *2 to convert to seconds
            
            scan_length = sum(subjInd==b)*rep_time/60;
            MPar(c,b)=MP_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
    

    
    
    PL=dwell(:,:,PL_start:tot);
    PL_count=zeros(numClusters,3,nscans); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    PLdt1=[];
    PLar1=[];
    [~,maxIndex,~]=size(dwell);
    for c=1:numClusters
        for b=1:nsubjs
            PL_count(c,1,b)=PL_count(c,1,b)+sum(PL(c,:,b));
            a=1;
            while a<=maxIndex
                if PL(c,a,b)~= 0
                    PL_count(c,2,b)=PL_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            PLdt(c,b)=PL_count(c,1,b)/PL_count(c,2,b)*rep_time; %total time/#appear *0.891 to convert to seconds
            scan_length = sum(subjInd==b)*rep_time/60;
            PLar(c,b)=PL_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
       
 
    %compute DT's and AR's averaged over clusters (does MP overall have
    %lower DT's and higher AR's?)
    DT=NaN(1,nsubjs*2);
    DT(1:nsubjs)=mean(MPdt,1);
    DT(nsubjs+1:nsubjs*2)=mean(PLdt,1);
    
    [~,pDT,~,tDT] = ttest(DT(1:nsubjs),DT(nsubjs+1:nsubjs*2));
    
    AR=NaN(1,nsubjs*2);
    AR(1:nsubjs)=mean(MPar,1);
    AR(nsubjs+1:nsubjs*2)=mean(PLar,1);
    
    [~,pAR,~,tAR] = ttest(AR(1:nsubjs),AR(nsubjs+1:nsubjs*2));
    
    
    
    %% Save
    
    clusters=char(clusterNames);
    save(fullfile(savedir,['ViolinData_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'MPfo', 'PLfo', 'MPdt', 'PLdt', 'MPar', 'PLar', 'clusters','DT','AR');%,'pDT','pAR','tDT','tAR')

end
