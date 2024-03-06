clc;clear
basedir = '/Users/yanw4/Documents/LNI_P1_MP_DFNC_202212/';

atlas_name = 'schae232'; % other options:'shen_atlas'
nparc = 232; % 200 cortical regions + 32 subcortical regions.
nscans = 74; % 37 subjects, 2 sessions per subject.
nsubjs = 37; % 37 subjects.
subj_scanInd = [1:nsubjs,1:nsubjs];
rep_time = 0.891;
% scan_length = 7.33
% split = 'main';

load([basedir,'results/demographic.mat']);
PL_time_series = [];
MP_time_series = [];

temp_MPsubjInd = [];
temp_PLsubjInd = [];

for i = 1:nsubjs
    
    if fMRI_session_PL_MP(i) == 1 
        session_for_PL = '1'
        session_for_MP = '2'
    else
        session_for_PL = '2'
        session_for_MP = '1'
    end

    %% Processing the placebo
    PL_filename_timeseries = strcat([basedir,'results/',atlas_name,'/tcs/sub-',char(control_finished_demogpraph_table(i,2).Variables),'_ses-',session_for_PL,'_atlas_schae232_roi_time_series_scrubbed.txt'])
    time_series = load (PL_filename_timeseries);
    time_series = time_series';

    time_series = time_series - mean(time_series, 1);

%     time_series_zscore = zscore(time_series);

    PL_time_series = [PL_time_series;time_series]; % concatenate the time series.
    temp_subjInd = i * ones(size(time_series,1),1);
    temp_PLsubjInd = [temp_PLsubjInd;temp_subjInd]; % concatenate the subject indexes.

    %% Processing the MP
    MP_filename_timeseries = strcat([basedir,'results/',atlas_name,'/tcs/sub-',char(control_finished_demogpraph_table(i,2).Variables),'_ses-',session_for_MP,'_atlas_schae232_roi_time_series_scrubbed.txt'])
    time_series = load (MP_filename_timeseries);
    time_series = time_series';

    time_series = time_series - mean(time_series, 1);

%     time_series_zscore = zscore(time_series);

    MP_time_series = [MP_time_series;time_series]; % concatenate the time series.
    temp_subjInd = i * ones(size(time_series,1),1); 
    temp_MPsubjInd = [temp_MPsubjInd;temp_subjInd]; % concatenate the subject indexes.

end

concTS = [MP_time_series;PL_time_series]; % 38480*48 single
MPsubjInd = [temp_MPsubjInd;zeros(size(temp_PLsubjInd))];
PLsubjInd = [zeros(size(temp_MPsubjInd));temp_PLsubjInd];
subjInd = [temp_MPsubjInd;temp_PLsubjInd];
TR = 520;

% Move the subcortical regions to the last to make it consistent with the 'data/sch232_to_yeo.csv'. (index: 201-232)
temp = concTS(:,1:32);
concTS(:,1:32) = [];
concTS = [concTS,temp];
save('data/main.mat','nparc','nscans','nsubjs','subj_scanInd','concTS','MPsubjInd','PLsubjInd','subjInd','TR','atlas_name','rep_time');




