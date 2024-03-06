clear; clc

% Basic demographic information
demograph_table = readtable('./data/SubjectList_17-0114.xlsx','Sheet','SignedConsent','VariableNamingRule','preserve','Format','auto');
demograph_table = demograph_table(1:111, :); % remove the redundant rows.
finished_demograph_table = demograph_table(find(demograph_table.Finished==1),:); % Select the subjects finished the pipelines;
control_finished_demogpraph_table = finished_demograph_table(strcmp(finished_demograph_table.Case_Control,'Control'),:); % Select only the Controls.

fMRI_session_PL_MP = zeros(size(control_finished_demogpraph_table,1), 1);

% To differenciate the Placebo session and MP session.
for i=1:size(control_finished_demogpraph_table,1)
    duration_time = datetime(control_finished_demogpraph_table.MP_Day(i))-datetime(control_finished_demogpraph_table.PO_day(i));
    if hours(duration_time) > 0
        fMRI_session_PL_MP(i)=1; % 1st session: PL, 2nd session: MP;
    else
        fMRI_session_PL_MP(i)=0; % 1st session: MP, 2nd session: PL;
    end
end


control_name_list = control_finished_demogpraph_table(:,2);
writetable(control_name_list,'./results/control_name.txt','WriteVariableNames',0);
read_me = '0_represent_first_session_MP_second_session_PL';
save('results/demographic.mat','control_finished_demogpraph_table','fMRI_session_PL_MP','read_me');