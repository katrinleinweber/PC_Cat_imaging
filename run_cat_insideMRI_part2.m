

%=================================
%%   get the files from dropbox
%=================================
experimentName = 'cat_between_MRI'; % hard-coded since true for all subjects
subjectNum = input('Subject number (only the digits): ');
subjectID = [experimentName '_' num2str(subjectNum)];

DropboxPath='~/Dropbox/experimentsoutput/shiran/';
FileFromDropbox=dir([DropboxPath  subjectID '*']);
for file =1:size(FileFromDropbox,1)
    CurrentFile=[FileFromDropbox(file).name];
    system(['cp ' DropboxPath CurrentFile ' Output/' ])
end




%=================================
%%   response to snacks - 4 repetitions
%=================================
responseToSnacks

% =========================================================================
%% training - 6 repetitions
% =========================================================================
training_imaging


%=================================
%   anatomical scans
%=================================
%MPRAGE & FLAIR



%=================================
%%   response to snacks - 2 repetitions
%=================================
responseToSnacks


%=================================
%%   probe - 4 repetitions
%=================================
probe_imaging



%=================================
%%   resting state - 1 repetition
%=================================
RestingState_imaging


%%
%=================================
%   copy the files to dropbox
%=================================

FileForDropbox=dir(['Output/*' subjectID '*']);
DropboxPath='~/Dropbox/experimentsoutput/shiran/';
for file =1:size(FileForDropbox,1)
    CurrentFile=['Output/' FileForDropbox(file).name];
    system(['cp ' CurrentFile ' ' DropboxPath])
end
%%
