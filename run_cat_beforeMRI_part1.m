function run_cat_beforeMRI_part1()



DoBDMdemo=1;
DoBDM=1;
DoSortBDM=1;
DoResponse=1;
DoTrainingDemo=1;
DoProbeDemo=1;
DoOrganizeProbe=1;
use_eyetracker=0;
DoSave=1;
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================= by Rotem Botvinik Nezer July 2016 ===============
% = = = = = = = = = =modified by Shiran Oren 01.2017 = = = = = = = = = = = = = = = = = = = = = =

% This function runs all the parts of the caau_snacks experiment
% DAY1: BDM, BDM sorting, training, BDM fractals, probe, recognition (with confidence levels), personal details, BDM resolve,
% probe resolvse.

% DAY30: probe, recognition (with confidence levels), BDM2, BDM resolve,
% probe resolve.

% The try-catch is because the mac caused an error with opening the screen
% from time to time, so we want to prevent it from failing.

% This version is for running only 40 items in training,
% 16 runs in the training session

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'BDM_Snacks'
% % %   'BDM_SnacksDemo'

% % %   'sort_BDM'


% % %   --- Other codes: ---
% % %   'CenterText'
% % %   'mysetdefaultbutton'
% % %   'myinputdlg'
% % %   'mygetnicedialoglocation'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FOLDERS REQUIRED TO RUN PROPERLY: ------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'Output': the folder for the output files- results.
% % %   'Stim': with the bmp files of all the stimuli for the cue-approach
% % %    task (old stimuli).
% % %   'Instructions': a folder with the png images of the instructions
% % %   for each part

tic

rng shuffle

% =========================================================================
% Get input args and check if input is ok
% =========================================================================

% %---dummy info for debugging purposes --------
% subjectNum =  '999'; % for debugging. Real subjects should start from 101
% and on
% subjectNum = '998'; % to test both order 1 and 2
% isMRI = 0; % 0- not MRI; 1- MRI experiment (and then test_comp = 1)
% sessionNum = 1; % if not follow-up
% mainPath = pwd;

% - - - - - - - - - - - - - - - - - - - - - - - -
% fixed variables for this behavioral experiment
% - - - - - - - - - - - - - - - - - - - - - - - -
experimentName = 'cat_between_MRI'; % hard-coded since true for all subjects
sessionNum = 1; % this is the code for day1 - the first session
isMRI = 0; % this is a behavioral study. change if switch to MRI...
mainPath = pwd;
outputPath = [mainPath '/Output'];

% get time and date
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

subjectID_ok = 0;

while subjectID_ok == 0
    subjectNum = input('Subject number (only the digits): ');
    while isempty(subjectNum)
        disp('ERROR: no value entered. Please try again.');
        subjectNum = input('Subject number (only the digits):');
    end
    % order number (same across all tasks\runs for a single subject. Should be
    % counterbalanced between 1,2 across subjects)
    % define order according to subjectNum being even/odd
    if mod(subjectNum, 2) == 0
        order = 2;
    else
        order = 1;
    end

    % Assign reward type
    % --------------------------
    % give order value of '1' or '2' for subjects with odd or even ID, respectively
    if subjectNum <200 % subject code is 100+
        CueType = 1;
    else % subject code is 200+
        CueType = 2;
    end
            
    subjectID = [experimentName '_' num2str(subjectNum)];
    disp(['subjectID is: ' num2str(subjectNum)]);
    disp(['order is: ' num2str(order)]);
    disp(['cue type is: ' num2str(CueType)]);
    disp(['session is: ' num2str(sessionNum)]);
    
    % read files from the output folder to make sure the subject number is
    % correct and there aren't any files for this subject (if session1) or
    % there are files from previous sessions and not for the current one (for
    % follow up sessions)
    subject_files = dir([outputPath '/' subjectID '*']);
    if ~isempty(subject_files)
        warning_msg = ['There are already ' num2str(length(subject_files)) ' files for this subject- ' subjectID '. Please make sure you entered the right number'];
        set(groot,'defaultUicontrolFontSize', 16);
        warning_answer = questdlg(warning_msg,'Warning!','it is OK', 'change subject number','it is OK');
    else
        subjectID_ok = 1;
    end
    
    if exist('warning_msg', 'var') && strcmp(warning_answer, 'it is OK')
        subjectID_ok = 1;
    end
    
end

% open a txt file for crashing logs
fid_crash = fopen([outputPath '/' subjectID '_crashingLogs' num2str(sessionNum) '_' timestamp '.txt'], 'a');

% =========================================================================
% BDM demo
% =========================================================================
if DoBDMdemo
    crashedDemoBDM = 0;
    keepTrying = 1;
    demo_again=1;
    while demo_again
        while keepTrying < 10
            try
                BDM_SnacksDemo(subjectID,sessionNum)
                keepTrying = 10;
            catch
                sca;
                crashedDemoBDM = crashedDemoBDM + 1;
                keepTrying = keepTrying + 1;
                disp('CODE HAD CRASHED - BDM DEMO!');
            end
        end
        fprintf(fid_crash,'BDM demo crashed:\t %d\n', crashedDemoBDM);

        set(groot,'defaultUicontrolFontSize', 16)
        demo_again = questdlg('Do you want more practice?','No','Yes','No','Yes');
        if strcmp(demo_again, 'No')
            demo_again = 0;
        else
            demo_again = 1;
            keepTrying=1;
        end
    end
end

% =========================================================================
% BDM full
% =========================================================================
if DoBDM
    crashedBDM = 0;
    keepTrying = 1;
    while keepTrying < 10
        try
            BDM_Snacks(subjectID,sessionNum)
            keepTrying = 10;
        catch
            sca;
            crashedBDM = crashedBDM + 1;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - BDM!');
        end
    end
    fprintf(fid_crash,'BDM crashed:\t %d\n', crashedBDM);
    
end
% =========================================================================
% Sort stimuli according to the BDM ranking (BDM = PART 1)
% =========================================================================
if DoSortBDM
    sort_BDM(subjectID,order,outputPath);
end

% =========================================================================
%% do demo's outside the scanner 
% =========================================================================

% are you ready for the demos?
    set(groot,'defaultUicontrolFontSize', 16)
    questdlg('Start demos?','No','Yes','No','Yes');
  
%=================================
%   response to snacks 
%=================================
if DoResponse
    crashedDemoResponse = 0;
    keepTrying = 1;
    demo_again=1;
    while demo_again
        while keepTrying < 10
            try
            responseToSnacks_demo(subjectID)
            keepTrying = 10;
            catch
                sca;
                crashedDemoResponse = crashedDemoResponse + 1;
                keepTrying = keepTrying + 1;
                disp('CODE HAD CRASHED - RESPONSE TO SNACKS DEMO!');
            end
         end
        fprintf(fid_crash,'training demo crashed:\t %d\n', crashedDemoResponse);

        set(groot,'defaultUicontrolFontSize', 16)
        demo_again = questdlg('Do you want more practice?','No','Yes','No','Yes');
        if strcmp(demo_again, 'No')
            demo_again = 0;
        else
            demo_again = 1;
            keepTrying=1;
        end
    end
end


% =========================================================================
% training demo
% =========================================================================

if DoTrainingDemo
    crashedDemoTraining = 0;
    keepTrying = 1;
    demo_again=1;
    while demo_again
        while keepTrying < 10
            try 
                trainingDemo(subjectID,order,mainPath,isMRI,use_eyetracker,CueType);
                keepTrying=10;
            catch
                sca;
                crashedDemoTraining = crashedDemoTraining + 1;
                keepTrying = keepTrying + 1;
                disp('CODE HAD CRASHED - TRAINING DEMO!');
            end
        end
        fprintf(fid_crash,'training demo crashed:\t %d\n', crashedDemoTraining);

        set(groot,'defaultUicontrolFontSize', 16)
        demo_again = questdlg('Do you want more practice?','No','Yes','No','Yes');
        if strcmp(demo_again, 'No')
            demo_again = 0;
        else
            demo_again = 1;
        end
    end
end



%=================================
%   probe_demo 
%=================================
if DoProbeDemo
    crashedDemoProbe=0;
    keepTrying = 1;
    demo_again=1;
    while demo_again
        while keepTrying < 10
            try
                probeDemo(subjectID, order, mainPath, isMRI, sessionNum, use_eyetracker);
                keepTrying = 10;
            catch
                sca;
                crashedDemoProbe = crashedDemoProbe + 1;
                keepTrying = keepTrying + 1;
                disp('CODE HAD CRASHED - PROBE DEMO!');
            end
        end
        fprintf(fid_crash,'probe demo crashed:\t %d\n', crashedDemoProbe);

        set(groot,'defaultUicontrolFontSize', 16)
        demo_again = questdlg('Do you want more practice?','No','Yes','No','Yes');
        if strcmp(demo_again, 'No')
            demo_again = 0;
        else
            demo_again = 1;
            keepTrying=1;
        end
    end
end

%=================================
%   organize probe 
%=================================

% Run blocks. Before each block, stimuli need to be organized in a list and
% divided to the required number of runs
numBlocks = 2; % Define how many blocks for the probe session. Each block includes all comparisons, one time each.
numRunsPerBlock = 1; % Define the required number of runs per block
if DoOrganizeProbe==1
    for ind = 1:numBlocks
        block = (sessionNum-1)*numBlocks + ind;
        % Organize the stimuli for the probe
        % ===================================
         organizeProbe(subjectID, order, mainPath, block, numRunsPerBlock);
    end % end for block = 1:numBlocks

end


if DoSave
    FileForDropbox=dir(['Output/*' subjectID '*']);
    DropboxPath='~/Dropbox/experimentsoutput/shiran/';
    for file =1:size(FileForDropbox,1)
        CurrentFile=['Output/' FileForDropbox(file).name];
        system(['cp ' CurrentFile ' ' DropboxPath])
    end
end

fclose(fid_crash);

end % end function
