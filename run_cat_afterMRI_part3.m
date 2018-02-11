function run_cat_afterMRI_part3()

GetFromDropbox=1;
DoRecognitiondemo=1;
DoRecognition=1;
DoPersonal=1;
DoAEBQ=1;
DoResolveBDM=1;
DoResolveProbe=1;
DoSave=1;
use_eyetracker=0;




% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================= by Rotem Botvinik Nezer July 2016 ===============
% = = = = = = = = = =modified by Shiran Oren 01.2017 = = = = = = = = = = = = = = = = = = = = = =

% This function runs all the parts of the cat_snacks experiment
% DAY1: recognition, personal details and winnings
 

% The try-catch is because the mac caused an error with opening the screen
% from time to time, so we want to prevent it from failing.

% This version is for running only 40 items in training,
% 16 runs in the training session

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'recognition_confidence_demo'
% % %   'recognition_confidence'

% % %   'disp_resolveBDM'
% % %   'disp_probeResolve'
% % %   'personal_details'


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
mainPath = pwd;
outputPath = [mainPath '/Output'];
isMRI=0;

% get time and date
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];


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
      



% open a txt file for crashing logs
fid_crash = fopen([outputPath '/' subjectID '_crashingLogs' num2str(sessionNum) '_' timestamp '.txt'], 'a');



%% =========================================================================
% get probe files from dropbox
% =========================================================================
if GetFromDropbox
    FileFromDropbox=dir(['Output/*' subjectID '*_probe_block_*']);
    DropboxPath='~/Dropbox/experimentsoutput/shiran/';
    for file =1:size(FileFromDropbox,1)
        CurrentFile=[DropboxPath , FileFromDropbox(file).name];
        system(['cp ' CurrentFile ' ' outputPath '/'])
    end
end




%% =========================================================================
% recognition demo
% =========================================================================
if DoRecognitiondemo
    crashedDemoRecognition = 0;
    keepTrying = 1;
    demo_again=1;
    while demo_again
        while keepTrying < 10
            try 
                recognition_confidence_demo(isMRI,mainPath,use_eyetracker)
                keepTrying = 10;
            catch
                sca;
                crashedDemoRecognition = crashedDemoRecognition + 1;
                keepTrying = keepTrying + 1;
                disp('CODE HAD CRASHED - Recognition DEMO!');
            end
        end
        fprintf(fid_crash,'recognition demo crashed:\t %d\n', crashedDemoRecognition);

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
% recognition full 
% =========================================================================
if DoRecognition
    crashedRecognition = 0;
    keepTrying = 1;
    while keepTrying < 10
        try 
            recognition_confidence(subjectID,isMRI,mainPath,order,sessionNum,use_eyetracker)
            keepTrying = 10;
        catch
            sca;
            crashedRecognition = crashedRecognition + 1;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - Recognition!');
        end
    end
    fprintf(fid_crash,'Recognition crashed:\t %d\n', crashedRecognition);
    
end








% =========================================================================
% personal details
% =========================================================================
if DoPersonal
    personal_details(subjectID, order, outputPath, sessionNum)
end



% =========================================================================
% Eating Questionaire
% =========================================================================
if DoAEBQ
    web('https://docs.google.com/forms/d/e/1FAIpQLSfQvgQN9S6DFUrh28dPFI7IFR1djcaPZb4xN57yrx6HTCmQVQ/viewform', '-browser')
end


% =========================================================================
% resolve probe
% =========================================================================

if DoResolveProbe
    crashedResolveProbe = 0;
    keepTrying = 1;
    while keepTrying < 10
        try 
             disp_probeResolve(subjectID, sessionNum, outputPath, 2)
            keepTrying=10;
        catch
            sca;
            crashedResolveProbe = crashedResolveProbe + 1;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - disp_probeResolve!');
        end
    end
    fprintf(fid_crash,'disp_probeResolve crashed:\t %d\n', crashedResolveProbe);
    
end




% =========================================================================
% resolve BDM
% =========================================================================

if DoResolveBDM
    crashedResolveBDM = 0;
    keepTrying = 1;
    while keepTrying < 10
        try 
            disp_resolveBDM(subjectID, mainPath, sessionNum)
            keepTrying=10;
        catch
            sca;
            crashedResolveBDM = crashedResolveBDM + 1;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - disp_resolveBDM!');
        end
    end
    fprintf(fid_crash,'disp_resolveBDM crashed:\t %d\n', crashedResolveBDM);
    
end





% =========================================================================
% save files to dropbox
% =========================================================================
if DoSave
    FileFromDropbox=dir(['Output/*' subjectID '*']);
    DropboxPath='~/Dropbox/experimentsOutput/shiran/';
    for file =1:size(FileFromDropbox,1)
        CurrentFile=['Output/' FileFromDropbox(file).name];
        system(['cp ' CurrentFile ' ' DropboxPath])
    end

    fclose(fid_crash);
end


end % end function
