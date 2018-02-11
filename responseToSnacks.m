% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% ================== Created by Rotem Botvinik July 2015 ==================
% ============== Edited to include eyetracking on July 2016 ===============
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the trained stimuli for scanning fMRI to compare
% activations and representations of the stimuli before and after the
% training.
% fixed ISI = 5
% each image is represented for 2 secs, and then 6 secs of fixation cross
% The "dummy" task is to count either "one" items or "several" items (Does
% the item contain one or several items inside a new closed bag?)

% The code saved onset times to txt file with onsets starting from 0
% instead of 2 (2 secs of fixation before the first image).
% it was fixed after:
% Experiment group: subject MRI_snacks_112 in session 1, subject
% MRI_snacks_111 in session 2.
% Control group: subject MRI_snacks2_110 in session 1, subject
% MRI_snacks2_106 in session 2.


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order%d.txt', order --> Created by the 'sortBdm_Israel' function
%    Misc\oneSeveral.mat  (in which there is a vector named 'oneSeveral' in
%    which there are 1s and 2s whether each item is 'one' (1) or 'several'
%    (2)(ordered by abc of item names)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [outputPath '/' subjectID '_responseToSnacks_session' num2str(sessionNum) '_' representationSession '_' timestamp '.txt']


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'testMRI';
% order = 1;
% test_comp = 1;
% sessionNum = 1;
% mainPath = 'D:\Rotem\Dropbox\Rotem\experiments\BMI_MRI_snacks_40\BMI_MRI_snacks_40';


clear all

Screen('Preference', 'SkipSyncTests', 1);


rng shuffle
% =========================================================================
% Get input args and check if input is ok
% =========================================================================

oksessionNum = [1 2 3 4 5];
%okComputer = [0 1 2 3 4 5];
%okOrder = [1 2];

okOneOrSeveral = [1 2];

subjectID = input('Subject code: ','s');
while isempty(subjectID)
    disp('ERROR: no value entered. Please try again.');
    subjectID = input('Subject code:','s');
end
% order number (same across all tasks\runs for a single subject. Should be
% counterbalanced between 1,2 across subjects)

subject_num = subjectID(end-2:end);
subject_num = str2double(subject_num);
if mod(subject_num,2)==0
    order = 2;
else
    order = 1;
end

sessionNum=1;

beforeAfter = input('Is it the representation before (1) or after (2) training?: ');
while isempty(beforeAfter)
    disp('ERROR: input must be 1 or 2 . Please try again.');
    beforeAfter = input('Is it the representation before (1) or after (2) training?: ');
end

RunNum = input('What is the run number (1-2)?: ');
while isempty(RunNum)
    disp('ERROR: input must be 1-4 . Please try again.');
    RunNum = input('What is the run number (1-2)?: ');
end

use_eyetracker = input('Do you want eyetracker (1-yes, 0-no)?: ');
while isempty(use_eyetracker)
    disp('ERROR: input must be 0 or 1 . Please try again.');
    use_eyetracker = input('Do you want eyetracker (1-yes, 0-no)?: ');
end


fprintf(['\n\nsubjectID is: ' num2str(subject_num) '\n']);
disp(['sessionNum is: ' num2str(sessionNum)]);
disp(['beforeAfter is: ' num2str(beforeAfter)]);
disp(['RunNum is: ' num2str(RunNum)]);
disp(['use eyetracker is: ' num2str(use_eyetracker)]);
fprintf('\n')

GoOn=input('are all variables ok? (1-yes, 0-no) ');
if GoOn==0
    error('please check you numbers and start again')
end

% decide whether it's a "one" counting run or a "several" counting run
switch order
    case 1
        if mod(randi(2),2)==1
            oneOrSeveral = 1;
        else
            oneOrSeveral = 2;
        end
    case 2
        if mod(randi(2),2)==1
            oneOrSeveral = 2;
        else
            oneOrSeveral = 1;
        end
end % end switch order

switch beforeAfter
    case 1
        representationSession = 'before';
    case 2
        representationSession = 'after';
end

% =========================================================================
% set the computer and path
% =========================================================================

test_comp = 1;
RepetitionNum=1;
% Set main path
mainPath = pwd;

outputPath = [mainPath '/Output'];

% - - - - - - - - - - - - - - - -
% check that subjectID is correct
% - - - - - - - - - - - - - - - -

subject_stopgolist = dir([outputPath '/*' subjectID '*StopGoList*']);
switch oneOrSeveral
    case 1
        subject_current_run_file = [outputPath '/' subjectID '_responseToSnacks_session' num2str(sessionNum) '_' representationSession '_one_*.txt'];
    case 2
        subject_current_run_file = [outputPath '/' subjectID '_responseToSnacks_session' num2str(sessionNum) '_' representationSession '_several_*.txt'];
end
subject_current_run_file = dir(subject_current_run_file);

if isempty(subject_stopgolist)
    warning_msg = ['There is no stopgo list for this subject- ' subjectID];
    set(groot,'defaultUicontrolFontSize', 16);
    warning_answer = questdlg(warning_msg,'Warning!','continue', 'stop','continue');
end

if exist('warning_msg', 'var') && strcmp(warning_answer, 'stop')
    error('Stopped following the user request');
end


HideCursor;

%==========================================================
%% 'INITIALIZE Screen variables to be used in each task'
%==========================================================

Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize = 32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;

sizeFactor = 0.8;
stimW = 576*sizeFactor;
stimH = 432*sizeFactor;
rect = [xcenter-stimW/2 ycenter-stimH/2 xcenter+stimW/2 ycenter+stimH/2];


% % Set the colors
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
Green = [0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% Set up screen positions for stimuli
[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;


% Text settings
theFont ='Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
if oneOrSeveral == 1
    Instructions = dir([mainPath '/Instructions/functional_response_to_snacks_one.JPG' ]);
else
    Instructions = dir([mainPath '/Instructions/functional_response_to_snacks_several.JPG' ]);
end
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);

Trigger = dir([mainPath '/Instructions/Trigger_image.png' ]);
Trigger_image = imread([mainPath '/Instructions/' Trigger(1).name]);

%---------------------------------------------------------------
%% 'GLOBAL VARIABLES'
%---------------------------------------------------------------
c = clock;
hr = sprintf('%02d', c(4));
min = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',min,'m'];

%-----------------------------------------------------------------
%% Initializing eye tracking system %
%-----------------------------------------------------------------
% use_eyetracker=1; % set to 1/0 to turn on/off eyetracker functions
if use_eyetracker
    dummymode=0;
    
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(w);
    % Disable key output to Matlab window:
    
    el.backgroundcolour = black;
    el.backgroundcolour = black;
    el.foregroundcolour = white;
    el.msgfontcolour    = white;
    el.imgtitlecolour   = white;
    el.calibrationtargetcolour = el.foregroundcolour;
    EyelinkUpdateDefaults(el);
    
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end;
    
    [v,ELversion]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', ELversion );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    edfFile='Res2stim.edf';
    Eyelink('Openfile', edfFile);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    
    %     % STEP 5
    %     % start recording eye position
    %     Eyelink('StartRecording');
    %     % record a few samples before we actually start displaying
    %     WaitSecs(0.1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish Initialization %
    %%%%%%%%%%%%%%%%%%%%%%%%%
end




%---------------------------------------------------------------
%%   'PRE-TRIAL DATA ORGANIZATION'
%---------------------------------------------------------------

%   'Reading in the sorted BDM list - defines which items should be shown'
% - - - - - - - - - - - - - - -
file = dir([outputPath '/*' subjectID '_stopGoList_trainingstim.txt']);
fid = fopen([outputPath '/' sprintf(file(length(file)).name)]);
vars = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
fclose(fid);


%   'Reading in the sorted BDM list - for all the stimNames to rightly shuffle the oneSeveral matrix'
% - - - - - - - - - - - - - - -
file = dir([outputPath '/*' subjectID sprintf('_stopGoList_allstim_order%d.txt', order)]);
fid2 = fopen([outputPath '/' sprintf(file(length(file)).name)]);
vars2 = textscan(fid2, '%s %d %d %f %d') ;% these contain everything from the sortbdm
fclose(fid2);
allStimNames = vars2{1};
[ABCallStimNames, ABCind] = sort(allStimNames);


LastFixation=7;

%---------------------------------------------------------------
%% 'Write output file header'
%---------------------------------------------------------------
if oneOrSeveral ==1
    fid1 = fopen([outputPath '/' subjectID '_responseToSnacks_session' num2str(sessionNum) '_run' num2str(RunNum) '_' representationSession '_one_' timestamp '.txt'], 'a');
    fprintf(fid1,'subjectID\torder\trun\tbefore_or_after\tsession(visit)\titemName\tbidInd\toneOrSeveral?\tonsettime\tfixationTime\n'); %write the header line
else
    fid1 = fopen([outputPath '/' subjectID '_responseToSnacks_session' num2str(sessionNum) '_run' num2str(RunNum) '_' representationSession '_several_' timestamp '.txt'], 'a');
    fprintf(fid1,'subjectID\torder\trun\tbefore_or_after\tsession(visit)\titemName\tbidInd\tisOne?\tonsettime\tfixationTime\n'); %write the header line
end


%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

Screen('TextSize',w, 40);

Screen('PutImage',w,Instructions_image);
Screen(w,'Flip');

noresp = 1;
while noresp,
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp,
        noresp = 0;
    end;
end

if test_comp == 1
    Screen('PutImage',w,Trigger_image);
    Screen('Flip',w);
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
end; % end if test_comp == 1
DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected


anchor = GetSecs;
repeat=1;
tic

%   baseline fixation cross
% - - - - - - - - - - - - -
prebaseline = GetSecs;
% baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
baseline_fixation_dur = 2; % Need to modify based on if first few volumes are saved or not
while GetSecs < prebaseline+baseline_fixation_dur
    %    Screen(w,'Flip', anchor);
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
    
end

postbaseline = GetSecs;
baseline_fixation = postbaseline - prebaseline;

% start recording eye position
%-----------------------------
if use_eyetracker
    
    Eyelink('StartRecording');
    WaitSecs(.05);
    
    %   Eyelink MSG
    % ---------------------------
    % messages to save on each trial ( trial number, onset and RT)
    Eyelink('Message',['SYNCTIME at run start:',num2str(GetSecs)]); % mark start time in file
end

%---------------------------------------------------------------
%% loop through num of repetition of all images
%---------------------------------------------------------------

% SHUFFLE the stimuli order for random presentation. because half of
% the items are presented in each run, do the shuffle in the odd run
% and save the order. is it is an even run read the list that was
% saved last run and present its second half.
stimNamesAll = vars{1};
[shuff_namesAll,shuff_indAll] = Shuffle(stimNamesAll);
bidIndex = vars{3};
shuff_bidIndexAll = bidIndex(shuff_indAll);

% read the oneSeveral matrix
load([mainPath '/Misc/oneSeveral.mat']);
oneSeveralAll(ABCind) = oneSeveral;
oneSeveralAll = oneSeveralAll(shuff_bidIndexAll);

shuff_names=shuff_namesAll;
shuff_bidIndex=shuff_bidIndexAll;
oneSeveral=oneSeveralAll;


numStimuli = length(shuff_names);


%     %   'load onsets'
%     %---------------------------
%     r = Shuffle(1:4);
%     load(['Onset_files/resp_onset_length_' num2str(numStimuli) '_' num2str(r(1)) '.mat']);
% for the fixed timing version, create the vector of onsets with 6sec
% gaps - 2 sec stimuli & 7 sec ITI
onsetlist=0:9:9*numStimuli;

actual_onset_time = zeros(numStimuli,1);
fixationTime = zeros(1,numStimuli);
stimLength = 2; % the length of each stimulus presentation
lastOnset = onsetlist(numStimuli-1) ;



%---------------------------------------------------------------
%% 'LOAD image arrays'
%---------------------------------------------------------------
imgArrays = cell(1, numStimuli);
for ind = 1:numStimuli
    imgArrays{ind} = imread([mainPath '/Stim/' shuff_names{ind}],'bmp');
end


%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------
runStart = GetSecs;

% for stimInd = 1:3 % for debugging
for stimInd = 1:numStimuli
    %-----------------------------------------------------------------
    % display image
    Screen('PutImage',w,imgArrays{stimInd},rect);
    image_start_time = Screen('Flip',w,runStart + onsetlist(stimInd)); % display an image each
    actual_onset_time(stimInd,1) = image_start_time - anchor + baseline_fixation;
    
    if use_eyetracker && runStart + onsetlist(stimInd) < GetSecs
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',['trial ',num2str(stimInd),' stim: ',shuff_names{stimInd},' time:',num2str(GetSecs), ' onset: ', num2str(actual_onset_time(stimInd,1))]);
    end
    
    %-----------------------------------------------------------------
    % show fixation ITI
    Screen('TextSize',w, 60);
    Screen('DrawText', w, '+', xcenter, ycenter, white);
    Screen('Flip',w,runStart+onsetlist(stimInd)+stimLength);
    fixationTime(stimInd) = GetSecs - anchor + baseline_fixation;
    
    if use_eyetracker && runStart + onsetlist(stimInd) < GetSecs
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',['ITI fixation cross time:',num2str(GetSecs), ' onset: ', num2str(fixationTime(stimInd))]);
    end
    
    %---------------------------------------------------------------
    %% 'Write to output file'
    %---------------------------------------------------------------
    fprintf(fid1,'%s\t %d\t %d\t %s\t %d\t %s\t %d\t %d\t %f\t %f\t\n', subjectID, order, repeat, representationSession, sessionNum, shuff_names{stimInd}, shuff_bidIndex(stimInd), oneSeveral(stimInd), actual_onset_time(stimInd,1),fixationTime(stimInd));
    
end % end for stimInd = 1:numStimuli


fclose(fid1);
WaitSecs(LastFixation); % for the last fixation to last like the ones before

display('run rnd');

toc


ShowCursor;


if use_eyetracker
    
    %---------------------------------------------------------------
    %%   Finishing eye tracking system %
    %---------------------------------------------------------------
    
    % STEP 7
    %---------------------------
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    WaitSecs(.1);
    Eyelink('CloseFile');
    
    
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', edfFile );
        rdf;
    end
    
    
    if dummymode==0
        oneOrSeveral_strings = {'one', 'several'};
        movefile(edfFile,[outputPath,'/', sprintf('%s_resopnseToSnacks_session_%d_%s_run_%d_%s.edf',subjectID,sessionNum,representationSession,RunNum,timestamp)]);
    end;
end


%   outgoing msg & closing
% ------------------------------
if RunNum == 2 % if this is not the last run
    CenterText(w,'Great Job. Thank you!',Green, 0,-270);
    Screen('Flip',w);
end
% clc;
WaitSecs(3);
Screen('CloseAll');
ShowCursor;

% Save variables to mat file
outfile = strcat(outputPath,'/', sprintf('%s_resopnseToSnacks_session_%d_%s_run_%d_%s.mat',subjectID,sessionNum,representationSession,RunNum,timestamp));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
run_info.script_name = mfilename;

% calculate oneSeveral information

if oneOrSeveral == 1
    howManyItems = sum(oneSeveral==1);
    run_info.ShumTask='one';
    fprintf('There were %d "one" items in this run\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n',howManyItems); % many 'enters' so the subject will not see the answer
else
    howManyItems = size(oneSeveral,2)-sum(oneSeveral==1);
    fprintf('There were %d "several" items in this run\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n',howManyItems); % many 'enters' so the subject will not see the answer
    run_info.ShumTask='several';
end
run_info.RealAnswer=howManyItems;


if RunNum<2
    fprintf(['\nyour next run is: ' num2str(RunNum+1),'\n\n']);
else
    disp('continue to next part\n\n');
end


KbQueueCreate;
KbQueueFlush;

SubAnswer= input('How many items did the participant count?: ');
while isempty(SubAnswer)
    disp('ERROR: input must be an integer . Please try again.');
    SubAnswer = input('How many items did the participant count?: ');
end

run_info.SubAnswer=SubAnswer;


save('outfile','run_info');





% end % end function