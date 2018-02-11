function responseToSnacks_demo(subjectID)

Screen('Preference', 'SkipSyncTests', 1);

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% ================== Created by Rotem Botvinik July 2015 ==================
% ============== Edited to include eyetracking on July 2016 ===============
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


% This function runs the trained stimuli for scanning fMRI to compare
% activations and representations of the stimuli before and after the
% training.
% fixed ISI = 7
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



tic

rng shuffle
% =========================================================================
% Get input args and check if input is ok
% =========================================================================
use_eyetracker = 0; % set to 1/0 to turn on/off eyetracker functions

% decide whether it's a "one" counting run or a "several" counting run
if mod(randi(2),2)==1 % this is an odd runNum
oneOrSeveral = 1;   
else % this is a pair runNum
oneOrSeveral = 2;   
end

% =========================================================================
% set the computer and path
% =========================================================================

% Set main path
mainPath = pwd;

outputPath = [mainPath '/Output'];


% - - - - - - - - - - - - - - - -
% check that subjectID is correct
% - - - - - - - - - - - - - - - -

subject_stopgolist = dir([outputPath '/' subjectID '*stopGoList*']);
if isempty(subject_stopgolist)
    warning_msg = ['There is no stopgo list for this subject- ' subjectID];
    set(groot,'defaultUicontrolFontSize', 16);
    warning_answer = questdlg(warning_msg,'Warning!','stop', 'continue','stop'); 
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
% green = [0 255 0];

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

Instructions = dir([mainPath '/Instructions/oneSeveralDemonstration.jpg' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);

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

% SHUFFLE the stimuli that the subject thought were old, for the goNoGo recognition task
SitmList=dir('stim/demo/*.bmp');
stimNames = {SitmList.name};
[shuff_names] = Shuffle(stimNames);
numStimuli = length(stimNames);

%---------------------------------------------------------------
%% 'LOAD image arrays'
%---------------------------------------------------------------
imgArrays = cell(1, numStimuli);
for ind = 1:numStimuli
    imgArrays{ind} = imread([mainPath '/Stim/demo/' shuff_names{ind}],'bmp');
end



%---------------------------------------------------------------
%% create matrices before loop
%---------------------------------------------------------------
actual_onset_time = zeros(numStimuli,1);
fixationTime = zeros(1,numStimuli);
ISI = 7;
stimLength = 2; % the length of each stimulus presentation
onsetInterval = ISI + stimLength;
lastOnset = (numStimuli-1) * onsetInterval;
onsetlist = 0:onsetInterval:lastOnset;


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

DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected

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
%% 'Run Trials'
%---------------------------------------------------------------

runStart = GetSecs;

%for stimInd = 1:3 % for debugging
for stimInd = 1:numStimuli
    
    %-----------------------------------------------------------------
    % display image
    Screen('PutImage',w,imgArrays{stimInd},rect);
    % Screen(w,'Flip');
    image_start_time = Screen('Flip',w,runStart + onsetlist(stimInd)); % display an image each
    actual_onset_time(stimInd,1) = image_start_time - runStart + baseline_fixation;
    % WaitSecs(2); % display each image for 2 seconds
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
    fixationTime(stimInd) = GetSecs - runStart + baseline_fixation;
    % WaitSecs(5);
    if use_eyetracker && runStart + onsetlist(stimInd) < GetSecs
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',['ITI fixation cross time:',num2str(GetSecs), ' onset: ', num2str(fixationTime(stimInd))]);
    end
        
end % end for stimInd = 1:numStimuli

WaitSecs(ISI); % for the last fixation to last like the ones before

% ending screen
Screen('TextSize',w, 40);
CenterText(w,'How many did you count?', white,0,0);
Screen(w,'Flip');
WaitSecs(3);

CenterText(w,'Thank you!', white,0,0);
Screen(w,'Flip');
WaitSecs(3);

toc

ShowCursor;
Screen('CloseAll');


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
        eyetracking_filename = [outputPath '/' subjectID '_responseToSnacks_session_EyeTracking_' num2str(sessionNum) '_' representationSession '_' oneOrSeveral_strings{oneOrSeveral} '_' timestamp '.edf'];
        movefile(edfFile,[outputPath,'/', sprintf('%s_resopnseToSnacks_EyeTracking_session_%d_%s_run_%d_%s.edf',subjectID,sessionNum,representationSession,runNum,timestamp)]);
    end;
end



 end % end function