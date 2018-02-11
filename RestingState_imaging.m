clear all

Screen('Preference', 'SkipSyncTests', 1);


rng shuffle
% =========================================================================
% Get input args and check if input is ok
% =========================================================================

okEyeTracker = [0 1];

subjectID = input('Subject code: ','s');
while isempty(subjectID)
    disp('ERROR: no value entered. Please try again.');
    subjectID = input('Subject code:','s');
end


subject_num = subjectID(end-2:end);
subject_num = str2double(subject_num);
sessionNum=1;


use_eyetracker = input('Do you want eyetracker (1-yes, 0-no)?: ');
while isempty(use_eyetracker) || sum(okEyeTracker==use_eyetracker)~=1
    disp('ERROR: input must be 0 or 1 . Please try again.');
    use_eyetracker = input('Do you want eyetracker (1-yes, 0-no)?: ');
end

fprintf(['\n\nsubjectID is: ' num2str(subject_num) '\n']);
disp(['sessionNum is: ' num2str(sessionNum)]);
disp(['use eyetracker is: ' num2str(use_eyetracker)]);
fprintf('\n')

GoOn=input('are all variables ok? (1-yes, 0-no)');
if GoOn==0
    error('please check you numbers and start again')
end
  
    
% =========================================================================
% set the computer and path
% =========================================================================

test_comp = 1;

% Set main path
mainPath = pwd;
outputPath = [mainPath '/Output'];
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
Instructions = dir([mainPath '/Instructions/Rest.JPG' ]);
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
    edfFile='Rest.edf';
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
    % escapeKey = KbName('space');
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
end; % end if test_comp == 1
DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected

tic


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
%% 'Run 8 min of black screen'
%---------------------------------------------------------------

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

tic
WaitSecs(480);
toc

Screen('TextSize',w, 40);
CenterText(w,'Thank you!', white,0,0);
Screen(w,'Flip');
WaitSecs(3);



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
    
    eyetracking_filename = [outputPath '/' subjectID '_RestingState_' num2str(sessionNum)  '_' timestamp '.edf'];
    movefile(edfFile,eyetracking_filename);
end

