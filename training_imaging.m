% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ====================== by Rotem Botvinik May 2015 =======================
% ------------- edited on March 2017 to include eye-tracking --------------
% ------------- edited on Jan 2018 by Shiran Oren --------------
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the boost (cue-approach) training session,
% in which the items are shown on the screen while some of them (GO items) are
% paired with a beep. The subject should press a predefined button as fast
% as possible after hearing the beep.
% This session is composed of one run with two repetitions for each item.
% This version is for training  40 items.
% this vrsion include a condiion of pavlovian cue - a money reward instead
% of a neutral cue.


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_trainingstim.txt' ---> The file for training 40 items,
%   created by the function sort_BDM_Israel
%   '/Onset_files/train_onset_' num2str(r(1)) '.mat''  where r=1-4
%   all the contents of 'stim/' food images
%   'Misc/soundfile.mat'
%   'CenterText.m'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'training_run_' sprintf('%02d',RunNum) '_' timestamp '.txt'
%


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = '909';
% subjectID = '908'; % to test both order 1 and 2
% mainPath = pwd;
% isMRI = 0;
% total_num_runs_training = 8;
% Ladder1IN = 750;
% Ladder2IN = 750;

clear all
rng shuffle

experimentName='cat_between_MRI';
mainPath=pwd;
outputPath = [mainPath '/Output'];
isMRI=1;



% input checkers
subjectID_ok = 0;
while subjectID_ok == 0
    subjectNum = input('Subject number (only the digits): ');
    while isempty(subjectNum)
        disp('ERROR: no value entered. Please try again.');
        subjectNum = input('Subject number (only the digits):');
    end
    subjectID = [experimentName '_' num2str(subjectNum)];
    
    % Assign order
    % --------------------------
    % give order value of '1' or '2' for subjects with odd or even ID, respectively
    if mod(subjectNum,2) == 1 % subject code is odd
        order = 1;
    else % subject code is even
        order = 2;
    end
    
    % Assign reward type
    % --------------------------
    % give order value of '1' or '2' for subjects with odd or even ID, respectively
    if subjectNum <200 % subject code is 100+
        CueType = 1;
    else % subject code is 200+
        CueType = 2;
    end
    
    
    RunNum = input('Which run are you at: ');
    while isempty(RunNum)
        disp('ERROR: invalid input. Please try again.');
        RunNum = input('Which run are you at: ');
    end
    
    okEyetracker = [1 0];
    ask_if_want_eyetracker = input('Do you want eyetracking (1 - yes, 0 - no): ');
    while isempty(ask_if_want_eyetracker) || sum(okEyetracker == ask_if_want_eyetracker) ~=1
        disp('ERROR: input must be 1 or 0. Please try again.');
        ask_if_want_eyetracker = input('Do you want eyetracking (1 - yes, 0 - no): ');
    end
    use_eyetracker=ask_if_want_eyetracker; % set to 1/0 to turn on/off eyetracker functions
    
    % read files from the output folder to make sure the subject number is
    % correct and there aren't any files for this subject (if session1) or
    % there are files from previous sessions and not for the current one (for
    % follow up sessions)
    subject_files = dir([outputPath '/' subjectID '*training_run_*']);
    if size(subject_files,1)>RunNum-1
        warning_msg = ['There are already ' num2str(length(subject_files)) ' files for this subject- ' subjectID '. Please make sure you entered the right number'];
        set(groot,'defaultUicontrolFontSize', 16);
        warning_answer = questdlg(warning_msg,'Warning!','it is OK', 'change subject number','it is OK');
    else
        subjectID_ok = 1;
    end
    
    if exist('warning_msg', 'var') && strcmp(warning_answer, 'it is OK')
        subjectID_ok = 1;
    end
    
    
    fprintf(['\n\nsubjectID is: ' num2str(subjectNum) '\n']);
    disp(['cue type is: ' num2str(CueType)]);
    disp(['use eyetracker is: ' num2str(use_eyetracker)]);
    fprintf('\n')
    
    GoOn=input('are all variables ok? (1-yes, 0-no) ');
    if GoOn==0
        error('please check you numbers and start again')
    end
end
%---------------------------------------------------------------
%%   'GLOBAL VARIABLES'
%---------------------------------------------------------------

% Set the total number of runs in the experiment
total_number_of_runs=6;
num_runs_per_block=1;

% about timing
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

% about ladders
Step = 50;

% about timing
image_duration = 1; 
baseline_fixation = 2;
afterrunfixation = 4;


% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
Instructions = dir([mainPath '/Instructions/Training_Reward_Between_cue' num2str(CueType) '.JPG' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);

Trigger = dir([mainPath '/Instructions/Trigger_image.png' ]);
Trigger_image = imread([mainPath '/Instructions/' Trigger(1).name]);

Continue = dir([mainPath '/Instructions/Continue.jpg' ]);
Continue_image = imread([mainPath '/Instructions/' Continue(1).name]);

% -----------------------------------------------
%% 'INITIALIZE SCREEN'
%---------------------------------------------------------------
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize=32;
%[w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

%   colors
% - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
Green = [0 255 0];

Screen('FillRect', w, black);
Screen('Flip', w);


%---------------------------------------------------------------
%%  visual cue definitions
%---------------------------------------------------------------

pixels = get(0,'ScreenSize'); % get screen size
pixels = pixels(3:4);
rect = [0 0 pixels(1) pixels(2)];

% define rewarding visual cue
image=dir('./Stim/*.bmp');
image_read=imread(['./Stim/',image(1).name]);
image_size=size(image_read);
% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', w);
% Set up alpha-blending
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Load cue image
CueImageLocation = [mainPath '/gabor_cue.png'];
[CueImage, ~, CueImageAlpha] = imread(CueImageLocation);

% Place the transparency layer of the foreground image into the 4th (alpha)
% image plane. This is the bit which is needed, else the tranparent
% background will be non-transparent
transparencyFactor=0.7;  %change this factor to manipulate cue's transparency
CueImage(:, :, 4) = CueImageAlpha*transparencyFactor    ;

% Get the size of the Cue
[s1, s2, ~] = size(CueImage);
scaleFactor = 2; %change this factor to manipulate cue's size

% Make the images into a texture
imageTextureCue = Screen('MakeTexture', w, CueImage);
dstRect = CenterRectOnPointd([0 0 image_size(2)/4 image_size(2)/4] .* scaleFactor,...
    screenXpixels / 2-3, screenYpixels / 2+25);


%   text
% - - - - - -
theFont = 'Arial';
Screen('TextSize',w,40);
Screen('TextFont',w,theFont);
Screen('TextColor',w,white);


HideCursor;

%-----------------------------------------------------------------
%% Initializing eye tracking system %
%-----------------------------------------------------------------
% use_eyetracker=0; % set to 1/0 to turn on/off eyetracker functions
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
    
    [~,vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    task = GenFlags.Training.str; % change to current task - with GenFlags
    edfFile='Training.edf';
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


%%---------------------------------------------------------------
%%  'FEEDBACK VARIABLES'
%%---------------------------------------------------------------
KbName('UnifyKeyNames');
blue = 'b';
yellow = 'y';



%---------------------------------------------------------------
%%   'PRE-TRIAL DATA ORGANIZATION'
%---------------------------------------------------------------

%   'Reading in the sorted BDM list - defines which items will be GO/NOGO'
% - - - - - - - - - - - - - - -
file = dir([outputPath '/' subjectID '_stopGoList_trainingstim.txt']);
fid = fopen([outputPath '/' sprintf(file(length(file)).name)]);
vars = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
fclose(fid);


% set cue values for all runs- in order to control the mean of reward for each 'go' item
VarNum=4; %number of different rewards
LowReward=21;
RewardVec=repmat(LowReward:LowReward+VarNum-1,size(vars{1,1},1),floor(num_runs_per_block/VarNum)+mod(num_runs_per_block,VarNum)*10);
AllRewards=Shuffle(RewardVec.').';

if CueType==1
    CueStr='0';
    AllRewards=zeros(size(AllRewards));
elseif CueType==2
    CueStr='+';
else
    CueStr='-';
end



%---------------------------------------------------------------
%% start run
%---------------------------------------------------------------

%   'Write output file header'
%---------------------------------------------------------------
fid1 = fopen([outputPath '/' subjectID '_training_run_' sprintf('%02d',RunNum) '_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\torder\tRunNum\titemName\tonsetTime\ttrialType\tRT\trespInTime\tCueTime\tresponse\tfixationTime\tladder1\tladder2\tbidIndex\tBDMtrialIndex\tbidValue\n'); %write the header line

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

KbQueueCreate;


    if RunNum ==1   % if this is the first run, present full instructions
        Screen('PutImage',w,Instructions_image);
    else
        Screen('PutImage',w,Continue_image);
    end
    Screen(w,'Flip');
    
    noresp = 1;
    while noresp,
        [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
        if keyIsDown && noresp,
            noresp = 0;
        end;
    end
    
    if isMRI == 1
        Screen('PutImage',w,Trigger_image);
        Screen('Flip',w);
        escapeKey = KbName('t');
        while 1
            [keyIsDown,~,keyCode] = KbCheck(-1);
            max(keyCode);
            if keyIsDown && keyCode(escapeKey)
                break;
            end
        end
    end; % end if isMRI == 1
    
    DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected
    tic
    anchor = GetSecs ; % (before baseline fixation) ;


KbQueueCreate;
KbQueueFlush;


if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    %   Eyelink MSG
    % ---------------------------
end

%---------------------------------------------------------------
%%  'TRIAL PRESENTATION'
%---------------------------------------------------------------
%
%   trial_type definitions:
% - - - - - - - - - - - -
% 11 = High-Value GO
% 12 = High-Value NOGO
% 22 = Low-Value GO
% 24 = Low-Value NOGO
% 0 = not trained


%   'load onsets'
%---------------------------
r = Shuffle(1:4);
load(['Onset_files/train_onset_length_80_' num2str(r(1)) '.mat']);



%%   'pre-trial fixation'
%% ---------------------------

        prebaseline = GetSecs;
        % baseline fixation - currently 2 seconds = 1 TR
        while GetSecs < prebaseline + baseline_fixation
            CenterText(w,'+', white,0,0);
            Screen('TextSize',w, 60);
            Screen(w,'Flip');
        end
        


if use_eyetracker
    Eyelink('Message',Eventflag(GenFlags.Fixation.str,task,RunNum,1,anchor));
end


%   Reading everying from the sorted StopGo file - vars has everything
%---------------------------
%do it two times because there are two repetitions for each scanning
%run

TrialData1=table(Shuffle(1:40).',vars{1},vars{2},vars{3},vars{4},vars{5},AllRewards(:,RunNum),'VariableNames',{'shuff_ind','names','trialType','bidIndex','bidValues','BDMtrialIndex','rewards'});
TrialData1=sortrows(TrialData1);

TrialData2=table(Shuffle(1:40).',vars{1},vars{2},vars{3},vars{4},vars{5},AllRewards(:,RunNum),'VariableNames',{'shuff_ind','names','trialType','bidIndex','bidValues','BDMtrialIndex','rewards'});
TrialData2=sortrows(TrialData2);

TrialData=[TrialData1;TrialData2];


%	pre-allocating matrices
%---------------------------
Cue_time(1:length(TrialData.trialType),1) = 999;
respTime(1:length(TrialData.trialType),1) = 999;
respInTime(1:length(TrialData.trialType),1) = 999;
keyPressed(1:length(TrialData.trialType),1) = 999;

%   reading in images
%---------------------------
food_items = cell(1, size(TrialData,1));
for i = 1:size(TrialData,1)
    food_items{i} = imread(sprintf('stim/%s',TrialData.names{i}));
end

%   Read in info about ladders
% - - - - - - - - - - - - - - -

if  RunNum==1 % if this is the very first run of the experiment, start ladders at 750
    Ladder1IN=750;
    Ladder2IN=750;
else % read the ladders from the previous run's txt file -----------------------------------------
    last_run=dir([outputPath '/' subjectID '_training_run_' sprintf('%02d',RunNum-1) '*.txt']);
    clear last_run_fid last_run_data
    last_run_fid=fopen([outputPath,'/',last_run(end).name]);
    last_run_data=textscan(last_run_fid,'%s %f %f %s %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',1);
    last_run_fid=fclose(last_run_fid);
    
    Ladder1IN=last_run_data{12}(end);
    Ladder2IN=last_run_data{13}(end);
end

Ladder1(1,1) = Ladder1IN;
Ladder2(1,1) = Ladder2IN;

runStartTime = GetSecs - anchor;
if use_eyetracker
    Eyelink('Message', Eventflag(GenFlags.RunStart.str,task,RunNum,1,anchor)); % mark start time in file
end

%% -----------------------------------
%%   start trial loop
%% -----------------------------------

%   for trialNum = 1:6; % shorter version for debugging
for trialNum = 1:size(TrialData,1)   % To cover all the items in one run.
    
    
    
    CueNum=TrialData.rewards(trialNum);
    CueStrComb=[CueStr num2str(CueNum)];
    if CueType==1
        CueStrComb='**';
    end
    
    Screen('PutImage',w,food_items{trialNum});
    Screen('Flip',w,anchor+onsetlist(trialNum)+runStartTime); % display images according to Onset times
    image_start_time = GetSecs;
    actual_onset_time(trialNum,1) = image_start_time - anchor;
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        % messages to save on each trial ( trial number, onset and RT)
        Eyelink('Message', Eventflag(GenFlags.TrialStart.str,task,RunNum,trialNum,anchor)); % mark start time in file
    end
    
    noresp = 1;
    notone = 1;
    
    KbQueueFlush; % this is important for the queue to be empty for the next trial (to prevent RT with negative values)
    KbQueueStart;
    
    %---------------------------------------------------
    %% 'EVALUATE RESPONSE & ADJUST LADDER ACCORDINGLY'
    %---------------------------------------------------
    while (GetSecs-image_start_time < image_duration)
        
        %   look for response
        [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
        if pressed && noresp
            firstKeyPressed = find(firstPress==min(firstPress(firstPress>0)));
            
            if length(firstKeyPressed)>=2 % In case two keys were pressed together on the exact micro-sec
                firstKeyPressed = firstKeyPressed(1);
            end
            
            respTime(trialNum,1) = firstPress(firstKeyPressed)-image_start_time;
            if use_eyetracker
                %   Eyelink MSG
                % ---------------------------
                Eyelink('Message', Eventflag(GenFlags.Response.str,task,RunNum,trialNum,anchor)); % mark start time in file
            end
            
            tmp = KbName(firstKeyPressed);
            if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                tmp = char(tmp);
            end
            keyPressed(trialNum,1) = tmp(1);
        end
        
        %High-Valued BEEP items
        %---------------------------
        if  TrialData.trialType(trialNum) == 11 && (GetSecs - image_start_time >= Ladder1(length(Ladder1),1)/1000) && notone %  TrialData.trialType contains the information if a certain image is a GO/NOGO trial
            % Beep!
            
            notone = 0;
            Cue_time(trialNum,1) = GetSecs-image_start_time;
            Screen('PutImage',w,food_items{trialNum});
            % Draw the Cue image
            Screen('DrawTextures', w, imageTextureCue, [], dstRect, 0);
            Screen('TextSize',w, 60);
            CenterText(w,CueStrComb, white,0,0);
            Screen('Flip',w);
            
            if use_eyetracker
                %   Eyelink MSG
                % ---------------------------
                Eyelink('Message', Eventflag(GenFlags.CueStart.str,task,RunNum,trialNum,anchor)); % mark start time in file
            end
            
            % check responses according to the pre defined keys
                if keyPressed(trialNum,1) == blue || keyPressed(trialNum,1) == yellow
                    noresp = 0;
                    if respTime(trialNum,1) < Ladder1(length(Ladder1),1)/1000
                        respInTime(trialNum,1) = 11; %was a GO trial with HV item but responded before SS
                    else
                        respInTime(trialNum,1)= 110; %was a Go trial with HV item but responded after SS within 1000 msec
                    end
                end
       
            
            %Low-Valued BEEP items
            %---------------------------
        elseif   TrialData.trialType(trialNum) == 22 && (GetSecs - image_start_time >= Ladder2(length(Ladder2),1)/1000) && notone % TrialData.trialType contains the information if a certain image is a GO/NOGO trial
            % Beep!
            
            notone = 0;
            Cue_time(trialNum,1) = GetSecs-image_start_time;
            Screen('PutImage',w,food_items{trialNum});
            % Draw the Cue image
            Screen('DrawTextures', w, imageTextureCue, [], dstRect, 0);
            Screen('TextSize',w, 60);
            CenterText(w,CueStrComb, white,0,0);
            Screen('Flip',w);
            
            if use_eyetracker
                %   Eyelink MSG
                % ---------------------------
                Eyelink('Message', Eventflag(GenFlags.CueStart.str,task,RunNum,trialNum,anchor)); % mark start time in file
            end
            
          % check responses according to the pre defined keys
                if keyPressed(trialNum,1) == blue || keyPressed(trialNum,1) == yellow
                    noresp = 0;
                    if respTime(trialNum,1) < Ladder2(length(Ladder2),1)/1000
                        respInTime(trialNum,1) = 22; %was a GO trial with LV item but responded before SS
                    else
                        respInTime(trialNum,1) = 220; %was a Go trial with LV item but responded after SS within 1000 msec
                    end
                end
            
            %No-BEEP 
            %---------------------------
        elseif   mod( TrialData.trialType(trialNum),11) ~= 0 && noresp % these will now be the NOGO trials
            % check responses according to the pre defined keys
                if keyPressed(trialNum,1) == blue || keyPressed(trialNum,1) == yellow
                    noresp = 0;
                    if  TrialData.trialType(trialNum) == 12
                        respInTime(trialNum,1) = 12; % a stop trial but responded within 1000 msec HV item - not good but don't do anything
                    else
                        respInTime(trialNum,1) = 24; % a stop trial but responded within 1000 msec LV item - not good but don't do anything
                    end
                end
          
        end %evaluate trial_type
    end %%% End big while waiting for response within 1000 msec
    
    
    
    %   Show fixation - after pic fixation - fixed of 1 sec which is
    %   the minimum jitter time for fixation after pics.
    %---------------------------
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip', image_start_time+1);
    fix_time(trialNum,1) = GetSecs ;
    fixcrosstime = GetSecs;
    if use_eyetracker
        Eyelink('Message', Eventflag(GenFlags.Fixation.str,task,RunNum,trialNum,anchor)); % mark start time in file
    end
    
    
    if noresp == 1
        %---------------------------
        % these are additional 500msec to monitor responses
        
        while (GetSecs-fix_time(trialNum,1) < 0.5)
            
            %   look for response
            [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
            if pressed && noresp
                firstKeyPressed = find(firstPress==min(firstPress(firstPress>0)));
                if length(firstKeyPressed)>=2 % In case two keys were pressed together on the exact micro-sec
                    firstKeyPressed = firstKeyPressed(1);
                end
                respTime(trialNum,1) = firstPress(firstKeyPressed)-image_start_time;
                tmp = KbName(firstKeyPressed);
                if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                    tmp = char(tmp);
                end
                keyPressed(trialNum,1) = tmp(1);
                if use_eyetracker
                    Eyelink('Message', Eventflag(GenFlags.Response.str,task,RunNum,trialNum,anchor)); % mark start time in file
                end
            end
            
                if keyPressed(trialNum,1) == blue || keyPressed(trialNum,1) == yellow
                    noresp = 0;
                    switch TrialData.trialType(trialNum)
                        case 11
                            if respTime(trialNum,1) >= 1
                                respInTime(trialNum,1) = 1100; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                            elseif respTime(trialNum,1) < 1
                                respInTime(trialNum,1) = 110;
                            end
                        case 22
                            if respTime(trialNum,1) >= 1
                                respInTime(trialNum,1) = 2200; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                            elseif respTime(trialNum,1) < 1
                                respInTime(trialNum,1) = 220;
                            end
                        case 12
                            respInTime(trialNum,1) = 12; % a stop trial and responded after 1000 msec  HV item - don't touch
                        case 24
                            respInTime(trialNum,1) = 24; % % a stop trial and  responded after 1000 msec HV item - don't touch
                    end
                end
                
            
        end % End while of additional 500 msec
    else % the subject has already responded during the first 1000 ms
    end  % end if noresp
    
    %%	This is where its all decided !
    %---------------------------
    if noresp
        switch TrialData.trialType(trialNum)
            case 11
                respInTime(trialNum,1) = 1; %unsuccessful Go trial HV - didn't press a button at all - trial too hard - need to decrease ladder
            case 22
                respInTime(trialNum,1) = 2; % unsuccessful Go trial LV - didn't press a button at all - trial too hard - need to decrease ladder
            case 12
                respInTime(trialNum,1) = 120; % ok NOGO trial didn't respond after 1500 msec in NOGO trial HV
            case 24
                respInTime(trialNum,1) = 240; % ok NOGO trial didn't respond after 1500 msec in NOGO trial LV
        end
    end
    
    
    switch respInTime(trialNum,1)
        case 1 % didn't respond even after 1500 msec on HV GO trial - make it easier decrease SSD by step
            if (Ladder1(length(Ladder1),1)<0.001)
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1);
            else
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1)-Step;
            end;
            
        case 2 % didn't respond even after 1500 msec on LV GO trial - make it easier decrease SSD by step
            if (Ladder2(length(Ladder2),1)<0.001)
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1);
            else
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1)-Step;
            end;
            
            
        case 1100 %  responded after 1500 msec on HV GO trial - make it easier decrease SSD by step
            if (Ladder1(length(Ladder1),1)<0.001)
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1);
            else
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1)-Step;
            end;
            
        case 2200 %  responded after 1500 msec on LV GO trial - make it easier decrease SSD by step
            if (Ladder2(length(Ladder2),1)<0.001)
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1);
            else
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1)-Step;
            end;
            
            
            
        case 11
            if (Ladder1(length(Ladder1),1) > 910); %was a GO trial with HV item but responded before SS make it harder - increase SSD by Step/3
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1);
            else
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1)+Step/3;
            end;
            
        case 22
            if (Ladder2(length(Ladder2),1) > 910); %was a GO trial with LV item but responded before SS make it harder - - increase SSD by Step/3
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1);
            else
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1)+Step/3;
            end;
            
        case 110 % pressed after Go signal but below 1000 - - increase SSD by Step/3 - these are the good trials!
            if (Ladder1(length(Ladder1),1) > 910);
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1);
            else
                Ladder1(length(Ladder1)+1,1) = Ladder1(length(Ladder1),1)+Step/3;
            end;
            
        case 220 % pressed after Go signal but below 1000 - - increase SSD by Step/3 - these are the good trials!
            if (Ladder2(length(Ladder2),1) > 910);
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1);
            else
                Ladder2(length(Ladder2)+1,1) = Ladder2(length(Ladder2),1)+Step/3;
            end;
            
    end % end switch respInTime(trialNum,1)
    
    
    %   'Save data'
    %---------------------------
    
    fprintf(fid1,'%s\t %d\t %d\t %s\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %.2f\t %.2f\t %.2f\t %d\t %.2f\t \n', subjectID, order, RunNum, TrialData.names{trialNum}, actual_onset_time(trialNum,1), TrialData.trialType(trialNum), respTime(trialNum,1)*1000, respInTime(trialNum,1), Cue_time(trialNum,1)*1000, keyPressed(trialNum,1),   fix_time(trialNum,1)-anchor, Ladder1(length(Ladder1)), Ladder2(length(Ladder2)), TrialData.bidIndex(trialNum), TrialData.BDMtrialIndex(trialNum,1),TrialData.bidValues(trialNum,1));
    
end; %	End the big trialNum loop showing all the images in one run.



KbQueueFlush;

Ladder1end = Ladder1(length(Ladder1));
Ladder2end = Ladder2(length(Ladder2));
correct(1) = 0;
% Correct trials are when the subject pressed the button on a go trial,
% either before (11,22) or after (110,220)the beep (but before the
% image disappeared)
correct(1) = length(find(respInTime == 11 | respInTime == 110 | respInTime == 22 | respInTime == 220 ));
numGoTrials(RunNum) = length(find(TrialData.trialType == 11 | TrialData.trialType == 22));
mean_RT = mean(respTime(respInTime == 110 | respInTime == 220));


if use_eyetracker
    Eyelink('Message', Eventflag(GenFlags.RunEnd.str,task,RunNum,trialNum,anchor)); % mark start time in file
end


%% after run fixation - in imaging, to capture the hrf of the last trial.
%because there are 2 secs after image anyway - add 4 to achive 6 sec for
%hrf
postexperiment = GetSecs;

if use_eyetracker
    Eyelink('Message', Eventflag(GenFlags.Fixation.str,task,RunNum,trialNum,anchor)); % mark start time in file
end

while GetSecs < postexperiment+afterrunfixation;
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
end


fclose(fid1);

toc


%---------------------------------------------------------------
%%   save data to a .mat file & close out
%---------------------------------------------------------------
outfile = strcat(outputPath, '/', subjectID,'_training_run', sprintf('%02d',RunNum),'_', timestamp,'.mat');
% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
clear food_items Instructions*;

save(outfile);


%%   Finishing eye tracking  %
%---------------------------------------------------------------
if use_eyetracker
    
    %---------------------------
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message',['Eyetracking_closeTime: ',num2str(GetSecs-anchor)]);
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
        movefile(edfFile,['./Output/', subjectID,'_',task,'_run' ,sprintf('%02d',RunNum),'_', timestamp,'.edf']);
    end;
end





%   outgoing msg & closing
% ------------------------------
if RunNum == total_number_of_runs % if this is not the last run
    CenterText(w,'Great Job. Thank you!',Green, 0,-270);
    Screen('Flip',w);
end

% clc;
WaitSecs(3);
KbQueueFlush;
Screen('CloseAll');
ShowCursor;

if RunNum<total_number_of_runs
    fprintf(['\nyour next run is: ' num2str(RunNum+1),'\n\n']);
else
    disp('continue to next part');
end

% clear all
sca;
