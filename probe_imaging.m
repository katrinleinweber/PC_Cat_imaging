% function probe(subjectID, order, mainPath, isMRI, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun, use_eyetracker)

% function probe_Israel(subjectID, order, mainPath, isMRI, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun, use_eyetracker)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik May 2015 ====================
% edited on March 2017 to include eye-tracking
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the probe session of the boost (cue-approach) task. In
% this session, the subject is presented with two items in each trial, and
% should choose which one he prefers within 1.5 seconds.
% This function runs one run each time. Each block is devided to
% 'numRunsPerBlock' runs. The stimuli for each run are shuffled, chosen and
% organized in the 'organizeProbe_Israel' function.

% This function is a version in which only the 40 of the items are included
% in the training


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''stopGoList_allstim_order*.txt'' --> created by sortBDM_Israel
%   ''stimuliForProbe_order%d_block_%d_run%d.txt'' --> Created by
%   organizeProbe_Israel
%   onset lists


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''probe_block_' block '_' timestamp '.txt''
%


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'Brev_snacks_999';
% order = 1;
% subjectID =  'Brev_snacks_998';
% order = 2;
% isMRI = 0;
% sessionNum = 1;
% mainPath = pwd;
% block = 1;
% numRun = 1;
% trialsPerRun = 8; % for debugging

clear all
Screen('Preference', 'SkipSyncTests', 1);

rng shuffle
% Screen('Preference', 'SkipSyncTests', 1);

% =========================================================================
% Get input args and check if input is ok
% =========================================================================

experimentName='cat_between_MRI';
oksessionNum = [1 2 3 4 5];

% input checkers
subjectID_num = input('Subject number (only the digits): ');
while isempty(subjectID_num)
    disp('ERROR: no value entered. Please try again.');
    subjectID_num = input('Subject number (only the digits):');
end
subjectID = [experimentName '_' num2str(subjectID_num)];


% Assign order
% --------------------------
% give order value of '1' or '2' for subjects with odd or even ID, respectively
if mod(subjectID_num,2) == 1 % subject code is odd
    order = 1;
else % subject code is even
    order = 2;
end



% =========================================================================
% set the computer and path
% =========================================================================
isMRI=1; % 1-MRI experiment
% Set main path
mainPath=pwd;
outputPath = [mainPath '/Output'];
numRunsPerBlock = 2;
sessionNum=1;


%==============================================
%% 'GLOBAL VARIABLES'
%==============================================

% 'block' indicates how many times the subject has completed all
% trials of the probe experiment. Therefore, the first time the subject completes
% a probe block on his/her 2nd session, it is actually probe block 3 for
% that person:

if isMRI == 1 % If it's an fMRI experiment, let the experimenter insert the block number, to control when to start
    okblock = [1 2];
    which_block = input('Enter Block number (1 or 2): ');
    while isempty(which_block) || sum(okblock == which_block) ~=1
        disp('ERROR: input must be 1 or 2. Please try again.');
        which_block = input('Enter Block (1 or 2): ');
    end
    block = (sessionNum-1)*2 + which_block;
    
    %     okRun = [1 2];
    %     which_run = input('Enter run number (1 or 2): ');
    %     while isempty(which_run) || sum(okRun == which_run) ~=1
    %         disp('ERROR: input must be 1 or 2. Please try again.');
    %         which_run = input('Enter run number (1 or 2): ');
    %     end
    %     numRun = which_run;
    numRun=1;
    
    okEyetracker = [1 0];
    ask_if_want_eyetracker = input('Do you want eyetracking (1 - yes, 0 - no): ');
    while isempty(ask_if_want_eyetracker) || sum(okEyetracker == ask_if_want_eyetracker) ~=1
        disp('ERROR: input must be 1 or 0. Please try again.');
        ask_if_want_eyetracker = input('Do you want eyetracking (1 - yes, 0 - no): ');
    end
    use_eyetracker=ask_if_want_eyetracker; % set to 1/0 to turn on/off eyetracker functions
    
    
    
    fprintf(['\n\nsubjectID is: ' num2str(subjectID) '\n']);
    disp(['sessionNum is: ' num2str(sessionNum)]);
    disp(['block is: ' num2str(block)]);
    %     disp(['numRun is: ' num2str(numRun)]);
    disp(['use eyetracker is: ' num2str(use_eyetracker)]);
    fprintf('\n')
    
    GoOn=input('are all variables ok? (1-yes, 0-no) ');
    if GoOn==0
        error('please check you numbers and start again')
    end
    
end

%   'timestamp'
% - - - - - - - - - - - - - - - - -
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

%   'set phase times'
% - - - - - - - - - - - - - - - - -
maxtime = 1.5;      % 1.5 second limit on each selection
baseline_fixation_dur = 2; % Need to modify based on if first few volumes are saved or not
afterrunfixation = 6;



% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
Instructions = dir([mainPath '/Instructions/probe.JPG' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);
Trigger = dir([mainPath '/Instructions/Trigger_image.png' ]);
Trigger_image = imread([mainPath '/Instructions/' Trigger(1).name]);
Continue = dir([mainPath '/Instructions/Continue.jpg' ]);
Continue_image = imread([mainPath '/Instructions/' Continue(1).name]);
%==============================================
%% 'INITIALIZE Screen variables'
%==============================================
Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize=32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize/4);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

HideCursor;


% Define Colors
% - - - - - - - - - - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
green = [0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);


% setting the text
% - - - - - - - - - - - - - - -
theFont = 'Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

% stack locations
% - - - - - - - - - - - - - - -
[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;

stackW = 576;
stackH = 432;

leftRect = [xcenter-stackW-300 ycenter-stackH/2 xcenter-300 ycenter+stackH/2];
rightRect = [xcenter+300 ycenter-stackH/2 xcenter+stackW+300 ycenter+stackH/2];

penWidth = 10;

HideCursor;

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
    task=GenFlags.Probe.str;
    edfFile='Probe.edf';
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

%==============================================
%% 'ASSIGN response keys'
%==============================================
KbName('UnifyKeyNames');

if isMRI == 1
    leftstack = 'b';
    rightstack = 'y';
    badresp = 'x';
else
    leftstack = 'u';
    rightstack = 'i';
    badresp = 'x';
end

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_stopGoList_allstim_order*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
stimName = data{1};
% bidIndex = data{3};
bidValue = data{4};
fclose(fid);

%   'read in organized list of stimuli for this run'
% - - - - - - - - - - - - - - - - - - - - - - - - - -

fid = fopen([outputPath '/' sprintf('%s_stimuliForProbe_order%d_block_%d_run%d.txt',subjectID,order,block,numRun)]);
stimuli = textscan(fid, '%d %d %d %d %s %s') ;% these contain everything from the organizedProbe
stimnum1 = stimuli{1};
stimnum2 = stimuli{2};
leftGo = stimuli{3};
pairType = stimuli{4};
leftname = stimuli{5};
rightname = stimuli{6};
fclose(fid);

trialsPerRun = size(stimuli{1},1);

%   'load image arrays'
% - - - - - - - - - - - - - - -
food_items = cell(1,length(stimName));
for i = 1:length(stimName)
    food_items{i} = imread([mainPath sprintf('/stim/%s',stimName{i})]);
end


%   'load onsets'
% - - - - - - - - - - - - - - -
r = Shuffle(1:4);
onsetlist = load([mainPath '/Onset_files/probe_onset_length_72_' num2str(r(1)) '.mat']);
onsetlist = onsetlist.onsetlist;

%-----------------------------------------------------------------
%% 'Write output file header'
%-----------------------------------------------------------------

fid1 = fopen([outputPath '/' subjectID sprintf('_probe_block_%02d_run%d_', block, numRun) timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\tscanner\torder\tblock\trun\ttrial\tonsettime\tImageLeft\tImageRight\tbidIndexLeft\tbidIndexRight\tIsleftGo\tResponse\tPairType\tOutcome\tRT\tbidLeft\tbidRight\t\n'); %write the header line


%==============================================
%% 'Display Main Instructions'
%==============================================


KbQueueCreate;
Screen('TextSize',w, 40);

if block ==1 && numRun==1   % if this is the first block and first run, present full instructions
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
    
    % escapeKey = KbName('space');
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


anchor=GetSecs;

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    WaitSecs(.05);
end

%   baseline fixation cross
% - - - - - - - - - - - - -
prebaseline = GetSecs;
% baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
while GetSecs < prebaseline+baseline_fixation_dur
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
end

postbaseline = GetSecs;
baseline_fixation = postbaseline - prebaseline;


%==============================================
%% 'Run Trials'
%==============================================

runtrial = 1;
runStart = GetSecs;

if use_eyetracker
    Eyelink('Message',Eventflag(GenFlags.RunStart.str,task,numRun,1,runStart));
end

%for trial = 1:6 % for debugging
for trial = 1:trialsPerRun
    
    
    % initial box outline colors
    % - - - - - - -
    colorLeft = black;
    colorRight = black;
    out = 999;
    
    
    %-----------------------------------------------------------------
    % display images
    %-----------------------------------------------------------------
    if leftGo(trial) == 1
        Screen('PutImage',w,food_items{stimnum1(trial)}, leftRect);
        Screen('PutImage',w,food_items{stimnum2(trial)}, rightRect);
    else
        Screen('PutImage',w,food_items{stimnum2(trial)}, leftRect);
        Screen('PutImage',w,food_items{stimnum1(trial)}, rightRect);
    end
    CenterText(w,'+', white,0,0);
    StimOnset = Screen(w,'Flip', runStart+onsetlist(runtrial)+baseline_fixation);
    
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.TrialStart.str,task,numRun,trial,runStart));
    end
    
    KbQueueFlush;
    KbQueueStart;
    
    
    %-----------------------------------------------------------------
    % get response
    %-----------------------------------------------------------------
    
    noresp = 1;
    goodresp = 0;
    while noresp
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck;
        
        if keyIsDown && noresp
            keyPressed = KbName(firstPress);
            if use_eyetracker
                Eyelink('Message',Eventflag(GenFlags.Response.str,task,numRun,trial,runStart));
                
            end
            if ischar(keyPressed) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed = char(keyPressed);
                keyPressed = keyPressed(1);
            end
            switch keyPressed
                case leftstack
                    respTime = firstPress(KbName(leftstack))-StimOnset;
                    noresp = 0;
                    goodresp = 1;
                case rightstack
                    respTime = firstPress(KbName(rightstack))-StimOnset;
                    noresp = 0;
                    goodresp = 1;
            end
        end % end if keyIsDown && noresp
        
        
        % check for reaching time limit
        if noresp && GetSecs-runStart >= onsetlist(runtrial)+baseline_fixation+maxtime
            noresp = 0;
            keyPressed = badresp;
            respTime = maxtime;
        end
    end % end while noresp
    
    %-----------------------------------------------------------------
    % determine what bid to highlight
    %-----------------------------------------------------------------
    
    switch keyPressed
        case leftstack
            colorLeft = green;
            if leftGo(trial) == 0
                out = 0;
            else
                out = 1;
            end
        case rightstack
            colorRight = green;
            if leftGo(trial) == 1
                out = 0;
            else
                out = 1;
            end
    end
    
    if goodresp==1
        if leftGo(trial)==1
            Screen('PutImage',w,food_items{stimnum1(trial)}, leftRect);
            Screen('PutImage',w,food_items{stimnum2(trial)}, rightRect);
        else
            Screen('PutImage',w,food_items{stimnum2(trial)}, leftRect);
            Screen('PutImage',w,food_items{stimnum1(trial)}, rightRect);
        end
        Screen('FrameRect', w, colorLeft, leftRect, penWidth);
        Screen('FrameRect', w, colorRight, rightRect, penWidth);
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip',runStart+onsetlist(trial)+respTime+baseline_fixation);
        
    else
        Screen('DrawText', w, 'You must respond faster!', xcenter-400, ycenter, white);
        Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+baseline_fixation);
        if use_eyetracker
            Eyelink('Message',Eventflag(GenFlags.RespondFaster.str,task,numRun,trial,runStart));
        end
    end % end if goodresp==1
    
    
    %-----------------------------------------------------------------
    % show fixation ITI
    %-----------------------------------------------------------------
    Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+.5+baseline_fixation);
    
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.Fixation.str,task,numRun,trial,runStart));
        
    end
    
    if goodresp ~= 1
        respTime = 999;
    end
    
    %-----------------------------------------------------------------
    % 'Save data'
    %-----------------------------------------------------------------
    if leftGo(trial)==1
        fprintf(fid1,'%s\t %d\t %d\t %s\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t %d\t %s\t %d\t %d\t %.2f\t %.2f\t %.2f\t \n', subjectID, isMRI, order, sprintf('%02d', block), numRun, runtrial, StimOnset-anchor, char(leftname(trial)), char(rightname(trial)), stimnum1(trial), stimnum2(trial), leftGo(trial), keyPressed, pairType(trial), out, respTime*1000, bidValue(stimnum1(trial)), bidValue(stimnum2(trial)));
    else
        fprintf(fid1,'%s\t %d\t %d\t %s\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t %d\t %s\t %d\t %d\t %.2f\t %.2f\t %.2f\t \n', subjectID, isMRI, order, sprintf('%02d', block), numRun, runtrial, StimOnset-anchor, char(leftname(trial)), char(rightname(trial)), stimnum2(trial), stimnum1(trial), leftGo(trial), keyPressed, pairType(trial), out, respTime*1000, bidValue(stimnum2(trial)), bidValue(stimnum1(trial)));
    end
    
    runtrial = runtrial+1;
    %     KbQueueFlush;
    
end % loop through trials
fclose(fid1);


Postexperiment = GetSecs;

while GetSecs < Postexperiment + afterrunfixation;
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
    
end

toc

if use_eyetracker
    Eyelink('Message',Eventflag(GenFlags.RunEnd.str,task,numRun,trial,runStart));
    
end


%-----------------------------------------------------------------
%	display outgoing message
%-----------------------------------------------------------------
Screen('FillRect', w, black);
Screen('TextSize',w, 40);

if block == 2 &&  numRun == numRunsPerBlock
    CenterText(w,sprintf('Thank you!') ,white,0,-270);
    Screen('Flip',w);
end % end if mod(block,2) == 1

Screen('CloseAll');
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
    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message', ['Eyetracking close time ',num2str(GetSecs-runStart)]); % mark start time in file
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
        movefile(edfFile,['./Output/', subjectID,'_Probe_eyetracking_block_',num2str(block),'_run_',num2str(numRun),'_' timestamp,'.edf']);
    end;
end


%---------------------------------------------------------------
% create a data structure with info about the run
%---------------------------------------------------------------
outfile = strcat(outputPath,'/', sprintf('%s_probe_block_%2d_run_%2d_%s.mat',subjectID,block,numRun,timestamp));

% create a data structure with info about the run
run_info.subject=subjectID;
run_info.date=date;
run_info.outfile=outfile;

run_info.script_name=mfilename;
clear food_items Instruction*;
save(outfile);

KbQueueFlush;


if block==1
    fprintf('\nyour next block is: 2 \n');
    %     fprintf(['\nyour next run is: 2' ,'\n']);
else
    
    disp('continue to next part');
    
end
%end % end function

