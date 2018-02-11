function BDM_Snacks(subjectID,sessionNum)
%=========================================================================
% BDM task 
%=========================================================================
rng shuffle
Screen('Preference', 'SkipSyncTests', 1);

tic
c=clock;
hr=num2str(c(4));
minute=num2str(c(5));
timestamp=[date,'_',hr,'h',minute,'m'];


%--------------------------------------------------------------------------
%% Locations, File Types & Names PARAMETERS
%--------------------------------------------------------------------------
% Stimuli:
StimLocation = './stim/';
StimFileType = 'bmp';

% for output file
OutputFolder = 'Output/';
StimuliName = 'snacks';

%--------------------------------------------------------------------------
%% PARAMETERS FOR THE RANKING AXIS
%--------------------------------------------------------------------------
% Parameters for RANKING RANGE
RankingMin = 0;
RankingMax = 10;

% Parameters for Creatining the RANKING AXIS:
RelativeSizeOfRankingAxisFromScreenSize = 1/4;
YaxisRelativeMiddlePoint = 0.9;
YaxisWidthFactor = 0.0102;

% Parameters for the MOVING INDICATOR on the ranking axis:
penWidth = 3;
AdditionalMovingIndicatorLengthFactor = 0.018; % Extension from Each side of the ranking axis.

% Parameters for FIGURES PRESENTATION:
TextSizeForFiguresOnAxis = 30;
DistanceOfFiguresFactor = 0.0065;

% Fixation cross fix:
FixForFixationCrossLocation = -33.5; % A fix for fixation cross on center text to be in center on the Y axis. Relevant for text size 60.

%---------------------------------------------------------------
%% 'INITIALIZE Screen variables'
%---------------------------------------------------------------

Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
screennum = min(Screen('Screens'));

pixelSize=32;
%[w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

% Here Be Colors
black=BlackIndex(w); % Should equal 0.
white=WhiteIndex(w); % Should equal 255.
green=[0 255 0];


% set up Screen positions for stimuli
[wWidth, wHeight]=Screen('WindowSize', w);
xcenter=wWidth/2;
ycenter=wHeight/2;

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% text stuffs
theFont='Arial';
Screen('TextFont',w,theFont);

instrSZ=40;
betsz=60;

%--------------------------------------------------------------------------
%% SETTINGS FOR THE RANKING AXIS
%--------------------------------------------------------------------------
% Settings for Creatining the RANKING AXIS:
AxisFromX = wWidth*RelativeSizeOfRankingAxisFromScreenSize; % Starting point on X axis
AxisToX = wWidth*(1-RelativeSizeOfRankingAxisFromScreenSize); % Ending point on X axis
AxisFromY = round(wHeight * (YaxisRelativeMiddlePoint - (YaxisRelativeMiddlePoint * YaxisWidthFactor))); % Starting point on Y axis
AxisToY = round(wHeight * (YaxisRelativeMiddlePoint + (YaxisRelativeMiddlePoint * YaxisWidthFactor))); % Ending point on Y axis

% Settings for the MOVING INDICATOR on the ranking axis:
CenterOfMovingIndicator = mean([AxisFromY AxisToY]);
AdditionToYAxisFromEachSide = wHeight*AdditionalMovingIndicatorLengthFactor;

% Settings for FIGURES PRESENTATION:
RankingIntegers = RankingMax - RankingMin + 1;
SpotsForIndicatorsOnAxis = linspace(AxisFromX, AxisToX, RankingIntegers);
DistanceOfFiguresFromAxis = round(wHeight * DistanceOfFiguresFactor) + AdditionToYAxisFromEachSide ;
FixForFiguresOnXaxis = round(wHeight * 0.0074);

%---------------------------------------------------------------
%% 'LOAD image arrays'
%---------------------------------------------------------------
stimuli_images=dir([StimLocation '*.' StimFileType]);

shuffledlist=Shuffle(1:length(stimuli_images));
imageArrays=cell(length(shuffledlist),1);
for i=1:length(shuffledlist)
    imageArrays{i}=imread([StimLocation stimuli_images(shuffledlist(i)).name]);
    ImageSize=size(imageArrays{i});
    RectArray{i}= CenterRectOnPointd([ 0 0 ImageSize(2) ImageSize(1)], xcenter, ycenter);

end


%---------------------------------------------------------------
%% 'Write output file header'
%---------------------------------------------------------------
fid1=fopen([OutputFolder subjectID '_' StimuliName '_BDM' num2str(sessionNum) '_' timestamp '.txt'], 'w');
fprintf(fid1,'subjectID\truntrial\tonsettime\tName\tBid\tRT\n'); %write the header line

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

Screen('TextSize',w, instrSZ);
CenterText(w,'Press any key to start', green,0,60);
HideCursor;
Screen('Flip', w);
WaitSecs(0.05); % prevent key spillover
noresp=1;
while noresp
    [keyIsDown] = KbCheck;
    if keyIsDown && noresp
        noresp = 0;
    end
end

Screen('TextSize',w, betsz);
CenterText(w,'+', white,0,FixForFixationCrossLocation);
Screen(w,'Flip');
WaitSecs(0.3);


%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------
runStart=GetSecs;

%make it so that user can press cntrl+c and end execution of function
KbQueueCreate;
KbQueueFlush;
KbQueueStart;

% for trial=1:5 % debugging
for trial=1:length(stimuli_images)
    
    bid=[];
    noresp=1;
    Screen('TextSize',w,TextSizeForFiguresOnAxis);
    eventTime=[];
    ShowCursor;
    SetMouse(xcenter,ycenter);
    while noresp
        % Track cursor movement and check for response
        [CurrentX,CurrentY,buttons] = GetMouse(w);
        if CurrentX >= AxisFromX && CurrentX <= AxisToX && CurrentY >= AxisFromY - AdditionToYAxisFromEachSide && CurrentY <= AxisToY + AdditionToYAxisFromEachSide
            Screen('PutImage',w,imageArrays{trial},RectArray{trial});
            Screen('FillRect', w ,[211 211 211] ,[AxisFromX, AxisFromY,  AxisToX, AxisToY]);
            for i = 1:length(SpotsForIndicatorsOnAxis) 
                DrawFormattedText(w, num2str(i-1), SpotsForIndicatorsOnAxis(i)-FixForFiguresOnXaxis, CenterOfMovingIndicator+DistanceOfFiguresFromAxis, [255 255 255]);
            end
            Screen('DrawLine', w ,[0 0 255], CurrentX, CenterOfMovingIndicator+AdditionToYAxisFromEachSide, CurrentX, CenterOfMovingIndicator-AdditionToYAxisFromEachSide ,penWidth);
            Screen(w,'Flip');
            if buttons(1) == 1
                bid = (CurrentX - AxisFromX) / (AxisToX - AxisFromX) * (RankingMax - RankingMin) + RankingMin; % Number of pixels from X axis beggining / Length of the axis * Units + Beggining of units.
                respTime = GetSecs - runStart - eventTime;
                noresp = 0;
                while any(buttons) % wait for release
                    [~,~,buttons] = GetMouse;
                end
                    

                 
            end
        else
            Screen('PutImage',w,imageArrays{trial},RectArray{trial});
            Screen('FillRect', w ,[211 211 211] ,[AxisFromX, AxisFromY,  AxisToX, AxisToY]);
            for i = 1:length(SpotsForIndicatorsOnAxis)
                DrawFormattedText(w, num2str(i-1), SpotsForIndicatorsOnAxis(i)- FixForFiguresOnXaxis, CenterOfMovingIndicator+DistanceOfFiguresFromAxis, [255 255 255]);
            end
            Screen(w,'Flip');
            if isempty(eventTime) % recording the presentation start time
                eventTime = GetSecs-runStart;
            end
        end
    end
    
    %-----------------------------------------------------------------
    % show fixation ITI
    Screen('TextSize',w, betsz);
    CenterText(w,'+', white,0,FixForFixationCrossLocation);
    Screen(w,'Flip');
    
     
    WaitSecs(0.3);
    
    %-----------------------------------------------------------------
    % write to output file
    
    fprintf(fid1,'%s\t%d\t%d\t%s\t%d\t%d\n', subjectID, trial, eventTime, stimuli_images(shuffledlist(trial)).name, bid, respTime);
end


HideCursor;

    CenterText(w,'Thank you!', green,0,20);
    Screen('Flip', w, 0,1);

WaitSecs(3); % prevent key spillover

fclose(fid1);
toc

ShowCursor;
Screen('closeall');


end %end function

