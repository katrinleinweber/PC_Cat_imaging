% function [averageDiff, maximumDiff, minimumDiff, diff] = createOnsetList(mean_t,min_t,max_t,interval,requiredLength,numTrials,numOnsetlists)
% -------------------------------------
% Created by Rotem Botvinik 4/12/14
% Edited by Rotem Botvinik 11/12/14
% -------------------------------------

% function [averageDiff, maximumDiff, minimumDiff, diff] = createOnsetList(mean_t,min_t,max_t,interval,requiredLength,numTrials,numOnsetlists)
% Creates and saves numOnsetlists onsetlists for the probe of the cue approach task.
% This function makes sure that the mean of diff (variable averageDiff) is close enough (+- 2%) to
% the requested mean value and that the length of the onsetlist (that is,
% the last onset time) is the same in all the lists (as requested in the
% requiredLength variable).


% --------
% Inputs:
% -------
% mean_t: the mean requested jitter (not including the 2 seconds of stimuli
% presentation)
% min_t: the minimum requested jitter (not including the 2 seconds of stimuli
% presentation)
% max_t: the maximum requested jitter (not including the 2 seconds of stimuli
% presentation)
% interval: the interval to ceil to (1 for integers)
% requiredLength: the requested length of all onsets, meaning the value of
% the last onset time
% numTrials: how many onset value should be in the list (what is the length
% of each onsetlist)
% numOnsetlists: how many onsetlists should the function create
% StimTime- time of stimuli presentation

% --------
% Outputs:
% --------
% averageDiff: the average jitter before adding the two seconds of the
% stimuli presentation
% maximumDiff: the maximum jitter before adding the two seconds of the
% stimuli presentation
% minimumDiff: the minimum jitter before adding the two seconds of the
% stimuli presentation
% diff: the jitter arrays after adding the two seconds of the
% stimuli presentation


% ----------------------------------
% Creating onsetlists for the probe:
% ----------------------------------
% For creaing onsetlists for the probe, we used jitter mean 3, min 1, max 12, interval 1; numTrials 68 (for sum of 136 comparisons), 
% requiredLength 336 (5*(68-1)+1 for pair number ---> mean length of one trial*(number of trials per run-1)+1 for paired number) (and then when adding 2 to the diff matrix we get mean 5, min 3, max 14, interval 1),
% lengthOnsetlist should match the number of comparisons in the probe / numRunsPerBlock (for example, 136/2=68)

% %for resp:
StimTime=2; mean_t=5; min_t=1; max_t=12; interval=1;  numTrials=20; numOnsetlists=8; 
requiredLength=(mean_t+StimTime)*(numTrials-1)+1; %5*(64-1)+1% mean leangth of trial*(number of trials per run-1)+1
task='resp';

% %for training:
% StimTime=1; mean_t=2; min_t=1; max_t=12; interval=1;  numTrials=80; numOnsetlists=6; 
% requiredLength=(mean_t+StimTime)*(numTrials-1)+1; %5*(64-1)+1% mean leangth of trial*(number of trials per run-1)+1
% task='train'; 
% 
% %for probe:
% StimTime=2; mean_t=3; min_t=1; max_t=12; interval=1;  numTrials=72; numOnsetlists=4; 
% requiredLength=(mean_t+StimTime)*(numTrials-1)+1; %5*(64-1)+1% mean leangth of trial*(number of trials per run-1)+1
% task='probe';

%% 
diff = zeros(numOnsetlists,numTrials-1);
maximumDiff = zeros(1,numOnsetlists);
minimumDiff = zeros(1,numOnsetlists);
averageDiff = zeros(1,numOnsetlists);

for listNum = 1:numOnsetlists
    onsetlist = zeros(1,numTrials); 
    while onsetlist(end) ~= requiredLength

         %     while  onsetlist(end) < requiredLength*0.99 || onsetlist(end) > requiredLength*1.01
    averageDiff(listNum) = 0;   
%         while averageDiff(listNum) ~= mean_t
        while averageDiff(listNum) < mean_t*0.98 || averageDiff(listNum) > mean_t*1.02
            diff(listNum,:) = zeros(1,numTrials-1);
            for i = 1:numTrials-1
                diff(listNum,i) = expsample(mean_t,min_t,max_t,interval);
            end % end for i = 1:numTrials-1
            maximumDiff(listNum) = max(diff(listNum,:));
            minimumDiff(listNum) = min(diff(listNum,:));
            averageDiff(listNum) = mean(diff(listNum,:));
            display(averageDiff(listNum))
        end % end while averageDiff(listNum) < mean_t*0.98 || averageDiff(listNum) > mean_t*1.02
        
        diff(listNum,:) = diff(listNum,:)+StimTime;
        
        for i = 2:numTrials
            onsetlist(i)=onsetlist(i-1)+diff(listNum,i-1);
        end % end for i = 2:numTrials
        
    end % end while onsetlist(end) ~= requiredLength
    
    
    save([task '_onset_length_' num2str(numTrials) '_' num2str(listNum)],'onsetlist');
    
end % end for listNum = 1:numOnsetlists


% end % end function