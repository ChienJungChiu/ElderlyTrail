function [results] = PTBDMSV2(stimulus,retention,probe,ITI,numSquares,numTrials,filename,XPixels,HBin,trackNum,DMSmode)
%UNTITLED3 Summary of this function goes here
Screen('Preference', 'SkipSyncTests', 0);


% timestamp pixel
horizontal_pix=1600;
vertical_pixel=1040;

folderNumber=1;
while 1
    DMSfolder = [filename '\DMS_' num2str(folderNumber)];
    if ~exist(DMSfolder)
        mkdir(DMSfolder)
        break
    end
    folderNumber=folderNumber+1;
end
filename = DMSfolder;
% stimulus = str2num(answer{1});
% retention = str2num(answer{2});
% probe = str2num(answer{3});
% ITI = str2num(answer{4});
% numSquares = str2num(answer{5});
results.numTrials = numTrials;
%% Get baseline and recovery time
if DMSmode == 2
    prompt = {'Baseline time (s)','Recovery time (s) [0 for no recovery]'};
    dlgtitle = 'DMS';
    dims = [1 40];
    definput = {'120','0'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if isempty(answer)
        baselineTime = 0;
        recoveryTime = 0;
    else
        baselineTime = str2double(answer{1});
        recoveryTime = str2double(answer{2});
    end
end
%%
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen number
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
Screen('Preference', 'SkipSyncTests', 1);
% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Make a base Rect of 200 by 200 pixels
dim = screenXpixels/(2*numSquares + 3);
%dim = 128;
baseRect = [0 0 dim dim];

% Make the coordinates for our grid of squares
gridCoords = (numSquares/2)-0.5;
[xPos, yPos] = meshgrid(-gridCoords:1:gridCoords, -gridCoords:1:gridCoords);

% Calculate the number of squares and reshape the matrices of coordinates
% into a vector
[s1, s2] = size(xPos);
numSquares = s1 * s2;
results.numSquares = numSquares;
xPos = reshape(xPos, 1, numSquares);
yPos = reshape(yPos, 1, numSquares);

% Scale the grid spacing to the size of our squares and centre
xPosLeft = xPos .* dim + screenXpixels * 0.25;
yPosLeft = yPos .* dim + yCenter;

xPosRight = xPos .* dim + screenXpixels * 0.75;
yPosRight = yPos .* dim + yCenter;

xPosCenter = xPos .* dim + screenXpixels * 0.50;
yPosCenter = yPos .* dim + yCenter;

% Make our rectangle coordinates
allRectsLeft = nan(4, 3);
allRectsRight = nan(4, 3);
allRectsCenter = nan(4,3);
for i = 1:numSquares
    allRectsLeft(:, i) = CenterRectOnPointd(baseRect,...
        xPosLeft(i), yPosLeft(i));
    allRectsRight(:, i) = CenterRectOnPointd(baseRect,...
        xPosRight(i), yPosRight(i));
    allRectsCenter(:, i) = CenterRectOnPointd(baseRect,...
        xPosCenter(i), yPosCenter(i));
end

% Here we set the size of the arms of our fixation cross
fixCrossDimPix = 40;

% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
crossxCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
crossyCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
crossallCoords = [crossxCoords; crossyCoords];

% Set the line width for our fixation cross
lineWidthPix = 4;

%%
% Begin DMS
Screen('FillRect',window,[0.5,0.5,0.5]);
Screen('Flip',window);
disp('Grey Screen: press anykey')
KbStrokeWait;
Screen('FillRect',window,[0,0,0]);

if DMSmode == 2 && baselineTime ~= 0%baseline
    Screen('DrawLines', window, crossallCoords,...
        lineWidthPix, white, [xCenter yCenter]);
    disp('Baseline: displaying cross')
    Screen('Flip',window);
    [retDMS] = StartAcquisition(); % Begin acquisition
    CheckWarning(retDMS);
    
    % add time stamps (modified by TY)
    tic;
    Screen('TextSize', window, 40);
    t=timer('ExecutionMode','fixedRate','Period',1,'TimerFcn',@timerCallback);
    start(t);
    pause(baselineTime);
%         WaitSecs(baselineTime);
    stop(t);
    delete(t);

    % Stop acquisition
    [retDMS] = AbortAcquisition();
    CheckWarning(retDMS);
    % Retrieve Images
    [retDMS, frameCount] = GetTotalNumberImagesAcquired();
    CheckWarning(retDMS);
    [I] = GetImagesFromCamera(XPixels,HBin,trackNum,frameCount);
    save([filename '\DMSbaseline.mat'],'I');
    Screen('FillRect',window,[0.5,0.5,0.5]);
    Screen('Flip',window);
    disp('Finished: press anykey')
    KbStrokeWait;
    Screen('FillRect',window,[0,0,0]);
end
results.TaskStartTime = clock;
DMSStartTime = GetSecs;

tic;
for trial = 1:results.numTrials
    %      % begin start trails  Pan Modified
    %     Screen('FillRect',window,[0,0,0]);
    %     Screen('DrawLines', window, crossallCoords,...
    %         lineWidthPix, white, [xCenter yCenter]);
    %     Screen('Flip',window);
    %     disp(['Start Trail ' num2str(trial)])
    
%     Screen('FillRect',window,[0,0,0]);
    Screen('DrawLines', window, crossallCoords,...
        lineWidthPix, white, [xCenter yCenter]);
    Screen('Flip',window);
    disp(['Start ' num2str(trial)])
    WaitSecs(1);
    
    %generate memory task
    [multiColors, patternL, patternR, LR] = genTask(numSquares);
    
    % Draw the rect to remember to the screen
    Screen('FillRect', window, [multiColors],...
        [allRectsCenter]);
     Screen('FrameRect',window, [1, 1, 1], [allRectsCenter],2);
    % Flip to the screen (begin stimulus)
    disp(['stimulus ' num2str(trial)])
    now=toc;
    time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
    DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);
    Screen('Flip', window);
    startTime = GetSecs; % Get starting time of current trial
    if DMSmode==2||DMSmode==3
        [retDMS] = StartAcquisition(); % Begin acquisition
        CheckWarning(retDMS);
    end
    WaitSecs(stimulus);
    % Clear Screen (begin retention)
%     Screen('FillRect',window,[0,0,0]);
    Screen('DrawLines', window, crossallCoords,...
        lineWidthPix, white, [xCenter yCenter]);
    now=toc;
    time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
    DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);
    Screen('Flip',window);
    disp(['retention ' num2str(trial)])
    WaitSecs(retention);
    
    % Draw the recall rects to the screen
    Screen('FillRect',window,[0,0,0]);
    Screen('FillRect', window, [patternR patternL],...
        [allRectsLeft allRectsRight]);
    Screen('FrameRect',window, [1, 1, 1], [allRectsLeft allRectsRight],2);
    now=toc;
    time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
    DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);
    
    % Flip to the screen
    Screen('Flip', window);
    disp(['probe ' num2str(trial)])
    answer = 1;
    tt=0;
    startSecs = GetSecs;
    keyCode = zeros(1,256);
    response = -1;
    keyIsDown = 0;
    while tt < probe
        endSecs = GetSecs;
        tt = endSecs - startSecs;
        if ~strncmpi(KbName(keyCode),'6',4) && ~strncmpi(KbName(keyCode),'4',4)
            [ keyIsDown, timeSecs, keyCode ] = KbCheck;
        end
        if answer && (strncmpi(KbName(keyCode),'6',4) || strncmpi(KbName(keyCode),'4',4))
            fprintf('"%s" typed at time %.3f seconds\n', KbName(keyCode), timeSecs - startSecs);
            Screen('DrawLines', window, crossallCoords,lineWidthPix, white, [xCenter yCenter]);
            now=toc;
            time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
            DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);
            Screen('Flip',window);
            answer = 0;
        end
        
    end
    now=toc;
    time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
    DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);
    Screen('Flip',window);
    
    trialStartTime = startTime - DMSStartTime;
    responseTime = timeSecs - startSecs;
    if ~isempty(KbName(keyCode))
        response = strncmpi(KbName(keyCode),'6',4);
    end
    results.TrialStartTime(trial) = trialStartTime;
    results.ResponseTime(trial) = responseTime;
    results.LR{trial} = LR;
    results.Response{trial} = response;
    results.Correct(trial) = LR==response;
    [~, results.Exposure, results.Accumulate, results.Kinetic] = GetAcquisitionTimings();
    
%     Screen('FillRect',window,[0,0,0]);
    Screen('DrawLines', window, crossallCoords,...
        lineWidthPix, white, [xCenter yCenter]);
    now=toc;
    time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
    DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);
    Screen('Flip',window);
    disp(['ITI ' num2str(trial)])
    WaitSecs(ITI);
    if DMSmode==2||DMSmode==3
        % Stop acquisition
        [retDMS] = AbortAcquisition();
        CheckWarning(retDMS);
        % Retrieve Images
        [retDMS, frameCount] = GetTotalNumberImagesAcquired();
        CheckWarning(retDMS);
        [I] = GetImagesFromCamera(XPixels,HBin,trackNum,frameCount);
        % Save data for this trial
        save([filename '\DMStrial' num2str(trial) '.mat'],'I');
    end
end

if DMSmode == 2 && recoveryTime ~= 0 %recovery
    Screen('FillRect',window,[0,0,0]);
    Screen('DrawLines', window, crossallCoords,...
        lineWidthPix, white, [xCenter yCenter]);
    disp('Starting recovery')
    Screen('Flip',window);
    [retDMS] = StartAcquisition(); % Begin acquisition
    CheckWarning(retDMS);
    WaitSecs(recoveryTime);
    % Stop acquisition
    [retDMS] = AbortAcquisition();
    CheckWarning(retDMS);
    % Retrieve Images
    [retDMS, frameCount] = GetTotalNumberImagesAcquired();
    CheckWarning(retDMS);
    [I] = GetImagesFromCamera(XPixels,HBin,trackNum,frameCount);
    save([filename '\DMSrecovery.mat'],'I');
    disp('Recovery Finished')
end

if DMSmode ~= 3
    save([filename '\DMSResults.mat'],'results')
end
% End of trials, wait for a key press
Screen('FillRect',window,[0.5,0.5,0.5]);
Screen('Flip',window);
disp('Finished: press any key')
KbStrokeWait;

% Clear the screen
sca;

%%
    function [multiColors, patternL, patternR, LR] = genTask(numSquares)
        yellows = ones(1,floor(numSquares/2)+randi([-2 2]));
        reds = zeros(1,numSquares-length(yellows));
        pattern = [reds yellows];
        pattern = pattern(randperm(length(pattern)));
        for patternInd = 1:length(pattern)
            if pattern(patternInd) == 0
                multiColors(:,patternInd) = [1;0;0];
            else
                multiColors(:,patternInd) = [1;1;0];
            end
        end
        patternL = multiColors;
        patternR = multiColors;
        
        findRed = find(pattern==0);
        findYellow = find(pattern==1);
        swapRed = datasample(findRed,1);
        swapYellow = datasample(findYellow,1);
        LR = randi([0 1],1);
        
        if LR == 0 %0 is L, 1 is R
            patternL(:,swapYellow) = multiColors(:,swapRed);
            patternL(:,swapRed) = multiColors(:,swapYellow);
        else
            patternR(:,swapYellow) = multiColors(:,swapRed);
            patternR(:,swapRed) = multiColors(:,swapYellow);
        end
        
        
    end

    function [I] = GetImagesFromCamera(XPixels,HBin,nTracks,frameCount)
        I = zeros(nTracks,XPixels/HBin,frameCount);
        currentSeries = 1;
        while(currentSeries <= frameCount)
            
            [ret, imageData] = GetOldestImage(XPixels/HBin * nTracks);
            
            if ret == atmcd.DRV_SUCCESS % data returned
                I(:,:,currentSeries) = flipdim(transpose(reshape(imageData, XPixels/HBin, nTracks)),1);
                currentSeries=currentSeries+1;
            end
        end
    end

    function timerCallback(~, ~)
        now=toc;
        time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
        Screen('DrawLines', window, crossallCoords,...
            lineWidthPix, white, [xCenter yCenter]);
        DrawFormattedText(window, time_str, horizontal_pix, vertical_pixel, white);  %, 1850, 1040
        Screen('Flip', window);
    end
end

