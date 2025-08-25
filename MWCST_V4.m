function [result]= MWCST(XPixels,HBin,trackNum,filename,CSmode)
%%%%%%%  Pan modified
%setting baseline length
baselineTime = 120; %sec %background:5, trail:120
numTrials = 48; %trials %background:3, trail:48

% make output folder
if CSmode == 2
    folderNumber=1;
    while 1
        CardSortingfolder = [filename '\CS_' num2str(folderNumber)];
        if ~exist(CardSortingfolder)
            mkdir(CardSortingfolder)
            break
        end
        folderNumber=folderNumber+1;
    end
    filename = CardSortingfolder;
end
%% start exp
try
  
    % Setup PTB with some default values
    Screen('Preference', 'SkipSyncTests', 1);%*****have to fix, 1
    PsychDefaultSetup(2);
    
    % Seed the random number generator. Here we use the an older way to be
    % compatible with older systems. Newer syntax would be rng('shuffle'). Look
    % at the help function of rand "help rand" for more information
    rand('seed', sum(100 * clock));
    
    % Set the screen number to the external secondary monitor if there is one
    % connected
    screenNumber = max(Screen('Screens'));
    
    % Define black, white and grey
    white = WhiteIndex(screenNumber);
    grey = white / 2;
    black = BlackIndex(screenNumber);
    
    % Open the screen
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], 32, 2);
    
    % Flip to clear
    Screen('Flip', window);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Set the text size
    Screen('TextSize', window, 60);
    
    % Query the maximum priority level
    topPriorityLevel = MaxPriority(window);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
      %Obtain some information of experiment
    %     [info]=infoFun;
    %     HideCursor
    temp1=imread('image/111.png');temp2=imread('image/222.png');temp3=imread('image/333.png');temp4=imread('image/444.png');
    para=struct(...
        'bt1',temp1,...
        'bt2',temp2,...
        'bt3',temp3,...
        'bt4',temp4...
        );
    load('image/cardArray','cardArray');
    background=zeros(screenYpixels,screenXpixels,3);
%      background=zeros(1080,1920,3);
%     background(241:440,621:745,:)=para.bt1;
    background(round(screenYpixels*1/5):round(screenYpixels*1/5)+199,round(screenXpixels*1/3):round(screenXpixels*1/3)+124,:)=para.bt1;
    background(round(screenYpixels*1/5):round(screenYpixels*1/5)+199,round(screenXpixels*1/3)+124+62:round(screenXpixels*1/3)+2*124+62,:)=para.bt2;
    background(round(screenYpixels*1/5):round(screenYpixels*1/5)+199,round(screenXpixels*1/3)+2*124+2*62:round(screenXpixels*1/3)+3*124+2*62,:)=para.bt3;
    background(round(screenYpixels*1/5):round(screenYpixels*1/5)+199,round(screenXpixels*1/3)+3*124+3*62:round(screenXpixels*1/3)+4*124+3*62,:)=para.bt4;
    background=uint8(background);
    % Set the blend funciton for the screen
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    
    %----------------------------------------------------------------------
    %                       Timing Information
    %----------------------------------------------------------------------
    
    % Interstimulus interval time in seconds and frames
    isiTimeSecs = 1;
    isiTimeFrames = round(isiTimeSecs / ifi);
    
    % Numer of frames to wait before re-drawing
    waitframes = 1;
    
    
    %----------------------------------------------------------------------
    %                       Keyboard information
    %----------------------------------------------------------------------
    
    % Define the keyboard keys that are listened for. We will be using the left
    % and right arrow keys as response keys for the task and the escape key as
    % a exit/reset key
    KbName('UnifyKeyNames')
    escapeKey = KbName('ESCAPE');
    z=KbName('4');
    x=KbName('5');
    c=KbName('6');
    v=KbName('+');
    %----------------------------------------------------------------------
    %                     Card and data
    %----------------------------------------------------------------------
    
    % We are going to use four colors for this test.
    colorList = {'Red','Green','Blue','Yellow'};
    shapeList = {'Circle','Triangle','Cross','Star'};
    numberList = {'1','2','3','4'};
    % Number of trials per condition. We set this to one for this demo, to give
    % us a total of 9 trials.
    continueCondition = 6;
    
    % Duplicate the condition matrix to get the full number of trials
    
    
    % Get the size of the matrix each of rule have 6 cards and two test all
    % have 48 cards
    
    %
    % Randomise the conditions
    %     shuffler = Shuffle(1:numTrials);
    %     condMatrixShuffled = condMatrix(:, shuffler);
    % construct T & F image
    F=imread('image/sadred.png');
    T=imread('image/happygreen.png');
    RIMG=background;
    RIMG(round(screenYpixels*4/7):round(screenYpixels*4/7)+251,round(screenXpixels*3/5):round(screenXpixels*3/5)+274,:)=F;
    FTexture = Screen('MakeTexture', window, RIMG);
    RIMG=background;
    RIMG(round(screenYpixels*4/7):round(screenYpixels*4/7)+251,round(screenXpixels*3/5):round(screenXpixels*3/5)+274,:)=T;
    TTexture = Screen('MakeTexture', window, RIMG);
    %----------------------------------------------------------------------
    %                     Make a response matrix
    %----------------------------------------------------------------------
    % This is a four row matrix the first row will record the word we present,
    % the second row the color the word it written in, the third row the key
    % they respond with and the final row the time they took to make there response.
    
    cardOrder = randperm(numTrials);   
    continue_right=0;
    
    %%%%%%%pan modified
    % Here we set the size of the arms of our fixation cross
    fixCrossDimPix = 40;

    % Now we set the coordinates (these are all relative to zero we will let
    % the drawing routine center the cross in the center of our monitor for us)
    crossxCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    crossyCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    crossallCoords = [crossxCoords; crossyCoords];
    lineWidthPix = 4;
    
    if CSmode == 2 && baselineTime ~= 0
        Screen('DrawLines', window, crossallCoords,...
            lineWidthPix, white, [xCenter yCenter]);
        disp('Baseline: displaying cross')
        Screen('Flip',window);
        [retCST] = StartAcquisition(); % Begin acquisition
        CheckWarning(retCST);
        
%         % add time stamps (modified by TY)
%         tic;
%         Screen('TextSize', window, 40);
%         t=timer('ExecutionMode','fixedRate','Period',1,'TimerFcn',@timerCallback);
%         start(t);
%         pause(baselineTime);
        WaitSecs(baselineTime);
%         stop(t);
%         delete(t);
%         Screen('TextSize', window, 60);
        
        % Stop acquisition
        [retCST] = AbortAcquisition();
        CheckWarning(retCST);
        % Retrieve Images
        [retCST, frameCount] = GetTotalNumberImagesAcquired();
        CheckWarning(retCST);
        [I] = GetImagesFromCamera(XPixels,HBin,trackNum,frameCount);
        save([filename '\CSTbaseline.mat'],'I');
        disp('Baseline data is saved!!!')
        Screen('FillRect',window,[0.5,0.5,0.5]);
        Screen('Flip',window);
        disp('Finished: press anykey')
        KbStrokeWait;
        Screen('FillRect',window,[0,0,0]);
    end
    
  %%%%%%%
  
    %----------------------------------------------------------------------
    %                       Experimental loop
    %----------------------------------------------------------------------
    imagelist=[];
    temp=cell(3,1);
    
    % Animation loop: we loop for the total number of trials
    %the answer is decide after press button and if the first three rule is
    %1 2 3 after rule is 1 2 3
    rule_decide=0;
    for trial = 1:numTrials
        % color and shape and number num randnum 1 to 4
        
        % The color and the shape and the number it is drawn in
        
        
        % Cue to determine whether a response has been made
        respToBeMade = true;
        
        % If this is the first trial we present a start screen and wait for a
        % key-press
        if trial == 1 && CSmode == 2
            DrawFormattedText(window, 'Press Any Key To Begin',...
                'center', 'center', white);
            Screen('TextSize', window, 60);
            Screen('Flip', window);
            KbStrokeWait;%any key to begin
            [retCS] = StartAcquisition(); % Begin acquisition
            CheckWarning(retCS);
            Screen('TextSize', window, 40);
            tic;
        end
        
        
        % Flip again to sync us to the vertical retrace at the same time as
        % drawing our fixation point
        Screen('DrawDots', window, [xCenter; yCenter], 10, white, [], 2);
        vbl = Screen('Flip', window);
        
        % Now we present the isi interval with fixation point minus one frame
        % because we presented the fixation point once already when getting a
        % time stamp
        for frame = 1:isiTimeFrames - 1
            
            % Draw the fixation point
            Screen('DrawDots', window, [xCenter; yCenter], 10, white, [], 2);
            
            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end
        
        % Now present the word in continuous loops until the person presses a
        % key to respond. We take a time stamp before and after to calculate
        % our reaction time. We could do this directly with the vbl time stamps,
        % but for the purposes of this introductory demo we will use GetSecs.
        %
        % The person should be asked to respond to either the written word or
        % the color the word is written in. They make thier response with the
        % three arrow key. They should press "Left" for "Red", "Down" for
        % "Green" and "Right" for "Blue".
        tStart = GetSecs;
        while respToBeMade == true
            theImage=imread(fullfile('image','question',[num2str(cardArray(1,cardOrder(trial)))...
                num2str(cardArray(2,cardOrder(trial))) num2str(cardArray(3,cardOrder(trial))) '.png']));
            cardAttributs=cardArray(:,cardOrder(trial));
            stiIMG=background;
            stiIMG(round(screenYpixels*5/7):round(screenYpixels*5/7)+199,round(screenXpixels*1/2)-62:round(screenXpixels*1/2)+62,:)=theImage;
            % Draw the word
            %             DrawFormattedText(window, char(theWord), 'center', 'center', theColor);
            imageTexture = Screen('MakeTexture', window, stiIMG);
            Screen('DrawTexture', window, imageTexture, [], [], 0);
            
            % add time stamps
            if CSmode==2
                now=toc;
                time_str=[num2str(floor(now/60)) ':' num2str(floor(rem(now,60)))];
                DrawFormattedText(window, time_str, 1600, 1040, white);  %, 1850, 1040
            end
                
            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            Screen('Close', imageTexture);   %20250822 add by CJ

            % Check the keyboard. The person should press the
            
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                %                 ShowCursor;
                [retCS] = AbortAcquisition();
                CheckWarning(retCS);
                sca;
                return
                
            elseif keyCode(z)
                response = 1;
                respToBeMade = false;
                %descide
                if rule_decide==0
                    [rule_decide,ruleOrder]=rule(response,cardAttributs);
                    if rule_decide==1
                        ruleOrder=[ruleOrder ruleOrder];
                        rule_num=1;
                        
                    else
                        rule_num=0;
                    end
                    
                end
                if rule_num ~= 0
                    rule_now=ruleOrder(rule_num);
                else
                    rule_now=0;
                end
                
                [R,C]=checkAns(response,rule_now,cardAttributs);
                Screen('DrawTexture', window,eval([R 'Texture']), [], [], 0);
                Screen('Flip', window);
                WaitSecs(1);
            elseif keyCode(x)
                response = 2;
                respToBeMade = false;
                %descide
                if rule_decide==0
                    [rule_decide,ruleOrder]=rule(response,cardAttributs);
                    if rule_decide==1
                        ruleOrder=[ruleOrder ruleOrder];
                        rule_num=1;
                        
                    else
                        rule_num=0;
                    end
                    
                end
                if rule_num ~= 0
                    rule_now=ruleOrder(rule_num);
                else
                    rule_now=0;
                end
                
                [R,C]=checkAns(response,rule_now,cardAttributs);
                Screen('DrawTexture', window,eval([R 'Texture']), [], [], 0);
                Screen('Flip', window);
                WaitSecs(1);
            elseif keyCode(c)
                response = 3;
                respToBeMade = false;
                %descide
                if rule_decide==0
                    [rule_decide,ruleOrder]=rule(response,cardAttributs);
                    if rule_decide==1
                        ruleOrder=[ruleOrder ruleOrder];
                        rule_num=1;
                        
                    else
                        rule_num=0;
                    end
                    
                end
                if rule_num ~= 0
                    rule_now=ruleOrder(rule_num);
                else
                    rule_now=0;
                end
                
                [R,C]=checkAns(response,rule_now,cardAttributs);
                Screen('DrawTexture', window,eval([R 'Texture']), [], [], 0);
                Screen('Flip', window);
                WaitSecs(1);
            elseif keyCode(v)
                response = 4;
                respToBeMade = false;
                %descide
                if rule_decide==0
                    [rule_decide,ruleOrder]=rule(response,cardAttributs);
                    if rule_decide==1
                        ruleOrder=[ruleOrder ruleOrder];
                        rule_num=1;
                        
                    else
                        rule_num=0;
                    end
                    
                end
                if rule_num ~= 0
                   
                    rule_now=ruleOrder(rule_num);
                else
                    rule_now=0;
                end
                
                [R,C]=checkAns(response,rule_now,cardAttributs);
                Screen('DrawTexture', window,eval([R 'Texture']), [], [], 0);
                Screen('Flip', window);
                WaitSecs(1);
            end
            
        end
        tEnd = GetSecs;
        rt = tEnd - tStart;
        if R=='T'
            continue_right=continue_right+1;
            
        elseif R=='F'
            continue_right=0;
        end
        %disp(['Number of correct response' num2str(continue_right)]);
        % Record the trial data into out data matrix
        data{trial,1} = rule_now ;
        data{trial,2} = cardAttributs;
        data{trial,3} = response;
        data{trial,4} = rt;
        data{trial,5} = C;
        if continue_right==continueCondition
            rule_num=rule_num+1;
            continue_right=0;
            disp('Finish one rule and change to another rule!')
        end
        if rule_num==7
            break
        end

        fid = fopen('ResponseRecord.txt', 'w');
        fprintf(['Number of correct response' num2str(continue_right) '\n']);
        %fprintf('Finish one rule and change to another rule!')
        fclose(fid);

    end
    if CSmode == 2
        % Stop acquisition
        [retCS] = AbortAcquisition();
        CheckWarning(retCS);
        Screen('TextSize', window, 60);
        % Retrieve Images
        [retCS, frameCount] = GetTotalNumberImagesAcquired();
        CheckWarning(retCS);
        [I] = GetImagesFromCamera(XPixels,HBin,trackNum,frameCount);
        save([filename '\CardSorting.mat'],'I');
    end

    header = {'Rule','Ans', 'Resp', 'RT', 'Correct'};
    data_table = cell2table(data, 'VariableNames', header);
    % Create a csv file to save data
    [~, results.Exposure, results.Accumulate, results.Kinetic] = GetAcquisitionTimings();
    save([filename '\results.mat'],'results');
    exp_data = strcat([filename '\'], 'exp_', date, '.csv');
    writetable(data_table, exp_data);
    result.datatable=data_table;
    result.Imagelist=imagelist;
    
    % End of experiment screen. We clear the screen once they have made their
    % responsev
    DrawFormattedText(window, 'Experiment Finished \n\n Press Any Key To Exit',...
    'center', 'center', white);
    Screen('Flip', window);
    KbStrokeWait;
    sca;
    
    disp('Succeed!');
catch ME
    %     ShowCursor
    
    rethrow(ME);
    disp('Error!')
end
% function [info]=infoFun
% prompt = {'Subject Number','Gender[1 = m, 2 = f]','Age'};
% title = 'Exp infor'; % The title of the dialog box
% definput = {'','','',''}; % Default input value(s)
% % Using inputdlg() to obtain information and save it to cell array
% info=inputdlg(prompt,title,[1, 50],definput);
% end
    function [rule_decide, answer_order]=rule(response,cardAttributes)
        if response==cardAttributes(1) % colorrule
            ruleNum=1;
            rule_decide=1;
        elseif response==cardAttributes(2)% shaperule
            ruleNum=2;
            rule_decide=1;
        elseif response==cardAttributes(3) % numrule
            ruleNum=3;
            rule_decide=1;
        else
            ruleNum=0;
            rule_decide=0;
        end
        switch ruleNum
            case 0
                answer_order=0;
            case 1
                order_rand=rand(1);
                if order_rand<0.5
                    answer_order=[1 2 3];
                else
                    answer_order=[1 3 2];
                end
            case 2
                order_rand=rand(1);
                if order_rand<0.5
                    answer_order=[2 1 3];
                else
                    answer_order=[2 3 1];
                end
            case 3
                order_rand=rand(1);
                if order_rand<0.5
                    answer_order=[3 2 1];
                else
                    answer_order=[3 1 2];
                end
        end
    end
    function [R,C]=checkAns(response,ruleAns,Attributes)
        if ruleAns~=0
            if response == Attributes(ruleAns,1)
                R='T';
                C=1;
            else
                R='F';
                C=0;
            end
        elseif ruleAns==0
            R='F';
            C=0;
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
        DrawFormattedText(window, time_str, 1600, 1040, white);  %, 1850, 1040
        Screen('Flip', window);
    end
end