function spatial_sampling
%% Bubbles spatial sampling in face gender and emotion recognition tasks
%  write descrption here


%   Author:     Simon Faghel-Soubeyrand
%   Date:       feb/2018
%   Version:    1 
%   Tested by:  _____________ 

% TO ADJUST DISTANCE
% Distance must be ajusted so that the width of the image is 6 deg of
% visual angle. e.g..
%
% stim_deg=6
% screen_width_cm = 30;                                   % screen width in cm on Simon's  laptop, is 60 cm for Gosselin's lab                                              
% pixels_cm = screen_width_cm/rect(3);                    % size of 1 pixel in cm
% stim_pix = 497;                                         % params.ImSize in pixels
% stimuli_width = stim_pix*pixels_cm;                     % to measure stim on screen in cm, for instance
% viewing_distance=  (stimuli_width/2)/tand(stim_deg/2)  % gives viewing distance in cm
% viewing_distance_cm = 76;                              % viewing distance in cm

close all; clear;

%% folder organisation parameters
startPath=pwd;
functionPath = fullfile(startPath,'functionsBubbles'); addpath(genpath(functionPath));
StimuliPath  = fullfile(startPath,'stimuli');
resultsPath  = fullfile(startPath,'results');
instrucPath  = fullfile(startPath,'instructions');


%% query experiment's details

% Participant identifier
prompt       = {'Enter the subject identifier:'};
DatString    = 'Participant setup';
numlines     = 1;
answer       = inputdlg(prompt,DatString,numlines);
params.sID   = answer{1};
params.name = params.sID;

% Enter the task
prompt       = {'Enter the experiment (type gender or emotion)'};
DatString    = 'Experiment setup';
numlines     = 1;
answer       = inputdlg(prompt,DatString,numlines);
params.condition  = answer{1};

% query the experiment details
prompt       = {'Enter the block number'};
DatString    = 'Experiment setup';
numlines     = 1;
answer       = inputdlg(prompt,DatString,numlines);
params.blockNb      = str2double(answer{1});


% define parameters
prompt       = {'Enter the Display''s CFG name (e.g. labo EEG)'};
DatString    = 'Display setup';
numlines     = 1;
answer       = inputdlg(prompt,DatString,numlines);
params.DispCFG  = answer{1};


switch params.condition
    case 'emotion' % emotion task ( happy versus fear)
        params.possibleKeys   = 'pq';
                
        neutralKey='p';
        happyKey ='q';
        
        file_instru='TestInstru_emotionSSVEP';
        
    case 'gender' % gender recognition task
        params.possibleKeys   = 'fh';
        file_instru='TestInstru_genderSSVEP';
        womenKey='f';
        menKey ='h';
        
    otherwise
        % Enter the task
        prompt       = {'woops! The condition argument must either be "gender" or "emotion"'};
        DatString    = 'Experiment setup';
        numlines     = 1;
        answer       = inputdlg(prompt,DatString,numlines);
        params.condition  = answer{1};
end


stimsize_deg    = 6; % degrees of visual angle
nRepetitions    =  5; % Number of blocks

switch params.DispCFG
    
    case 'laboEEG'
        
        %------------------------------------------------------------------------
        % DISPLAY : Physical parameters
        %-----------------------------------------------------------------------
        screen_width_cm = 60; %was 60                                                                         % screen width in cm
        viewing_distance_cm = 76;                                                                       % viewing distance in cm
        %                 pixels_cm = screen_width_cm/rect(3);                                                           % size of 1 pixel in cm`
        
    case 'simons'
        %------------------------------------------------------------------------
        % DISPLAY : Physical parameters
        %-----------------------------------------------------------------------
        screen_width_cm = 30; %was 60                                                                         % screen width in cm
        viewing_distance_cm = 76;                                                                       % viewing distance in cm
        %                 pixels_cm = screen_width_cm/rect(3);                                                           % size of 1 pixel in cm`
    otherwise
        
        % define parameters
        prompt       = {'woops!Please enter the Display''s CFG name again (e.g. labo EEG)'};
        DatString    = 'Display setup';
        numlines     = 1;
        answer       = inputdlg(prompt,DatString,numlines);
        params.DispCFG  = answer{1};
end


CFG    =   params.DispCFG;


% update params struct
params.cfg            = CFG;
params.nruns          = nRepetitions;
params.StimuliPath    = StimuliPath;
params.resultsPath    = resultsPath;
params.ImSize         = 497;
params.displ          = zeros(params.ImSize);
params.nTrials        = 100;
params.QuitKeys       = KbName('escape');
params.instrupath     = instrucPath;
params.fixrange       = (params.ImSize/2)-3:(params.ImSize/2)+3;



StimPath = fullfile(params.StimuliPath,'cc_ims.mat');
load (StimPath,'cc_ims');

im = cc_ims;
load (strcat(params.StimuliPath,'/maskellipse.mat'));
facemask=logical(squeeze(double(facemask(:,:,1))));

% cid-file-already-exists test
file_name = sprintf('Bubbles_SSVEP_%s_%s_%d.mat', params.name, params.condition,params.blockNb);
fullPath = fullfile(params.resultsPath,file_name);
fid = fopen(fullPath);
if fid>0
    fclose(fid);
    error('File already exists.')
end

KbName('UnifyKeyNames');
if params.blockNb==1
% number of bubbles
prompt       = {'Enter the number of bubbles to show at the begining of the block'};
DatString    = 'Experiment setup';
numlines     = 1;
answer       = inputdlg(prompt,DatString,numlines);
params.qteBulles  = str2double(answer{1});
else
previousfile_name = sprintf('Bubbles_SSVEP_%s_%s_%d.mat', params.name, params.condition,params.blockNb-1);
fullPath = fullfile(params.resultsPath,previousfile_name); 
lastBlck=load(fullPath);
params.qteBulles=ceil(mean(lastBlck.cid.data(5,(end-20):end)));
params.blockNb
end


%%

% -------------------------------------------------------------------------------------------------------------------
%
% TEST VALIDITY OF INPUT ARGUMENTS AND MAKE APPROPRIATE INITS
% 
% -------------------------------------------------------------------------------------------------------------------


% -------------------------------------------------------------------------------------------------------------------
%
% SET THE EXPERIMENTAL PARAMETERS AND INITIALIZE OTHER VARIABLES
%
% -------------------------------------------------------------------------------------------------------------------
params.info = 'comparing individual visual strategies (Face Gender or EXNEX) with Bubbles and SSVEPs ';

try

seed_0 = round(sum(100*clock));
cid.noise = sprintf('initial seed value = %d', seed_0);
cid.participant=params.name;
cid.dataLabels = 'gender_face1 face1_nb face2_nb right_left_flip nb_bubbles overlap_target response RT accuracy simi2target';
nbColumns = 11;

seed_1=rng('shuffle')%  rand('state', round(sum(100*clock)));

cid.seed=rng;


nTrials=params.nTrials;
nBubbles=params.qteBulles;

nObjects = 264; %%%%% % nb of different stimuli
nIDs = 66;
stimMat = zeros(nObjects,5);
stimMat(:,1) = 1:nObjects; % object
aa = repmat(1:nIDs, [2 1]);
stimMat(:,2) = repmat(aa(:), [2 1]); % ID
stimMat(:,3) = [2*ones(nIDs,1); ones(nIDs,1); 2*ones(nIDs,1); ones(nIDs,1)]; % gender 1=H 2=F
stimMat(:,4) = repmat([2; 1], [nIDs*2 1]); % expr 1=N 2=H
stimMat(:,5) = (2*stimMat(:,3) + stimMat(:,4))-2; % types 1=HN 2=HH 3=FN 4=FH

cid.data       = zeros(nbColumns, nTrials);
idx_stims=randperm(264);
cid.data(1, :) = stimMat(idx_stims(1:nTrials),3);% gender of face 1: 1 = male, 2 = female
cid.data(2, :) = stimMat(idx_stims(1:nTrials),1);%
cid.data(3, :) = stimMat(idx_stims(1:nTrials),4);%
cid.data(4, :) = rand(1,nTrials)>.5; % 0 = no flip, 1 = right-left flip
cid.data(5, 1) = nBubbles;
accuracies = zeros(1,nTrials);

% ----------------------------------------------------------------------------------------------------------
% INITIALIZE QUEST
% -----------------------------------------------------------------------------------------------------------
bubbleStDev=12;
[tGuess,minBulle,surfaceBulle]= bubbles_questGuest(params.qteBulles,params.ImSize,bubbleStDev);

tGuessSd=0.05; % sd of Gaussian before clipping to specified range
pThreshold=0.75;
beta=3;delta=0.01;gamma=0.5;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma);
surface=QuestMean(q);


% Display controls
bkg=127;

GetChar;
HideCursor;
% ----------------------------------------------------------------------------------------------------------
% INSTRUCTIONS
% ----------------------------------------------------------------------------------------------------------
fileInstru=fullfile(params.instrupath,file_instru);
XpInstruction(fileInstru,'m',params.name)                         % Get .m text instruction and display on different screen before begining task


% ----------------------------------------------------------------------------------------------------------
% INITIALIZE ON- AND OFF-SCREENS
% -----------------------------------------------------------------------------------------------------------
Screen('Preference', 'SkipSyncTests', 0);
AssertOpenGL;
screens=Screen('Screens');
screenNumber=max(screens);
[window,rect]=Screen('OpenWindow',screenNumber,bkg);


Screen(window,'FillRect', bkg);
Screen('Flip', window);

% ----------------------------------------------------------------------------------------------------------
%
% EXPERIMENTAL LOOP
%
% -----------------------------------------------------------------------------------------------------------

% ellipse=SmoothCi(imresize(ellipse,[params.ImSize  params.ImSize]),2);
% ellipse=fabriquer_ellipse(params.ImSize,.6*pi,.8*pi);
masque=facemask;%ellipse;% max(sum((ellipse==1))) for width of face in pixels.
quitXP=0;
KbWait;
for trial=1:nTrials
    
    % ----------------------------------------------------------------------------------------------------------
    % MAKE STIMULUS
    % -----------------------------------------------------------------------------------------------------------
    
    % Generate the bubbles mask
     qteBulles = cid.data(5,trial);
     mask=bubMask_spatial(qteBulles,params.ImSize,masque,bubbleStDev);
     
    % load images
    if cid.data(1,trial)==1   % male face
        imDisp=double(im{cid.data(2,trial)});
    elseif cid.data(1,trial)==2  % female face
        imDisp= double(im{cid.data(2,trial)});
    end
    
    % right-left flip images1
    if cid.data(4,trial)==1
        for ii = 1:3
            imDisp(:,:,ii) = fliplr(imDisp(:,:,ii));
        end
    end
    
    % reveal partially with the bubbles
    imDisp=uint8((imDisp-bkg).*mask+bkg);
   
    
    texturePtr = Screen('MakeTexture', window, imDisp);
    
    % PUT THE IMAGE ON THE SCREEN
    % 	- Fixation cross (.75sec)
    Screen(window,'DrawLine', [0 0 0], rect(3)/2,rect(4)/2-10,rect(3)/2,rect(4)/2+10);
    Screen(window,'DrawLine', [0 0 0], rect(3)/2-10,rect(4)/2,rect(3)/2+10,rect(4)/2);
    Screen('Flip', window);
    WaitSecs(.75);
    % 	- Blank (.25sec)
    Screen(window,'FillRect', bkg)
    Screen('Flip', window);
    WaitSecs(.25);
    % 	- The image, until response
    %    Screen('DrawTexture', window, texturePtr);
    destRect = [round(rect(3)/2)-(params.ImSize/2) round(rect(4)/2)-(params.ImSize/2) round(rect(3)/2)+(params.ImSize/2) round(rect(4)/2)+(params.ImSize/2)]; % makes image 256 x 256
    Screen('DrawTexture', window, texturePtr, [], destRect);
    Screen('Flip', window);
    
    % GET RESPONSE FROM KEYBOARD
    FlushEvents('keyDown');
    secs1 = GetSecs;
    while 1
        [keyIsDown, secs2, keyCode, ~] = KbCheck;
        if keyIsDown
            temp = KbName(keyCode);
            whichKey = findstr(temp, params.possibleKeys);
            if ~isempty(whichKey)
                break;
            end
            
            if findstr(temp, params.QuitKeys)
                quitXP=1;
                break;
            end
            
        end
    end
    
    
 
    % RECORD RSPONSE IN CID STRUCT
    cid.data(7,trial) = KbName(keyCode); % RESPONSE
    cid.data(8,trial) = secs2-secs1; % RT
    
    % ACCURATE OR NOT?
    accuracy = 0;
    switch params.condition
        case 'gender'
            if  ((cid.data(7,trial)==menKey) && (cid.data(1,trial)==1)) || ((cid.data(7,trial)==womenKey) && (cid.data(1,trial)==2))
                accuracy = 1;
            end
        case 'emotion'
            if  ((cid.data(7,trial)==neutralKey) && (cid.data(3,trial)==1)) || ((cid.data(7,trial)==happyKey) && (cid.data(3,trial)==2))
                accuracy = 1;
            end
    end
    cid.data(9,trial) = accuracy;
    
    
    % UPDATE QUEST
    q=QuestUpdate(q, surface, accuracy);surface=QuestMean(q);
    qteBulles=max(round(((surface+minBulle)*params.ImSize^2)/surfaceBulle),1);
    if trial<nTrials
        cid.data(5,trial+1)=qteBulles;
    end
    
    
    if quitXP==1
        break;
    end
    

end
GetChar
ShowCursor
sca

% SAVE DATA (CID STUCT) IN RESULTS
dataPath = fullfile(params.resultsPath,file_name);
save(dataPath,'cid');

% some computations
nb_bubbles = cid.data(5,:);
accuracies(1,:) = cid.data(9,:);

s_accuracies = conv2(accuracies, ones(1,10)/10, 'valid');


figure,
subplot(1,2,1),plot(nb_bubbles)
title('Number of bubbles as a function of trials')
xlabel('trials')
ylabel('Number of bubbles')
axis([1 nTrials 0 100])
subplot(1,2,2), plot(s_accuracies)
title('Accuracy as a function of trials')
xlabel('trials')
ylabel('Accuracy (proportion)')
axis([1 nTrials .5 1])

catch ME
    sca;
    ShowCursor;
end

end

% ----------------------------------------------------------------------------------------------------------
%
% FUNCTIONS
%
% % ----------------------------------------------------------------------------------------------------------
% function XpInstruction(fileName,format)
% % %  XpInstruction displays a given format file on Screen until keypress
% %
% %   Author:     Simon Faghel-Soubeyrand
% %   Date:       October/2015
% %   Version:    1
% %   Tested by:  Frederic Gosselin
% 
% Screen('Preference', 'SkipSyncTests', 0); % uncomment if necessary
% 
% % Psychophysics inits
% AssertOpenGL;
% screens=Screen('Screens');
% screenNumber=max(screens);
% [w, wRect]=Screen('OpenWindow',screenNumber, 128,[],32,2);
% 
% %  Center coordinates
% [xCenter,yCenter]   = RectCenter(wRect);
% 
% %
% % % Open window with default settings:
% % w=Screen('OpenWindow', screenNumber,128);
% % Select specific text font, style and size:
% Screen('TextFont',w, 'Helvetica');
% Screen('TextSize',w, 23);
% Screen('TextStyle', w, 1);
% 
% WatSecs(2);
% % Read some text file, if existing. .m are prefered...
% fid = fopen(sprintf('%s.%s', fileName,format),'r');
% 
% if fid==-1
%     error('Could not open text file: it appears inexistant or unavailable');
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % read every line of the file and put it in "text" string
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% text = '';
% tl = fgets(fid);
% lcount = 0;
% try
%     while 1
%         
%         if ~ischar(tl), break, end
%         if(~isempty(line)&&line~=-1)
%             text = [text tl];
%             tl = fgets(fid);
%             % tl=textscan(tl,'%s');
%             lcount = lcount + 1;
%         end
%     end
% catch e
%     e.getReport
%     ttt
%     keyboard
% end
% 
% fclose(fid);
% text = [text newline];
% 
% % Draw centered text inside frame
% [nx, ny, bbox] = DrawFormattedText(w, text, 'center', 'center', 0);
% 
% Screen('FrameRect', w, [10 10 10], bbox+[-20 -20 +20 +20],3);
% 
% Screen('DrawText', w, 'Appuyez une touche sur le clavier pour d�buter', xCenter-250, ny+50, [0, 100, 100, 255]);
% Screen('Flip',w);
% WaitSecs(.5);
% KbWait;
% % clear Screen screens
% sca
% 
% end




function [BestGuess,minBulle,surfaceBulle]= bubbles_questGuest(qteBulles,spaceSize,bubStd)

% AJUSTEMENTS QUEST
% - On calcule la valeur a ajuster (variant de 0 a 1),
%  la qte de pixel revelee par les bulles.
bulle=bubble(bubStd);
tmp=rand(spaceSize ^2,1);
[y,~]=sort(tmp);
tmp=reshape(tmp<=y(1),spaceSize,spaceSize);
tmp=filter2(bulle,tmp);
tmp = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
surfaceBulle=sum(tmp(:)); % number of pixels in a single bubble
minBulle=(surfaceBulle/spaceSize.^2);
BestGuess=((qteBulles*surfaceBulle)/spaceSize.^2)-minBulle;


end

function bulle=bubble(bubStd)
% a single bubble
maxHalfSize = 6*bubStd;
[y,x] = meshgrid(-maxHalfSize:maxHalfSize,-maxHalfSize:maxHalfSize);
bulle = exp(-(x.^2/bubStd^2)-(y.^2/bubStd^2));
clear x y tmp
end


function mask=bubMask_spatial(qteBulles,spaceSize,masque,bubStd)


% Generate the bubbles mask
bulle=bubble(bubStd);
prob_tmp = qteBulles/sum(masque(:));
tmp=rand(spaceSize^2,1) .* masque(:);
tmp=reshape(tmp>=(1-prob_tmp),spaceSize,spaceSize); % makes the criteria probabilistic
masque2D=filter2(bulle,tmp);
masque2D = min(max(masque2D, 0), 1); % this is better

mask= repmat(masque2D,[1 1 3]);
end



