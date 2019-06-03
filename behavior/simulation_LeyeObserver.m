task=1;

path_all='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/';
if task==1
task_flder='emotion/emotion_CIs/';
else
task_flder='gender/gender_CIs/';
end
fll_path=fullfile(strcat(path_all,task_flder));
data=dir(fullfile(fll_path,'*.mat'));
nSubs=length(data);
subCount=0;

clear cid
% load CIs for every participant, correlate between tasks
clear CI_all all_use_ROI
subject=1
SUB=load(fullfile(data(subject).folder,data(subject).name));
load ('maskellipse.mat');
facemask=logical(squeeze(double(facemask(:,:,1))));

cid.data(5,1)=40;
bubStd=12;

cid.seed_0=round(sum(100*clock));
cid.seed=rng(cid.seed_0,'twister');% 
%rng('shuffle')

cid.seed=rng;
spaceSize=497;

[tGuess,minBulle,surfaceBulle]= bubbles_questGuest(cid.data(5,1),spaceSize,bubStd);

tGuessSd=0.05; % sd of Gaussian before clipping to specified range
pThreshold=0.75;
beta=3;delta=0.01;gamma=0.5;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma);
surface=QuestMean(q);
nTrials=200;
%

for trial=1:nTrials
    if mod(trial,100)==0
        sprintf('// Trial %d of simulated observer //',trial)
    end
    qteBulles=cid.data(5,trial);
    [posit,bubbles_mask]=bubMask_spatial(qteBulles,spaceSize,facemask,bubStd);
    
    cid.masks{trial}(:)=find(logical(posit(:)));
    
    prop_Leye=sum(posit(SUB.masks{1}==1));%/(spaceSize*spaceSize);%+sum(bubbles_mask(SUB.masks{2}==1))/(spaceSize*spaceSize);
    % ACCURATE OR NOT?
    accuracy = 0;
    if  prop_Leye>=4;
        accuracy = 1;
    end
 
    cid.data(9,trial) = accuracy;
    
    % UPDATE QUEST
    q=QuestUpdate(q, surface, accuracy);surface=QuestMean(q);
    qteBulles=max(round(((surface+minBulle)*spaceSize^2)/surfaceBulle),1);
    if trial<nTrials
        cid.data(5,trial+1)=qteBulles;
    end
    
end
 save('simulated_LeyeUser.mat','cid') 
%%
nBlocks=5;
sigma=20;
subjectN=0;
subjectN=subjectN+1;
for block=1:nBlocks

clear cid q

addpath('/home/adf/faghelss/Documents/TOOLBOXES/Stat4CiToolbox/Stat4CiToolbox')
%load('simulated_LeyeUser.mat','cid')

load(strcat('Bubbles_SSVEP_','CamilleG_emotion_',num2str(block),'.mat'))


X_all=[];
y_all=[];

nTrials=200;

spaceSize=497;
X           = zeros(nTrials, spaceSize^2);
y           = zeros(1, nTrials);


%cid.seed=rr;
%rng(cid.seed_0,cid.seed.Type)% SHUFFLE takes the current time, so if you
% find the time, you should be good to go. check in the rng function what
% happens when rng('shuffle') without specified seed
% rng(RandStream.shuffleSeed,'twister')
rng(cid.seed.Seed,cid.seed.Type)% SHUFFLE takes the current time, so if you




randperm(264);
rand(1,200)>.5; % 0 = no flip, 1 = right-left flip


[tGuess,minBulle,surfaceBulle]= bubbles_questGuest(cid.data(5,1),spaceSize,bubStd);
% 
tGuessSd=0.05; % sd of Gaussian before clipping to specified range
pThreshold=0.75;
beta=3;delta=0.01;gamma=0.5;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma);
surface=QuestMean(q);

clear reproducted_masks
% Reproduce the noise for each trial and put it in a matrix
for trial = 1:200
    
    if mod(trial,100)==0
        sprintf('// Trial %d  //',trial)
    end
    
    % Creation de bruit
    qteBulles = cid.data(5,trial);
    [posit,stim]=bubMask_spatial(qteBulles,spaceSize,facemask,bubStd);
%     mask_plain=zeros(spaceSize);
%     mask_plain(cid.masks{trial})=1;
%     mask_bulles=mask_plain;    
%     reproducted_masks{trial}=find(logical(posit(:)));
    tmp=posit;
    X(trial,:) =tmp(:);%mask_bulles(:);%cid.masks{trial};% tmp(:);
    
    % UPDATE QUEST
    q=QuestUpdate(q, surface, accuracy);surface=QuestMean(q);
    qteBulles=max(round(((surface+minBulle)*spaceSize^2)/surfaceBulle),1);
    
    
    % Variable dependente
    y(1,trial)  = cid.data(9,trial); % ACCURACY
    
    if find(isnan(X(trial,:)))
        toto
        X(trial,:)=zeros(1,sizeX*sizeX);
    end
    
end
% Standardisation de la VD
y       = (y - mean(y)) / std(y);

X_all = [X_all; X;];
y_all = [y_all y];

%On combine les blocs pour chq sujets..
b_all = (y_all) * X_all;
b_2D_all = reshape(b_all, spaceSize,spaceSize);
% Bootstrapping the accuracry vector.

index = randperm(size(y_all,2));
b_boot = (y_all(index)) * X_all;
b_boot_2D_all = reshape(b_boot, spaceSize,spaceSize);

mCi=mean(b_2D_all(:));
ci=SmoothCi(b_2D_all-mCi,sigma)+mCi;

mCi_boot=mean(b_boot_2D_all(:));
ci_boot=SmoothCi(b_boot_2D_all-mCi_boot,sigma)+mCi_boot;

% Standardize sci with sci_boot
SSP = (ci - mean(ci_boot(:))) / std(ci_boot(:));

CI_all(subjectN,block,:,:)=SSP;

%figure, imagesc(SSP)
    
end

figure, imagesc(squeeze(sum(CI_all)/sqrt(nBlocks)))



