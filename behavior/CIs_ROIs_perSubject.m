% Constants
nBlocks = 5;
mean_ROIs = [];
max_ROIs = [];
DataMat = [];
info_parts= {
    {'CamilleG'}
    {'CarolaneD'}
    {'CarolineS'}
    {'EmilieR'}
    {'GuillaumeB'}
    {'LaurenceC'}
    {'MarietheQ'}
    {'MelanieP'}
    {'SVT'}
    {'YanniM'}
    {'Ziad'}
    };
nSubs = length(info_parts);
sizex = 497;

addpath('/home/adf/faghelss/Documents/TOOLBOXES/Colormaps')
addpath('/home/adf/faghelss/Documents/TOOLBOXES/Stat4CiToolbox/DISPLAY_RESULTS')
addpath('/home/adf/faghelss/Documents/TOOLBOXES/Stat4CiToolbox')
ims=load('/home/adf/faghelss/Documents/miscelllaneous/cc_ims.mat');ff=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/maskellipse.mat');
%%

subjectN=0;
%for
 task=2
tasks_folders={'emotion','gender'};
task_str=tasks_folders{task};

nBlocks=5;
sigma=15;
bubStd=12;
nTrials=200;
spaceSize=497;

%load('simulated_LeyeUser.mat','cid')

facemask=squeeze(ff.facemask(:,:,1));

clear CI_all
for subject=1:nSubs
    subjectN=subjectN+1;
    
    sprintf('// subject %d  //',subjectN)
    name=info_parts{subject}{1};
    for block=1:nBlocks
        
        clear cid q
        
        sprintf('/ block %d /',block)
        
        
        load(strcat('Bubbles_SSVEP_',name,'_',task_str,'_',num2str(block),'.mat'))
        
        X_all=[];
        y_all=[];
        
        X           = zeros(nTrials, spaceSize^2);
        y           = zeros(1, nTrials);
        
        rng(cid.seed.Seed,cid.seed.Type)
        
        % % this is required for this specific experiment; the random number
        % generator won't be "at the right place" otherwise.
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
            

            
            % Creation de bruit
            qteBulles = cid.data(5,trial);
            [posit,stim]=bubMask_spatial(qteBulles,spaceSize,facemask,bubStd);
            tmp=posit;
            X(trial,:) =tmp(:);%mask_bulles(:);%cid.masks{trial};% tmp(:);
            
            % UPDATE QUEST
            q=QuestUpdate(q, surface, y(1,trial));surface=QuestMean(q);
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
        
    end
    
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
    
    CI_all{task}(subjectN,:,:)=SSP;
    
        
    figure, imagesc(squeeze((CI_all{task}(subjectN,:,:)))), title(sprintf(' subject %s, task %s',name,task_str))
    
end
%% % % % %  % % % % %  % % % % %  % % % % % % % % % % % % % % 
% % % % % compute ROI usage, save per participant % % % % % 
% load ROI
mask_file=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/interestregions.mat')
tasks_folders={'emotion','gender'};
task_str=tasks_folders{task};
% Link CI and ROIs
masks{1}  = logical(mask_file.mask(:,:,1));
masks{2} = logical(mask_file.mask(:,:,2));
masks{3} = logical(mask_file.mask(:,:,3));
masks{4} = ~(masks{1}+masks{2}+masks{3});

if task==1
    path_data='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/emotion/raw_data/';
    res_path='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/emotion/emotion_CIs';
else
    path_data='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/gender/raw_data/';
    res_path='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/gender/gender_CIs';
end
subjectN=0;
clear all_meanROI_usage all_maxROI_usage
for subject=1:nSubs
    subjectN=subjectN+1;
    
    name=info_parts{subject}{1};
    filename = sprintf('CI_%s_%s', name,task_str);
    
    CIs_sub=squeeze((CI_all{task}(subject,:,:)));
    clear cid q
    % Calculate de mean and max
    mean_ROIs = [nanmean(CIs_sub(masks{1})), nanmean(CIs_sub(masks{2})), nanmean(CIs_sub(masks{3}))];
    max_ROIs = [max(CIs_sub(masks{1})), max(CIs_sub(masks{2})), max(CIs_sub(masks{3}))];
    
    
    figure, imagesc(squeeze((CI_all{task}(subjectN,:,:)))), title(sprintf(' subject %s, task %s',name,task_str))
    
    save(fullfile(res_path,filename), 'CIs_sub', 'mean_ROIs', 'max_ROIs','masks')%,'mACC','mRT');
end

show_ci=squeeze(sum((CI_all{task})))/sqrt(subjectN);
figure, imagesc(show_ci),colormap(jet)
%%
task=2;


ims=load('/home/adf/faghelss/Documents/miscelllaneous/cc_ims.mat');ff=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/maskellipse.mat');
tasks_folders={'emotion','gender'};
task_str=tasks_folders{task};
facemask=squeeze(double(ff.facemask(:,:,1)));
if task==1
    path_data='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/emotion/raw_data/';
    res_path='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/emotion/emotion_CIs';
    face_img_show=stretch(double(squeeze(ims.cc_ims{11}))).*facemask;
else
    path_data='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/gender/raw_data/';
    res_path='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/gender/gender_CIs';
    face_img_show=stretch(double(squeeze(ims.cc_ims{11}))).*facemask;
end


p = .00001;
tC = 3.0;
counter=0;
sigma=15;
counter=counter+1;
FWHM(counter) = sigma * sqrt(8*log(2));

[volumes, N] = CiVol(sum(facemask(:)), 2);
[tP, k,tC1] = stat_threshold(volumes, N, FWHM, Inf, p, tC)

subjectN=0;

for subject=1:nSubs
    subjectN=subjectN+1;
    
    name=info_parts{subject}{1};
    filename = sprintf('CI_%s_%s', name,task_str);
    ind=load (fullfile(res_path,filename), 'CIs_sub', 'mean_ROIs', 'max_ROIs','masks')
    CI_all{task}(subjectN,:,:)=ind.CIs_sub;
    % Calculate de mean and max
    mean_ROIs = [nanmean(ind.CIs_sub(masks{1})), nanmean(ind.CIs_sub(masks{2})), nanmean(ind.CIs_sub(masks{3}))];
    max_ROIs = [max(ind.CIs_sub(masks{1})), max(ind.CIs_sub(masks{2})), max(ind.CIs_sub(masks{3}))];
    
    all_meanROI_usage{task}(subjectN,:)=(mean_ROIs);all_maxROI_usage{task}(subjectN,:)=(max_ROIs);
    
end
figure, scatter3(all_meanROI_usage{task}(:,1),all_meanROI_usage{task}(:,2),all_meanROI_usage{task}(:,3),'filled'),xlabel('left eye usage'),...
    ylabel('right eye usage'),zlabel('mouth usage')
hold on
 scatter3(all_meanROI_usage{1}(:,1),all_meanROI_usage{1}(:,2),all_meanROI_usage{1}(:,3),'filled'),legend({'gender task','emotion task'})
%
show_ci=squeeze(sum(CI_all{task}))/sqrt(subjectN);

%,4.3);

[signif_regions]=find_sig_cluster(show_ci,k,tC);
tmp_range=[-3 20];%[min(show_ci(:)) max(show_ci(:))];

show_CI=showDatCi(signif_regions,1,show_ci,jet,face_img_show,tmp_range);
[complete, thresh]=overlay_pixel(face_img_show,show_ci,tP);

figure, imshow(thresh)
imwrite(thresh,strcat('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/figures/ClassificationImages/','averageCI_',task_str,'.png'))
