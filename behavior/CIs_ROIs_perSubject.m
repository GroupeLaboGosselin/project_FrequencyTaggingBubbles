% Constants
nBlock = 5;
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

task=2;
tasks_folders={'emotion','gender'}
task_str=tasks_folders{task};

% load ROI
load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/interestregions.mat')

% Load CI, example % Bubbles_SSVEP_CarolaneD_emotion_4.mat
%path_data = fullfile('/Users/labogosselin/Desktop/RESULTS_BUBBLES_SSVEP/bubbles_31January2019/');
if task==1
path_data='/home/adf/faghelss/Documents/project_FrequencyTaggingBubblesP/behavior/emotion/raw_data/';
res_path='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/emotion/emotion_CIs';
else
path_data='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/gender/raw_data/';
res_path='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/gender/gender_CIs';
end


for subject = 1:nSubs
    %initialize 
    CIs_sub = zeros(sizex);
    name=info_parts{subject}{1};    
    filename = sprintf('CI_%s_%s', name,task_str);
    blcks=1:5;
 
    CIs_sub = quickCI_BubblesSSVEP(path_data,name, task_str,blcks); 
        
    %CIs_sub=CIs_sub/sqrt(nBlock);
    
    % Link CI and ROIs
    masks{1}  = logical(mask(:,:,1));
    masks{2} = logical(mask(:,:,2));
    masks{3} = logical(mask(:,:,3));
    masks{4} = ~(masks{1}+masks{2}+masks{3});

    % Calculate de mean and max
    mean_ROIs = [nanmean(CIs_sub(masks{1})), nanmean(CIs_sub(masks{2})), nanmean(CIs_sub(masks{3}))];
    max_ROIs = [max(CIs_sub(masks{1})), max(CIs_sub(masks{2})), max(CIs_sub(masks{3}))];
    
    
   figure, imagesc(CIs_sub),title(sprintf('subject %s',name)),colorbar,colormap(jet)
   save(fullfile(res_path,filename), 'CIs_sub', 'mean_ROIs', 'max_ROIs','masks')%,'mACC','mRT');
end




