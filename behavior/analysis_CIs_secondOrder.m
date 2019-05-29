% As part of the frequency tagging / use of information project.
% 2nd order analysis on the use of facial information for emotion (happy vs
% neutral) and gender tasks.

% S. Faghel-Soubeyrand, may 2019


%%

task=2;

maskROI=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/interestregions.mat');
ff=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/maskellipse.mat');
images=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/cc_ims.mat');

im_avg=zeros(497,497);for im=1:264,im_avg=im_avg+rgb2gray(double(images.cc_ims{im}));end
im_avg=im_avg/264;


facemask=squeeze(double(ff.facemask(:,:,1)));
masque=logical(facemask);

path_all='/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/';
if task==1
task_flder='emotion/emotion_CIs/';
else
task_flder='gender/gender_CIs/';
end
fll_path=fullfile(strcat(path_all,task_flder));
data=dir(fullfile(fll_path,'*.mat'))
nSubs=length(data);
subCount=0;

% load CIs for every participant, correlate between tasks
clear CI_all all_use_ROI
for subject=1:nSubs
    subCount=subCount+1;
    SUB=load(fullfile(data(subject).folder,data(subject).name));
    
    CI_all(subCount,:,:)=SUB.CIs_sub;
    all_use_ROI(subCount,:)=SUB.mean_ROIs;
    %figure, imagesc(SUB.CIs_sub)
end
show_CI=squeeze(sum(CI_all)/sqrt(subCount));

figure, imagesc(show_CI.*masque),colorbar
hold on
contour(show_CI,[5 7 10],'--w')