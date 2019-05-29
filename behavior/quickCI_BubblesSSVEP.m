function [CI]=quickCI_BubblesSSVEP(path_data,name,task,blocks)


%load('/Users/labogosselin/Documents/LaboUmontreal/SSVEP_Bubbles/stimuli/cc_ims.mat');
ff=load('/home/adf/faghelss/Documents/project_FrequencyTaggingBubbles/behavior/maskellipse.mat');
facemask=squeeze(double(ff.facemask(:,:,1)));
masque=logical(facemask);
sizeX=size(masque,1);

nBins=10;
nTrials=100;
nTrialsPerBin   = nTrials/nBins;
thr             = .75; %Threshold for bubbles computation



p = .05;
tC = 3;
counter=0;
sigma=15;
counter=counter+1;
FWHM(counter) = sigma * sqrt(8*log(2));

[volumes, N] = CiVol(sum(masque(:)), 2);
[tP, k,tC1] = stat_threshold(volumes, N, FWHM, Inf, p, tC);  % the actual statistical tests
tP


n=0;
h = waitbar(0, 'ICs : 0 % complete');

CI_RTall=zeros(128);

%ims=load('images.mat'); % load images

% Bulles
sigma_bulles       = 12;
bulle=bubble(sigma_bulles);

params.ImSize         = 497;

mask=load ('maskellipse.mat');
facemask=logical(squeeze(double(mask.facemask(:,:,1))));
spaceSize=params.ImSize;

nblock=length(blocks)
for subject = 1
    X_all = [];
    y_all = [];
    RT_all = [];
    CID.DATA = [];
    contrast = [];
    
    n=n+1;
    b_boot_2D_ind=zeros(sizeX);
    b_2D_ind=zeros(sizeX);
    bub_stdev=12;
    
    
    for block = blocks,
        
        
        
        fname=fullfile(path_data,sprintf('Bubbles_SSVEP_%s_%s_%d.mat',name,task,block));
        load(fname); %e.g :
        %cid.data=cid.data;
        cid.DATA    = [CID.DATA cid.data];
        
        % Initialize vectors and matrices
        nTrials     = size(cid.data, 2);
        X           = zeros(nTrials, sizeX^2);
        y           = zeros(1, nTrials);
        
        % Get the seed from the cid and initialize the rand function
        %temp        = sscanf(cid.noise, '%s%s%s%s%d');
        %seed_0      = temp(end);

        rng(cid.seed)%initializing random number generator
        
        
        % Reproduce the noise for each trial and put it in a matrix
        for trial = 1:nTrials
            
            
            % Creation de bruit
            qteBulles = cid.data(5,trial);
            
            [posit,stim]=bubMask_spatial(qteBulles,spaceSize,facemask,bub_stdev);
            %             masque2D = (masque2D - min(masque2D(:)))/(max(masque2D(:)) - min(masque2D(:)));
            % mask= repmat(masque2D,[1 1 3]);% was repmat(masque2D, 1, 1, 3)
            %tmp = tmp2;
            tmp=posit;
            X(trial,:) = tmp(:);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%
%             X(trial,:) = (X(trial,:)-mean(X(trial,:)))/std(X(trial,:));
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Variable dependente
            y(1,trial)  = cid.data(9,trial); % ACCURACY
            RT(1,trial)  = cid.data(8,trial); % Response Time
            
            
            if find(isnan(X(trial,:)))
                X(trial,:)=zeros(1,sizeX*sizeX);
                %%%%%%%%%%%%%%
                % SOME PARTICIPANTS HAD 1 OR 0 BUBBLES,  made NAN WHEN WE DIVISE BY 0
                %%%%%%%%%%%%%%%
            end
            
        end
        % Standardisation de la VD
        y       = (y - mean(y)) / std(y);
        RT       = -1 * ((RT - mean(RT)) / std(RT));
        
        X_all = [X_all; X;];
        y_all = [y_all y];
        
        %RT_all = [RT_all RT];
        
        
%         b_part = y * X;
%         b_2D = reshape(b_part, sizeX, sizeX);
%         
%         b_2D_ind=b_2D_ind+b_2D;
%         % Bootstrapping the accuracry vector.        
%         index = ceil(size(y,2)*rand(size(y))); % avec remise, preferable
%         b_boot = y(index) * X;
%         b_boot_2D = reshape(b_boot, sizeX, sizeX);
%         
%         b_boot_2D_ind=b_boot_2D_ind+b_boot_2D;
        

    end
     h=waitbar(block/nblock, h, sprintf('Image de Classification %s : %3.2f %% complete', name,block/nblock));
     
    %On combine les blocs pour chq sujets..
    b_all = (y_all) * X_all;
    b_2D_all = reshape(b_all, sizeX, sizeX);
  % Bootstrapping the accuracry vector.
  
    index = randperm(size(y_all,2));
    b_boot = zscore(y_all(index)) * X_all;
    b_boot_2D_all = reshape(b_boot, sizeX, sizeX);
 
    mCi=mean(b_2D_all(:));
    ci=SmoothCi(b_2D_all-mCi,sigma)+mCi;
    
    mCi_boot=mean(b_boot_2D_all(:));
    ci_boot=SmoothCi(b_boot_2D_all-mCi_boot,sigma)+mCi_boot;
    
    % Standardize sci with sci_boot
    SSP = (ci - mean(ci_boot(:))) / std(ci_boot(:));
    

end
delete(h)

CI=SSP;

end


function bulle=bubble(bubStd)
% a single bubble 
maxHalfSize = 6*bubStd;
[y,x] = meshgrid(-maxHalfSize:maxHalfSize,-maxHalfSize:maxHalfSize);
bulle = exp(-(x.^2/bubStd^2)-(y.^2/bubStd^2));
clear x y tmp
end


function [posit,mask]=bubMask_spatial(qteBulles,spaceSize,masque,bubStd)
% Generate the bubbles mask
bulle=bubble(bubStd);
prob_tmp = qteBulles/sum(masque(:));
tmp=rand(spaceSize^2,1) .* masque(:);
tmp=reshape(tmp>=(1-prob_tmp),spaceSize,spaceSize); % makes the criteria probabilistic
posit=tmp;
masque2D=filter2(bulle,tmp);
masque2D = min(max(masque2D, 0), 1); % this is better

mask= repmat(masque2D,[1 1 3]);
end

