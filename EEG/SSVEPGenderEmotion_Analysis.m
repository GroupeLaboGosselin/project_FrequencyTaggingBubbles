%% SSVEP features analysis
%
% Simon Faghel-Soubeyrand,2016, last edit : feb 2019;

clear all


Subjects= {
    

% {'PaulLene','Paul_SSVEP_Gender.cnt','Gender','paul_lene_pilote',5}

{'RG','SSVEP_BubblesGender_Rose_16aout.cnt','Gender','RG',6} % there is only 19 trials for this task
{'RG','SSVEP_Bubbles_Emotions_RM_16_aout_2018.cnt','Emotions','RG',6}

{'LAET','SSVEPGender_Bubbles_Laetitia_10aout.cnt','Gender','Laetitia',6} % FACE TAG
{'LAET','SSVEPEmotions_Bubbles_Laetitia_10aout.cnt','Emotions','Laetitia',6} % FACE TAG

{'CAS','CarolineS_SSVEP_Gender_15aout.cnt','Gender','CarolineS',6} % FACE TAG
{'CAS','CarolineS_SSVEP_Emotion_15aout.cnt','Emotions','CarolineS',6} % FACE TAG

{'MEP','MelanieP_SSVEP_Gender14aout.cnt','Gender','MelanieP',6} % FACE TAG
{'MEP','MelanieP_SSVEP_Emotion_14aout.cnt','Emotions','MelanieP',6} % FACE TAG

{'MTQ','Mariethe_SSVEP_Bubbles_Gender_22aout.cnt','Gender','MTQ',6} %?	Probl�me avec l'�lectrode p3 durant le dernier essai du deuxi�me block. L'imp�dance est revenue � 256. En posant le casque on a eu des difficult�s avec les �lectrodes PO3 et POZ.
{'MTQ','Mariethe_SSVEP_Bubbles_Emotion_22aout.cnt','Emotions','MTQ',6} % FACE TAG

{'ER','EmilieR_SSVEP_Bubbles_Gender_21aout.cnt','Gender','ER',6} % N.B. Une image � fig�e durant le deuxi�me block de la t�che gender et l'exp�rience a bogu�e. On a �t� obliger de recommencer le block au complet.
{'ER','EmilieR_SSVEP_Bubbles_Emotion_21aout.cnt','Emotions','ER',6} % FACE TAG

{'YM','YanniM_SSVEP_Bubbles_Gender_24aout.cnt','Gender','YM',6} % FACE TAG
{'YM','YanniM_SSVEP_Bubbles_Emotion_24aout.cnt','Emotions','YM',6} % FACE TAG

{'SVT','SSVEP_Bubbles_Gender_SVT_31_aout_2018.cnt','Gender','SVT',6} % Aucun probl�me a priori
{'SVT','SVT_SSVEP_Bubbles_Emo_31aout.cnt','Emotions','SVT',6} % Aucun probl�me a priori

{'GB','GuillaumeB-SSVEP_Bubbles_Gender_29aout.cnt','Gender','GB',6} % probl�me �lectricit� + bug voir .docx
% {'GB','GuillaumeB_SSVEP_Bubbles_Gender1_29aout.cnt','Gender','GB',6} % probl�me �lectricit� + bug voir .docx
% {'GB','GuillaumeB_SSVEP_Bubbles_Gender2Bloc4et5_29aout.cnt','Gender','GB',6} % probl�me �lectricit� + bug voir .docx
{'GB','GuillaumeB_SSVEP_Bubbles_Emotion_29aout.cnt','Emotions','GB',6} %Manquait de concentration penant le 4eme bloc emotions.
% 
{'CG','CamilleG_SSVEP_Bubbles_Gender_20aout.cnt','Gender','CamilleG',6} %
{'CG','CamilleG_SSVEP_Bubbles_Emotion_20aout.cnt','Emotions','CamilleG',6} 

{'LC','LaurenceC_SSVEP_Bubbles_Gender_25aout.cnt','Gender','LC',6} % Pour Laurence l'image avait fig�e d�s qu'elle a r�pondu au premier visage du 6em bloc de la t�che �motion, donc j'ai pas eu le choix de recommencer le bloc. C'�tait le tout premier visage du 6em bloc.
{'LC','LaurenceC_SSVEP_Bubbles_Emotion_25aout.cnt','Emotions','LC',6} 


% Pour Laurence l'image avait fig�e d�s qu'elle a r�pondu au premier visage du 6em bloc de la t�che �motion, donc j'ai pas eu le choix de recommencer le bloc. C'�tait le tout premier visage du 6em bloc.

};

addpath /home/adf/faghelss/Documents/TOOLBOXES/fieldtrip-20190419
ft_defaults
%% Block 2 :  trial Creation
% http://www.fieldtriptoolbox.org/tutorial/networkanalysis

%res_path='
for subject=11%:22
    subject
    
    cfg                     = [];
    cfg.dataset             = Subjects{subject}{2};
    
    % % % % % % % For this here, you should put your path ! ! ! ! ! !
    cfg.dataFile            = ['/home/adf/faghelss/Documents/DATA/DATA_EEG/Data_SSVEPBubbles/SSVEP', cfg.dataset];
    % % % % % % % For this here, you should put your path ! ! ! ! ! !
    cfg.hdr               = ft_read_header(cfg.dataFile);
    cfg.trialdef.eventtype  = 'trigger';
    
    cfg.demean              = 'yes';
    cfg.trialdef.pre        = [.5 0];
    cfg.trialdef.post       = [69];
    
%     cfg.baselinewindow      = [-.400 -.100];   % baseline?
    cfg.baselinewindow      = [-.2 0];   % baseline?
    cfg.hpfilter            = 'yes';
    cfg.hpfiltertype        = 'but';% apply HIGHpass filter
    cfg.hpfreq              =  0.8;        % Highpass at X Hz
    cfg.lpfilter            = 'yes';     % apply' LOWpass filter
    cfg.lpfreq              =  80;      % lowpass at X Hz.
    cfg.hpfiltdf            =  5;
    cfg.dftfreq             =  [60]; %  Line noise removal (fourier)
    cfg.dftfilter           = 'yes';
    cfg.derivative          = 'no';
    cfg.trialfun            = @trialfunction_SSVEP_Features;%InterResponseSSVEP_trialfunctionWINK_GOODInterResponseSSVEP_trialfunctionBLINK%@InterResponseSSVEP_trialfunction @trialfunction_SSVEP_Features_Inverted;%InterResponseSSVEP_trialfunction
    cfg                     = ft_definetrial(cfg);
    
    %     Choosing which electrodes to work with ? maybe in freq analysis only
    all_electrodes=1:66;
    freq_electrodes_names = cfg.hdr.label(all_electrodes);
    freqcfg.channel=freq_electrodes_names;%freq_electrod
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %     RAW "PRE-PROCESSED" DATA.
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    cfg.reref         =  'yes';  
    cfg.refmethod     = 'avg';     % average referencing..look at Rellecke Shommer & Schacht (2013) for a discussion on rereferencing electrodes.
    cfg.refchannel    = {'M1','M2'};%;'all';'all'
    raw_data                = ft_preprocessing(cfg);
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %     WE RESAMPLE RAW DATA (1024Hz --> 256Hz) Take this part and
    %     recompute it if problem
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    cfg=[];
    cfg.resamplefs  =  256;%frequency at which the data will be resampled (default = 256 Hz)
    cfg.detrend     = 'no';%, detrend the data prior to resampling (no default specified, see below)
    [raw_data] = ft_resampledata(cfg, raw_data)
%     raw_data_1=raw_data
%     raw_data_2=raw_data;
%     raw_data_3=raw_data;
%     raw_data = ft_appenddata(cfg,raw_data_1,raw_data_2,raw_data_3)
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %     HERE, WE SAVE DATA
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        save(strcat('Trialed_avg_data_SSVEP_',Subjects{subject}{1},Subjects{subject}{3}),'raw_data')
end
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% PRE-PROCESSING DATA WITH RUN ICA.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Subj_cond2=[9 13 15 17 19 21 23 25]

Subj_cond1=[8 12 14 16 18 20 22 24]

clear freqcfg freq
for subject=61
    eval(sprintf('load(''Trialed_data%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{5}))
    
    cfg        = [];
    cfg.method = 'fastica';% runica  is the default and uses the implementation from EEGLAB
    comp = ft_componentanalysis(cfg, raw_data);
    
    eval(sprintf('save(''Comp_data%d_%s'',''raw_data'',''comp'')',subject,Subjects{subject}{5}))
end

%%
clear raw_data cfg comp
subject=61; % WATCH OUT FOR SUBJECT 19!!! rejected comp 11 but not sure..
eval(sprintf('load(''Comp_data%d_%s.mat'')',subject,Subjects{subject}{5}))

prepare layout
load('quickcap64.mat') % load 2D layout of electrodes

elec_montage = ft_read_sens('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc');
plot the components for visual inspection
figure
cfg = [];
cfg.component = [1:15];       % specify the component(s) that should be plotted
cfg.elec    = elec_montage; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)



cfg=[];
cfg.layout    = lay;
cfg.viewmode = 'component'
okok=ft_databrowser(cfg, comp)

[2 3 28]; %Subject 6 cond 2
%%

% remove the bad components and backproject the data
cfg = [];
all=1:66;
fftgood=[19,13,22,50,54,65,49,35,41,47,46,27,20,42];
all(fftgood)=[];
cfg.component = [1 4 9 10 14 6 2 3]; % to be removed component(s) THIS IS FOR INTER WINK
cfg.component = [1 14 6 10 3 9 24 17]; % to be removed component(s) THIS IS FOR INTER RESPONSE
cfg.component = [ 2 3 5]; % to be removed component(s) THIS IS FOR ALEX 12 trials WINK/BLINK
cfg.component =all; % to be removed component(s) THIS IS FOR 70 sec TRIAL 8 14 15
%
all(66)=[]
cfg.component = [2 17 43 51 52]; % Catherine blck1-3

cfg.component = [1 2 3 16 15 24:30 42:45 51 56 65] % FG COND1

cfg.component =   [3];




data_postica = ft_rejectcomponent(cfg, comp, raw_data)

cfg                     = [];
cfg.dataset             = Subjects{subject}{2};
cfg.dataFile            = ['/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/', cfg.dataset];
cfg.hdr                 = ft_read_header(cfg.dataFile);

electrodes=cfg.hdr.label(1:64);


%  eval(sprintf('save(''Trialed_PostIcadataHEAVYRej%s_%s'',''data_postica'',''comp'',''raw_data'',''cfg'')',Subjects{subject}{1},Subjects{subject}{5}))
% eval(sprintf('save(''Trialed_PostIcadataSMOOTHRej%d_%s'',''data_postica'',''comp'',''raw_data'',''cfg'')',subject,Subjects{subject}{5}))
%%  frequency of Components.
clear freqcfg
TrialTime=raw_data.time{1}(end);


hz_intermod=intermod([cell2mat(data{1}.freq)],[0.5 20]); % thre trial_features=cell2mat(data{trialBlckCounter}.freq);

freqcfg.method='mtmfft' ;%'wavelet';
freqcfg.tapsmofrq    = 1/TrialTime;
if strcmp(freqcfg.method,'mtmfft')
    freqcfg.pad = 'nextpow2';
    freqcfg.foilim       = [1 50];
    freqcfg.taper    = 'dpss';
else strcmp(freqcfg.method,'mtmconvol')
    freqcfg.pad = 'nextpow2';
    freqcfg.aper        = 'hanning';
    freqcfg.foi          = 2:freqcfg.tapsmofrq:20;                         % analysis 2 to 30 Hz in steps of 2 Hz
    freqcfg.t_ftimwin    = ones(length(freqcfg.foi),1).*0.5;   % length of time window = 0.5 sec
    freqcfg.toi          = -0.5:1:35;
end
electrodes=cfg.hdr.label(1:64)
posterior_elec = {'T7','T8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'};
freqcfg.channel=ft_channelselection({'all','-EOG'},posterior_elec);%,'-POz','-P5','-P6','-P1','-P2'},posterior_elec)%,'-M1','-M2','-AF7','-AF8','-AF3','-AF4','-F5','-F1','-F2','-FT7','-FT8','-F6','-FC3','-FCz','-FC4', '-Fp1', '-EOG'},freqcfg.channel)

[freq{subject}] = ft_freqanalysis(freqcfg,dataTime{1})%compraw_dataraw_data
load('quickcap64.mat') % load 2D layout of electrodes
freqcfg.layout=lay;


plotcfg.layout=lay;
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';

figure, plot(freq{subject}.freq, abs(freq{subject}.powspctrm));
xlabel('frequency (Hz)')
ylabel('phase-locked amplitude (a.u.)')



%%
%     
    cfg                     = [];
    cfg.dataset             = Subjects{1}{2};
    cfg.dataFile            = ['/Volumes/Simon Backup 2TB/DATA_EEG/', cfg.dataset];
    cfg.hdr                 = ft_read_header(cfg.dataFile);
    electrodes = ft_read_sens('/Volumes/Simon Backup 2TB/DATA_EEG/standard_waveguard64.elc');
    cfg.elec=electrodes;


for  subject=13:15
    clear behav
    taggs=[];
    counter=0;
    
    eval(sprintf('load(''Trialed70_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
%     eval(sprintf('load(''TrialedInter_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
    
    
    for block=1:Subjects{subject}{5}
        behav{block}=eval(sprintf('load(''/Volumes/Simon Backup 2TB/DATA_EEG/DATA_SSVEPBubbles/DataBehav/%s_SSVEP%s_%d.mat'')',Subjects{subject}{4},Subjects{subject}{3},block));
        for trial=1:4,counter=counter+1;for tag=1:3,taggs(counter,tag)=behav{block}.data{trial}.freq{tag};end;end
    end
    for trial=1:counter
        for ii=1:6
            if  taggs(trial,:)==Perm(ii,:)
                trial_perms{subject}(trial)=ii;
            end
        end
    end
    for whichPerm=1:6
        if ~any(trial_perms{subject}==whichPerm)
            freq{subject}{whichPerm}=freq{subject}{1};
            freq{subject}{whichPerm}.powspctrm=zeros(size(freq{subject}{whichPerm}.powspctrm));
            continue
        end
            
          
        freqcfg.baseline     = [];%[-.5 -.25];
        freqcfg.baselinetype = [];%'relchange'; %'db'
        TrialTime=raw_data.time{1}(end);
        freqcfg.trials=find(trial_perms{subject}==whichPerm);
        
        freqcfg.method='mtmfft' ;%'wavelet';
        freqcfg.tapsmofrq    = 1/TrialTime;
        %     freqcfg.baselinemethod='absolute' ;%'wavelet';
        %     freqcfg.baseline=[-1 .75];
        if strcmp(freqcfg.method,'mtmfft')
            freqcfg.pad = 'nextpow2';
            freqcfg.foilim       = [1 20];
            freqcfg.taper    = 'dpss';
            freqcfg.pad = 'nextpow2';
        else strcmp(freqcfg.method,'mtmconvol')
            freqcfg.pad = 'nextpow2';
            freqcfg.taper        = 'hanning';
            freqcfg.foi          = 1:freqcfg.tapsmofrq:20;                         % analysis 2 to 30 Hz in steps of 2 Hz
            freqcfg.t_ftimwin    = ones(length(freqcfg.foi),1).*0.5;   % length of time window = 0.5 sec
            freqcfg.toi          = 'all'%:1:35;
        end
        freqcfg.keeptrials='no';
        
        % electrodes=cfg.hdr.label(1:64)
        %     posterior_elec = {'T7','T8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'};
        %     freqcfg.channel=ft_channelselection({'all','-EOG','-M1','-M2'},posterior_elec);%,'-POz','-P5','-P6','-P1','-P2'},posterior_elec)%,'-M1','-M2','-AF7','-AF8','-AF3','-AF4','-F5','-F1','-F2','-FT7','-FT8','-F6','-FC3','-FCz','-FC4', '-Fp1', '-EOG'},freqcfg.channel)
        freqcfg.channel=ft_channelselection({'all','-EOG','-M1','-M2'},raw_data.label);
        % [freq] = ft_freqanalysis(freqcfg, raw_data)
        [freq{subject}{whichPerm}] = ft_freqanalysis(freqcfg,raw_data)%compraw_dataraw_datadata_postica
        freqcfg.elec=electrodes;
        
        
%                 figure, plot(freq{subject}{whichPerm}.freq,squeeze(mean(squeeze(abs(freq{subject}{whichPerm}.powspctrm(:,:,:))))));
%                 xlabel('frequency (Hz)')
%                 ylabel('phase-locked amplitude (a.u.)')
%                 legend(freqcfg.channel)
    end
end

perm=1;
Feat_taggs=[3.5294    8.0000    6.3158];
Perm=perms(Feat_taggs);
clear ALL_GOODindexes
hz_intermod=(intermod(Feat_taggs',[0.5 20]));
hz_intermod(8)=abs(2*Perm(perm,1)-2*Perm(perm,2));
hz_intermod(9)=abs(3*Perm(perm,1)-2*Perm(perm,2));
hz_intermod(10)=abs(2*Perm(perm,1)-3*Perm(perm,2));
all_taggs=[Feat_taggs Feat_taggs*2  hz_intermod hz_intermod*2];

all_taggs((all_taggs>20)&(all_taggs<.9))=[];
for ii=1:length(all_taggs)
    FeatFreq_indx = abs(bsxfun(@minus,all_taggs(ii)',freq{2}{whichPerm}.freq));
    [~,ALL_GOODindexes(ii)] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2))); 
end

noiseIndexes=1:size(FeatFreq_indx,2);
noiseIndexes=find(~ismember(noiseIndexes,ALL_GOODindexes));
IM_Indexes=find(ismember(noiseIndexes,ALL_GOODindexes));

%%

[grandavg_show{whichPerm}] = ft_freqgrandaverage(cfg, freq{3}{whichPerm})

    figure, plot(grandavg_show{whichPerm}.freq, abs(grandavg_show{whichPerm}.powspctrm));
    xlabel('frequency (Hz)')
    ylabel('phase-locked amplitude (a.u.)')
    legend(freqcfg.channel)
    hold on
    %     for ii=1:3 bar(Perm(whichPerm,ii),.03,0.1,'k'); end
    bar(Perm(whichPerm,1),.03,0.2, 'r')
    bar(Perm(whichPerm,2),.03,0.2, 'b')
    bar(Perm(whichPerm,3),.03,0.2, 'g')
    hold off
%%

for cond=1:2
% [freq_decscr] = ft_freqdescriptives(cfg, freq{1}{whichPerm})
for whichPerm=1:6
    if cond==1
    [grandavgOR{cond}{whichPerm}] = ft_freqgrandaverage(cfg, freq{2}{whichPerm},freq{4}{whichPerm},freq{6}{whichPerm},freq{8}{whichPerm},freq{10}{whichPerm},freq{12}{whichPerm},freq{14}{whichPerm})
    
    [grandavgOR{cond}{whichPerm}] = ft_freqgrandaverage(cfg, freq{2}{whichPerm})
    elseif cond==2
     [grandavgOR{cond}{whichPerm}] = ft_freqgrandaverage(cfg, freq{3}{whichPerm},freq{5}{whichPerm},freq{7}{whichPerm},freq{9}{whichPerm},freq{11}{whichPerm},freq{13}{whichPerm},freq{15}{whichPerm})
%     [grandavgOR{cond}{whichPerm}] = ft_freqgrandaverage(cfg, freq{3}{whichPerm},freq{5}{whichPerm},freq{7}{whichPerm})    
    end
    
%     figure, plot(grandavgOR{cond}{whichPerm}.freq, abs(grandavgOR{cond}{whichPerm}.powspctrm));
%     xlabel('frequency (Hz)')
%     ylabel('phase-locked amplitude (a.u.)')
%     legend(freqcfg.channel)
%     hold on
%     %     for ii=1:3 bar(Perm(whichPerm,ii),.03,0.1,'k'); end
%     bar(Perm(whichPerm,1),.03,0.2, 'r')
%     bar(Perm(whichPerm,2),.03,0.2, 'b')
%     bar(Perm(whichPerm,3),.03,0.2, 'g')
%     hold off
end
end
%%   Signal to Noise Ratio
% clear SNR_all
for whichPerm=1:6
    for subject=1:2
        
        %     [grandavgOR] = ft_freqgrandaverage(cfg, freq{subject}{whichPerm})
        % Left_eye=(Perm(2,1));
        % Right_eye=(Perm(2,2));
        % Mouth=(Perm(2,3));
        grandavg=grandavgOR{subject}{whichPerm};%freq{subject}{whichPerm};%grandavgOR;
        grandavg.powspctrm=grandavgOR{subject}{whichPerm}.powspctrm;%freq{subject}{whichPerm}.powspctrm;%grandavgOR.powspctrm;%freq{subject}.powspctrm;%-freq{1}.powspctrm;%abs(abs(grandavgOR.powspctrm));%-abs(grandavgBLINK.powspctrm))
        
        %     subject=2% grandavg.powspctrm(grandavg.powspctrm<0)=0;
        %  hz_intermod=(intermod([Perm(2,:)]',[0.5 20]))
        % grandavg=grandavgBLINK;
        clear SNR_1
        clear SNR_2 noise_pre noise_post SNR_1 SNR_2 SNR_Elec_powSpectrum
        
        interv=100;
        
        
        if (size(grandavg.powspctrm,1))>30
            nTrials=1;
            nb_elec=(size(grandavg.powspctrm,1));
        else
            nTrials=size(grandavg.powspctrm,1);
            nb_elec=(size(grandavg.powspctrm,2));
        end
        for trial= 1:nTrials
            for electrodes=1:nb_elec,
                
                %     Elec_powSpectrum=abs(grandavg.powspctrm(elec,:));
                if (size(grandavg.powspctrm,1))>30
                    Elec_powSpectrum=squeeze(grandavg.powspctrm(electrodes,:));
                    range_freq=max(size(Elec_powSpectrum));
                else
                    Elec_powSpectrum=squeeze(grandavg.powspctrm(trial,electrodes,:));
                    range_freq=max(size(Elec_powSpectrum));
                end
                
                
                for HzBin= interv+1:range_freq-(interv+2),
                    
                    noise_pre = (Elec_powSpectrum((HzBin-interv):(HzBin-10)));
                    noise_post= (Elec_powSpectrum((HzBin+10):(HzBin+interv)));
                    % end
                    noise= (mean(noise_pre)+mean(noise_post))/2;
                    Zs=(Elec_powSpectrum(HzBin)-mean([noise_pre; noise_post;]))/nanstd([noise_pre; noise_post;]);
                    
                    %                 SNR_1(HzBin)= snr(Elec_powSpectrum(HzBin),noise);
                    SNR_2(HzBin)=Zs;%Elec_powSpectrum(HzBin)/(noise);
                    
                    
                end
                SNR_Elec_powSpectrum(trial,electrodes,:)=SNR_2;
                %     SNR_1(SNR_1<0)=0;
                
            end
        end
        
        SNR_all{subject}{whichPerm}=SNR_Elec_powSpectrum;
        %     figure,
        figure, plot(grandavg.freq((1):size(SNR_2,2)),(squeeze((SNR_Elec_powSpectrum))));

        %     xlabel('frequency (Hz)')
        %     ylabel('SNR')
        %     hold on
        %     legend(grandavg.label)
        % bar(Left_eye(1),3,0.2, 'r')
        % % bar(Left_eye(2),1,0.1, 'r')
        % bar(Right_eye(1),3,.2, 'b')
        % % bar(Right_eye(2),1,0.1, 'b')
        % % bar(Mouth(1)+Right_eye(1),1,0.1, 'c')
        % % bar(Mouth(1)Right_eye(1),1,0.1, 'c')
        % bar(Mouth(1),3,0.2, 'g')
        % bar(Mouth(2),1,0.1, 'g')
        %  bar(Mouth(1)+Left_eye(1)+Right_eye(1),1,0.1, 'k')
        % for ii=1:2 bar(hz_intermod(ii),1,0.1,'m'); end
        % hold off
        
        % SNR{subject}=SNR_Elec_powSpectrum;
        % SNR_Elec_powSpectrum=SNR{1}-SNR{2};
        
        % NEED TO SEE SNROF DIFF ELECTRODES
%         
%             SNR_electrodes=freq{subject}{whichPerm};
%         
%             SNR_electrodes.powspctrm=SNR_Elec_powSpectrum;
%             SNR_electrodes.freq=freq{subject}{whichPerm}.freq(1:size(SNR_2,2));
%             figure,
%             for feat=1:3
%         
%                 plotcfg.channel=freqcfg.channel;
%                 plotcfg.interpolatenan ='yes';
%         
%                 plotcfg.interplimits = 'head';
%                 plotcfg.interpolation ='linear';%'nearest'
%                 plotcfg.colorbar='East';
%                 plotcfg.style='straight';
%                 plotcfg.layout='EEG1010.lay';
%                 plotcfg.parameter= 'powspctrm';
%                 plotcfg.showlabels ='yes';
%                 whichFreq=[Perm(whichPerm,feat)];%1.695%8%3.531;%6.312;%1.648%3.531% 6.312 8  %  9.8440 %8%hz_intermod(5);
%                 plotcfg.xlim=[whichFreq-0.00015 whichFreq+0.00015];
%         
%                 subplot(1,3,feat),ft_topoplotER(plotcfg,SNR_electrodes),colormap(jet)
%             end
    end
end

%%

plotcfg.channel=freqcfg.channel;
plotcfg.interpolatenan ='yes';

plotcfg.interplimits = 'head';
plotcfg.interpolation ='linear';%'nearest'
plotcfg.colorbar='no';%'East';
plotcfg.style='straight';
plotcfg.layout='EEG1010.lay';
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';
whichFreq=[Perm(whichPerm,1)];%1.695%8%3.531;%6.312;%1.648%3.531% 6.312 8  %  9.8440 %8%hz_intermod(5);
plotcfg.xlim=[whichFreq-0.00015 whichFreq+0.00015];
load ElectrodesMatrices.mat
hz_intermod=(intermod(Feat_taggs',[0.5 20]));
for ii=1:2 eyes_vs_mouth{ii}=zeros(62,1);end
for subject=[1]
    
    clear SNR_feat temp1 temp2
    whichPerms=1:6;
    nbFeats=3;
    figure,
    Feat_taggs_harm=[Feat_taggs Feat_taggs*2];
    for feat=1:nbFeats
        temp1{feat}=zeros(1,62);
        temp2{feat}=zeros(1,62);
        
        for whichPerm=whichPerms
            if any(isnan(SNR_all{subject}{whichPerm}))
                SNR_all{subject}{whichPerm}=zeros(size(SNR_all{subject}{whichPerm}));
            end
            if any(isnan(SNR_all{subject+1}{whichPerm}))
                SNR_all{subject+1}{whichPerm}=zeros(size(SNR_all{subject+1}{whichPerm}));
            end
            % finding fundamental frequency attached to trial feature and adding it to feature
            % topographical matrix
            Distag=Feat_taggs_harm;%intermod(Perm(whichPerm,:)',[.5 20]);
            Distag=Distag(feat)
            Freq_indx_1 = abs(bsxfun(@minus,Distag,freq{subject+1}{whichPerm}.freq(1:size(SNR_all{subject}{whichPerm}(:,:,:),3))));
            [~,index_freqTo] = min(Freq_indx_1(:,1:size(Freq_indx_1,2)));
            
            for trial=1:size(SNR_all{subject}{whichPerm},1)
                temp1{feat}=temp1{feat}+squeeze(SNR_all{subject}{whichPerm}(trial,:,index_freqTo));
            end
            for trial=1:size(SNR_all{subject+1}{whichPerm},1)
                temp2{feat}=temp2{feat}+squeeze(SNR_all{subject+1}{whichPerm}(trial,:,index_freqTo));
            end
            
            %     Freq_indx_1 = abs(bsxfun(@minus,Perm(whichPerm,feat)*2,freq{1}{whichPerm}.freq(1:size(SNR_Elec_powSpectrum(:,:,:),3))));
            %     [~,index_freqTo] = min(Freq_indx_1(:,1:size(Freq_indx_1,2)));
            %     for trial=1:size(SNR_all{1}{whichPerm},1)
            %     temp1{feat}=temp1{feat}+squeeze(SNR_all{1}{whichPerm}(trial,:,index_freqTo));
            %     end
            %     for trial=1:size(SNR_all{2}{whichPerm},1)
            %     temp2{feat}=temp2{feat}+squeeze(SNR_all{2}{whichPerm}(trial,:,index_freqTo));
            %     end
           
        end
        
        temp1{feat}=temp1{feat}/(numel(whichPerms)*size(SNR_all{subject}{whichPerm},1));
        temp2{feat}=temp2{feat}/(numel(whichPerms)*size(SNR_all{subject+1}{whichPerm},1));
        SNR_feat{feat}=freq{subject+1}{whichPerm};
        SNR_feat{feat}.freq=freq{subject+1}{whichPerm}.freq(1:max(size(SNR_all{subject}{whichPerm})));
        if size(SNR_all{subject}{whichPerm},1)>1
            for trial=1:size(SNR_all{subject}{whichPerm},1)
                SNR_feat{feat}.powspctrm(trial,1:62,1:1:length(SNR_feat{feat}.freq))=repmat(temp1{feat}',1,max(size(SNR_Elec_powSpectrum)),1);
            end
        else
            SNR_feat{feat}.powspctrm=repmat(temp1{feat}',1,max(size(SNR_Elec_powSpectrum)),1);
        end
        [h,p1]=ttest(temp1{feat}(:),temp2{feat}(:));
        
        subplot(3,nbFeats,feat),ft_topoplotER(plotcfg,SNR_feat{feat}),colormap(jet),title(sprintf(' emo vs gender : p=%.3f',p1))
        set(gca,'Clim',[2.5 8])
        if size(SNR_all{subject+1}{whichPerm},1)>1
            for trial=1:size(SNR_all{subject+1}{whichPerm},1)
                SNR_feat{feat}.powspctrm(trial,1:62,1:1:length(SNR_feat{feat}.freq))=repmat(temp2{feat}',1,max(size(SNR_all{subject}{whichPerm})),1);
            end
        else
            SNR_feat{feat}.powspctrm=repmat(temp2{feat}',1,max(size(SNR_all{subject}{whichPerm})),1);
        end
        subplot(3,nbFeats,feat+(nbFeats)),ft_topoplotER(plotcfg,SNR_feat{feat}),colormap(jet)
        set(gca,'Clim',[2.5 8])
        subplot(3,nbFeats,feat+(nbFeats*2)),bar([mean(temp1{feat}(:)) mean(temp2{feat}(:))])
        
        
    end
    
    elecs=find(elec_post_ind)-2;
%     eyes_vs_mouth{1}=(((temp1{1}(elecs))+(temp1{2}(elecs)))./2)./(temp1{3}(elecs));
%     eyes_vs_mouth{1}=eyes_vs_mouth{1}+(((temp1{1}(:))+(temp1{2}(:)))./2);
    eyes_vs_mouth{1}=eyes_vs_mouth{1}+(temp1{3}(:));
    eyes_vs_mouth{2}=eyes_vs_mouth{2}+(temp2{3}(:));
%     % eyes_vs_mouth{1}=((((temp1{1}(:)+temp1{4}(:))+(temp1{2}(:)+temp1{5}(:)))./2)./(temp1{3}(:)+temp1{6}(:)))./2;
%     eyes_vs_mouth{2}=(((temp2{1}(elecs))+(temp2{2}(elecs)))./2)./(temp2{3}(elecs));
    eyes_vs_mouth{2}=eyes_vs_mouth{2}+(((temp2{1}(:))+(temp2{2}(:)))./2);
%     eyes_vs_mouth{2}=eyes_vs_mouth{2}+(temp2{3}(:));%
    % eyes_vs_mouth{2}=((((temp2{1}(:)+temp2{4}(:))+(temp2{2}(:)+temp2{5}(:)))./2)./(temp2{3}(:)+temp2{6}(:)))./2;
%     
    [h,p]=ttest(eyes_vs_mouth{1}(:),eyes_vs_mouth{2}(:))
    figure, bar([mean(eyes_vs_mouth{1}) mean(eyes_vs_mouth{2})]),axis([0 3 .5 mean(eyes_vs_mouth{1}(:)+.5)]),title(sprintf(' p=%.2f',p))
% 
% %

% for feat=1:nbFeats
%     if size(SNR_all{2}{whichPerm},1)>1
%         for trial=1:size(SNR_all{2}{whichPerm},1)
%             SNR_feat{feat}.powspctrm(trial,1:62,1:1:length(SNR_feat{feat}.freq))=repmat(((temp2{feat}+temp2{feat})./2)',1,max(size(SNR_Elec_powSpectrum)),1);
%         end
%     else
%         SNR_feat{feat}.powspctrm=repmat(((temp2{feat}+temp2{feat})./2)',1,max(size(SNR_all{subject}{whichPerm}),1));
%     end
%        figure,  ft_topoplotER(plotcfg,SNR_feat{feat}),colormap(jet)
%        set(gca,'Clim',[1 5])
% end
figure, histogram(temp1{1})
end 
  %%
  
  
[h,p1]=ttest(temp1{1}(:),temp2{1}(:))
ys_emo=[mean(temp1{1}(:)) mean(temp1{2}(:)) mean(temp1{3}(:))];
std_emo=std([(temp1{1}(:)) (temp1{2}(:)) (temp1{3}(:))])/sqrt(length(temp2{1}(:)));
figure,
shadedErrorBar(1:3,ys_emo,std_emo,{'-or','markerfacecolor',[1,0.2,0.2]},1)
ys_emo=[mean(temp2{1}(:)) mean(temp2{2}(:)) mean(temp2{3}(:))];
std_emo=std([(temp2{1}(:)) (temp2{2}(:)) (temp2{3}(:))])/sqrt(length(temp2{1}(:)));
hold on
shadedErrorBar(1:3,ys_emo,std_emo,{'-ob','markerfacecolor',[0.2,.2,1]},1)

% figure, scatter(1:62,temp1{3}(:),'or')
% hold on
% scatter(1:62,temp2{3}(:),'ob')
%      
     %% Cluster tests
     
cfg = [];
cfg.channel          = 'all'
cfg.frequency        = 8;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 200;
% prepare_neighbours determines what sensors may form clusters
cfg.planarmethod = 'sincos';
cfg_neighb.method    = 'distance';
% cfg_neighb.method
cfg_neighb.layout    = 'EEG1010.lay';
cfg_neighb.channel= 'all';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freq{1}{1});

design = zeros(1,size(freq{2}{2}.powspctrm,1) + size(freq{1}{2}.powspctrm,1));
design(1,1:size(freq{2}{2}.powspctrm,1)) = 1;
design(1,(size(freq{2}{2}.powspctrm,1)+1):(2*size(freq{2}{2}.powspctrm,1)))=2; 
design

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, freq{1}{6}, freq{2}{6});

%% fourier on X electrodes, then see topo for each specific features
% RIGHT ELECTRODES OF INTEREST  64 (PO8), 62 (TP8) 18 (T8), 28(P8)  ::::: [64 62 18 28]
% Left ELECTRODES OF INTEREST  63 (PO8), 61 (TP7) 14 (T7), 24 (P7) ::::: [63 61 14 24]


load('FredPilote2_SSVEP_0_0_0.mat')


Subj_cond1=[ 8 12 14 16 18 20 22 24 26 28 30 33 35 38 41 43 45 48]

Subj_cond2=[ 9 13 17 19 21 23 25 27 29 31 36 39 46 49]

Subj_cond3=[32 34 37 40 42 44 47 50]


Subj_cond4=[51 53 55 57 59 61 65 68 71 74]
Subj_cond5=[52 54 56 58 60 62 66 69 72 75]

cfg                     = [];
cfg.dataset             = Subjects{1}{2};
cfg.dataFile            = ['/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/', cfg.dataset];
cfg.hdr                 = ft_read_header(cfg.dataFile);


clear freqcfg freq
for subject=Subj_cond1((end-5):end)
    %     eval(sprintf('load(''Trialed_PostIcadata%s_%s.mat'')',Subjects{subject}{1},Subjects{subject}{5}))
    %OR
    %     eval(sprintf('load(''Trialed_PostIcadataSMOOTHRej%d_%s'')',subject,Subjects{subject}{5}))
    
    eval(sprintf('load(''Trialed_data%s_%s.mat'')',Subjects{subject}{1},Subjects{subject}{5}))
    
    TrialTime=size(cell2mat(raw_data.time(1)),2)/raw_data.fsample;
    
    clear freq freqcfg
    harm=1;
    for frequency=1:6 % or frequency
        
        
        electrodes=[14:19 23:64]; %electrodes =[63 61 14 24];%[30 31 32];%[14 18 29 30 31 32 55 56 57 58 61 62];[29 30 31 32 55 56 57 58 59 60 61 62 65 66]%[55 56 57 58 59 60]% [14 18 29 30 31 32]%[29 30 31 32 55 56 57 58 59 60 61 62 65 66]%[30 31 32 61 62];%[30 31 32];%[29 30 31 32 55 56 57 58 63 64 19 66];    % nb of electrodes to include in analysis
        
        posterior_elec = {'T7','T8','P7','P3','Pz', 'P4', 'P8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8','Cz','CPz'};
        electrodes=[14:19 23:64]; %electrodes = Left tempo[63 61 14 24];% Occ [30 31 32];%
        
        freqcfg.channel=ft_channelselection({'all','-EOG','Cz','-M1','-M2','-AF7','-AF8','-AF3','-AF4','-F5','-F1','-F2','-F3','-F4','-F5','-F7','-F8','-FT7','-FT8','-F6','-FC3','-FCz','-FC4', '-Fp1', '-EOG'},posterior_elec)
        %     freqcfg.channel='all'%ft_channelselection({'all','-EOG'},posterior_elec)%,
        freqcfg.method= 'mtmfft';%'wavelet';
        
        if strcmp(freqcfg.method,'wavelet');
            cfg.pad = 'nextpow2';
            cfg.toi
        else
            if mod(frequency,2)==0
                
                freqcfg.foi   = cell2mat(data{1}.freq(harm))*2 %round(cell2mat(data{1}.freq(features)),2)*2 round(cell2mat(data{1}.freq(features)),2)*3
                harm=harm+1;
            else
                freqcfg.foi   = cell2mat(data{1}.freq(harm))
            end
            
            freqcfg.taper ='dpss';
            freqcfg.tapsmofrq  = 1/TrialTime;
        end
        load('quickcap64.mat') % load 2D layout of electrodes
        freqcfg.layout=lay;
        freqcfg.pad ='nextpow2';
        
        freqcfg.keeptrials='yes';
        [freq_F{subject}{frequency}] = ft_freqanalysis(freqcfg, raw_data)%data_postica % or CleanerData ( with no blinks and muscular artifacts
        
        %     figure,ft_topoplotER(plotcfg, freq_F{frequency})
    end
end
TrialTime=size(cell2mat(raw_data.time(1)),2)/1024;

%%

cfg                     = [];
cfg.dataset             = Subjects{2}{2};
cfg.dataFile            = ['/Volumes/Seagate Thin Plus/DATA_EEG/', cfg.dataset];
cfg.hdr                 = ft_read_header(cfg.dataFile);



posterior_elec = cfg.hdr.label(1:64);%{'T7','T8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'};
freqcfg.channel=ft_channelselection({'all','-EOG'},posterior_elec);%,'-POz','-P5','-P6','-P1','-P2'},posterior_elec)%,'-M1','-M2','-AF7','-AF8','-AF3','-AF4','-F5','-F1','-F2','-FT7','-FT8','-F6','-FC3','-FCz','-FC4', '-Fp1', '-EOG'},freqcfg.channel)
freqcfg.keeptrials='yes';
freqcfg.method='mtmfft' ;%'wavelet';
if strcmp(freqcfg.method,'mtmfft')
    freqcfg.pad = 'nextpow2';
    freqcfg.foilim       = [.8 80];
    freqcfg.taper    = 'dpss';
    freqcfg.pad = 'nextpow2';
end
 
h = waitbar(0,'Waiting for Fourier...');
% cd '/Volumes/Simon Backup 2TB/DATA_EEG/SSVEP_Bubbles_DATA/ProcessedData'
tasks={'Gender','Emotions'};
for subject=12:22

%     eval(sprintf('load(''/Volumes/Simon Backup 2TB/DATA_EEG/Trialed70_avg_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
    eval(sprintf('load(''/Volumes/Seagate Thin Plus/DATA_EEG/Data_SSVEPBubbles/SSVEP/Trialed_avg_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
    freqcfg.trials='all';

    TrialTime=raw_data.time{1}(end);
    
    freqcfg.tapsmofrq    = 1/TrialTime;
    
    [freq{subject}] = ft_freqanalysis(freqcfg,raw_data);
    Freq_sub=freq{subject};
    
%     eval(sprintf('save(''frequencySpctrm2_Response_SSVEP_%s_%d.mat'',''Freq_sub'',''freqcfg'')',tasks{mod(subject,2)+1},subject))
    eval(sprintf('save(''/Volumes/Seagate Thin Plus/DATA_EEG/Data_SSVEPBubbles/SSVEP/processedEEG/frequencySpctrm_SSVEP_%s_%d.mat'',''Freq_sub'',''freqcfg'')',tasks{mod(subject,2)+1},subject))
       waitbar(subject/15)
end
delete(h)


%%     hz_intermod=intermod([cell2mat(data{1}.freq)],[0.5 20]); % thre trial_features=cell2mat(data{trialBlckCounter}.freq);

% SNR FOR FEATURES



% Subj_cond2=[5 9 13 15 17 19 21 23 25 27]
% Subj_cond1=[5 9 13 15 17 19 21 23 25 27]
%  load SNR_subjectxtrialCond2_POSTICA.mat
CONDITION=2;

Subj_cond1=[ 8 12 14 16 18 20 22 24 26 28 30 33 35 38 41 43 45 48]

Subj_cond2=[ 9 13 17 19 21 23 25 27 29 31 36 39 46 49]

Subj_cond3=[32 34 37 40 42 44 47 50]
Subj=Subj_cond2;

clear SNR_1
clear SNR_2 noise_pre noise_post SNR_1 SNR_2 SNR_Elec_powSpectrum

% interval=[5:195]; % PO4 and PO6 seems to have better IM signal (e.g. 9.8)
interv=5;
% load('quickcap64.mat') % load 2D layout of electrodes
% freqcfg.layout=lay;


h = waitbar(0,'Waiting for SNR...');
%  Subj_cond1=[2 3 4 8 12 14 16 18]
%  Subj_cond2=[6 9 13 15 17 19]

subCounter=0;
for subject=Subj
    subCounter=subCounter+1;
    subject
    nTrials=size(freq{subject}.powspctrm,1);
    for trial=1:nTrials
         
        clear SNR_2 noise_pre noise_post SNR_1 Elec_powSpectrum
        
        nb_elec=size(freq{subject}.powspctrm,2);
        
        Mean_elec_powSpectrum= squeeze(mean(abs(freq{subject}.powspctrm(trial,:,:)),1));
        for electrodes=1:nb_elec,
            
            Elec_powSpectrum=squeeze(abs(freq{subject}.powspctrm(trial,electrodes,:)));
            
            for HzBin= interv+1:size(Elec_powSpectrum,1)-(interv+2),
                
                noise_pre = mean(Elec_powSpectrum((HzBin-interv):(HzBin-3)));
                noise_post= mean(Elec_powSpectrum((HzBin+3):(HzBin+interv)));
                noise= (noise_pre+noise_post)/2;
                SNR_2(HzBin)= Elec_powSpectrum(HzBin)/(noise);
                
            end
            SNR_Elec_powSpectrum{subject}(trial,electrodes,:)=SNR_2;
            
            
        end
        %                 figure, plot(freq{subject}.freq(1:size(SNR_2,2)), squeeze(SNR_Elec_powSpectrum{subject}(trial,:,:)));
        %                 xlabel('frequency (Hz)')
        %                 ylabel('SNR')
        %                 hold on
        %                 for ii=1:length(hz_intermod)  bar(hz_intermod(ii),1,0.1,'m'); end
        %                 hold off
        %
    end
    
    
    
    waitbar(subCounter/length(Subj))
end
close(h)
%%
% sound(R,40000);
switch (CONDITION)
    
    case 1
        eval(sprintf('save(''SNR_subjectxtrialCond%d_POSTICA.mat'',''SNR_Elec_powSpectrum'',''Subj_cond1'',''freq'') ',CONDITION))
    case 2
        eval(sprintf('save(''SNR_subjectxtrialCond%d_POSTICA_2.mat'',''SNR_Elec_powSpectrum'',''Subj_cond2'',''freq'') ',CONDITION))
    case 3
        eval(sprintf('save(''SNR_subjectxtrialCond%d.mat'',''SNR_Elec_powSpectrum'',''Subj_cond3'',''freq'') ',CONDITION))
end

% if CONDITION==1
%  eval(sprintf('save(''SNR_subjectxtrialCond%d_POSTICA.mat'',''SNR_Elec_powSpectrum'',''AllSub_Features'',''Subj_cond1'',''freq'') ',CONDITION))
% else
%  eval(sprintf('save(''SNR_subjectxtrialCond%d_POSTICA.mat'',''SNR_Elec_powSpectrum'',''AllSub_Features'',''Subj_cond2'',''freq'') ',CONDITION))
% end


cfg=[];
% cfg.layout=lay;
% for feature=1:3
%     figure,ft_topoplotER(plotcfg, Feature_topo{subject}{feature}), title(sprintf('%s',featureLabel{feature}));
% end


%% SNR FEATURE TOPOGRAPHY

% Really good intra-subject reliability.. suggest veryvery low noise.

featureLabel={'Left eye', 'Right eye', 'Mouth','Right versus Left Eye','Right eye versus Mouth','Both eyes vs Mouth'};
electrodes = ft_read_sens('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc');
freqcfg.elec=electrodes;

Subj_cond1=[ 12 14 16 18 20 22 24 26 28 30 33 35 38 41 43 45 48] % 8 has only 18 electrodes..

Subj_cond2=[ 9 13 17 19 21 23 25 27 29 31 36 39 46 49]

Subj_cond3=[32 34 37 40 42 44 47 50]

SUBJ_All=[Subj_cond1 Subj_cond2 Subj_cond3];
% all 3 conditions;
Subj_cond1_all=[30   35    38    41    43    45     48];

Subj_cond2_all=[31   36    39    13        46    49];

Subj_cond3_all=[32   37    40    42    44    47    50];

Subj_Attended_all=[9 13 17 19 21 25 29 31 36 39 46 49]; % check out 15 and 17, 25 for mismatch between conditions
% Subj_Attended_all=[8 12 14 16 18 20 28]; % check out 15 and 17, 25 for mismatch between conditions

for ii=1:5
    
    Feature_topoALL{ii}=AllSub_Features{1};
    Feature_topoALL{ii}.label=freq{BOGUSFREQ}.label;
    Feature_topoALL{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
    LeftATT_topoALL{ii}=Feature_topoALL{ii};
    RightATT_topoALL{ii}=Feature_topoALL{ii};
    
end
%%

% clear freq
CONDITION=2;
% counter=0;
% switch (CONDITION)
%     case 1
%     %      load SNR_subjectxtrialCond1.mat
%     load SNR_subjectxtrialCond1_POSTICA.mat
%     %     Subj_cond1=[2 3 4 8 12 14 16 18 20 22 24]
% %     Subj_cond1=[ 12 14 16 18 20 22 24 26 28] % Cond 1 only left out 14
%     Subj_cond1=[ 8 12 14 16 18 20 22 24 26 28 30 33 35 38 41 43 45 48]
%
%     %     Subj_cond1=[8 12 16 20 22 24 26 28] % 14 not sure of tag, 18 has very low SNR
%     %     Subj_cond1=[24] % NEED TO CHECK CHLOE's (20) comp data.
% %     Subj=Subj_cond1((end-7):end);
% Subj=Subj_cond1_all;
%     case 2
% %     load SNR_subjectxtrialCond2_POSTICA.mat
%  load SNR_subjectxtrialCond2_POSTICA_2.mat
%     %     load SNR_subjectxtrialCond2_POSTICA.mat
% %     Subj_cond2=[ 13 15 17 19 21 23 25 27 29] % Cond 2 only left out 15 PROBLEM WITH 27 and 29 subject
%     Subj_cond2=[ 9 13 17 19 21 23 25 27 29 31 36 39 46 49]
%     %    Subj_cond2=[ 9 13 17 21 23 25 27 29]
%     %     Subj_cond2=25; % FSAVOI (19) weird, CATHERINE (9) same pattern as COND1 (8), 15 a bit weird
% %     Subj=Subj_cond2;
% Subj=Subj_cond2_all;
%     case 3
%         load SNR_subjectxtrialCond3.mat
% %         Subj_cond3(3)=[];
% %         Subj_cond3(3)=[];
% %         Subj=Subj_cond3;
%         Subj=Subj_cond3_all;
%
% end

Subj=Subj_Attended_all%Subj_cond1;
BOGUSFREQ=Subj(1);

nSubjects=length(Subj);


rng('default')
clear plotcfg
featureLabel={'Left eye', 'Right eye', 'Mouth','Left versus Right Eye','Right eye versus Mouth','Both eyes vs Mouth'};
plotcfg.channel= freq{BOGUSFREQ}.label;
plotcfg.colorbar='East';
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='East';
load('quickcap64.mat') % load 2D layout of electrodes
plotcfg.elec=electrodes;
plotcfg.lay=lay;

plotcfg.showlabels ='yes';


for subject=Subj%13%Subj
    Lcounter=0;
    Rcounter=0;
    counter=counter+1;
    for ii=1:3
        
        Feature_topo{subject}{ii}=AllSub_Features{1};
        Feature_topo{subject}{ii}.label=freq{BOGUSFREQ}.label;
        Feature_topo{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
        Right_eyeAttended_topo{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
        Left_eyeAttended_topo{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
    end
    
    
    blockCounter=1;
    trialBlckCounter=0;
    
    nTrials=size(SNR_Elec_powSpectrum{subject},1);
    if subject==39
        nTrials=12;
    end
    for trial=1:nTrials
        
        if subject==22
            eval(sprintf('load(''%s_blck%d_SSVEP_0_1_1.mat'') ',Subjects{subject}{3} ,Subjects{subject}{4}(blockCounter)))
            
        else
            if ismember(subject,Subj_cond1)
                CONDITION=1;
            elseif ismember(subject,Subj_cond2)
                CONDITION=2;
            elseif ismember(subject,Subj_cond3)
                CONDITION=3;
            end
            
            %             if CONDITION==1
            eval(sprintf('load(''%s_blck%d_SSVEP_0_1_%d.mat'') ',Subjects{subject}{3} ,blockCounter,CONDITION))
            %             else
            %                 eval(sprintf('load(''%s_blck%d_SSVEP_0_1_2.mat'') ',Subjects{subject}{3} ,blockCounter,))
            %             end
        end
        
        if  trial== 5 || trial== 9 || trial== 13 || trial== 17
            blockCounter=blockCounter+1
            
            if subject==22 && CONDITION==1
                eval(sprintf('load(''%s_blck%d_SSVEP_0_1_1.mat'') ',Subjects{subject}{3} ,Subjects{subject}{4}(blockCounter)))
            else
                if ismember(subject,Subj_cond1)
                    CONDITION=1;
                elseif ismember(subject,Subj_cond2)
                    CONDITION=2;
                elseif ismember(subject,Subj_cond3)
                    CONDITION=3;
                end
                
                %                 if CONDITION==1
                eval(sprintf('load(''%s_blck%d_SSVEP_0_1_%d.mat'') ',Subjects{subject}{3} ,blockCounter,CONDITION))
                %                 else
                %                     eval(sprintf('load(''%s_blck%d_SSVEP_0_1_2.mat'') ',Subjects{subject}{3} ,blockCounter))
                %                 end
            end
            trialBlckCounter=1;
        end
        
        
        for feature=1:3
            if mod(trial,4)~=0
                trial_feature(feature)=cell2mat(data{mod(trial,4)}.freq(feature));
                harmon_feat(feature)=cell2mat(data{mod(trial,4)}.freq(feature))*2;
            else
                trial_feature(feature)=cell2mat(data{4}.freq(feature));
                harmon_feat(feature)=cell2mat(data{4}.freq(feature))*2;
            end
            
            % finding fundamental frequency attached to trial feature and adding it to feature
            % topographical matrix
            FeatFreq_indx = abs(bsxfun(@minus,trial_feature(feature)',freq{BOGUSFREQ}.freq));
            [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
            %             Feature_topo{subject}{feature}.powspctrm(:)=Feature_topo{subject}{feature}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
            Feature_topo{subject}{feature}.powspctrm(:)=Feature_topo{subject}{feature}.powspctrm(:)+squeeze(freq{subject}.powspctrm(trial,:,index_freqToFeatures))';
            if CONDITION==2 && strcmp(Subjects{subject}{5+blockCounter},'Right')
                if feature==1
                    Rcounter=Rcounter+1;
                end
                %                 Right_eyeAttended_topo{subject}{feature}.powspctrm(:)=Right_eyeAttended_topo{subject}{feature}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
                Right_eyeAttended_topo{subject}{feature}.powspctrm(:)=Right_eyeAttended_topo{subject}{feature}.powspctrm(:)+squeeze(freq{subject}.powspctrm(trial,:,index_freqToFeatures))';
                
            elseif CONDITION==2 && strcmp(Subjects{subject}{5+blockCounter},'Left')
                if feature==1
                    Lcounter=Lcounter+1;
                end
                %                 Left_eyeAttended_topo{subject}{feature}.powspctrm(:)=Left_eyeAttended_topo{subject}{feature}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
                Left_eyeAttended_topo{subject}{feature}.powspctrm(:)=Left_eyeAttended_topo{subject}{feature}.powspctrm(:)+squeeze(freq{subject}.powspctrm(trial,:,index_freqToFeatures))';
            end
            %            Feature_pow{feature}(counter,trial,:)=squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
            % same with first harmonic
            FeatHarmFreq_indx = abs(bsxfun(@minus,harmon_feat(feature),freq{BOGUSFREQ}.freq));
            [~,index_freqToFeatures] = min(FeatHarmFreq_indx(:,1:size(FeatHarmFreq_indx,2)));
            Feature_topo{subject}{feature}.powspctrm(:)=Feature_topo{subject}{feature}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
            
            if CONDITION==2 && strcmp(Subjects{subject}{5+blockCounter},'Right')
                
                Right_eyeAttended_topo{subject}{feature}.powspctrm(:)=Right_eyeAttended_topo{subject}{feature}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
            elseif CONDITION==2 && strcmp(Subjects{subject}{5+blockCounter},'Left')
                
                Left_eyeAttended_topo{subject}{feature}.powspctrm(:)=Left_eyeAttended_topo{subject}{feature}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
                
            end
            %
            
        end
        
        
    end
    
    
    
    %     for ii=1:3,Stats_all{ii}.powspctrm(counter,:)=Feature_topo{subject}{ii}.powspctrm(:)/nTrials;end% endRight_eyeAttended_topo{subject}{ii}.powspctrm(:)/(Rcounter/3);
    for ii=1:3,Stats_all{ii}.powspctrm(counter,:)=Right_eyeAttended_topo{subject}{ii}.powspctrm(:)/(Rcounter);end
    %
    for feature=1:3,
        Feature_topo{subject}{feature}.powspctrm(:)=Feature_topo{subject}{feature}.powspctrm(:)/nTrials;
        Right_eyeAttended_topo{subject}{feature}.powspctrm(:)=Right_eyeAttended_topo{subject}{feature}.powspctrm(:)/(Rcounter);
        Left_eyeAttended_topo{subject}{feature}.powspctrm(:)=Left_eyeAttended_topo{subject}{feature}.powspctrm(:)/(Lcounter);
        Feature_topoALL{feature}.powspctrm(:)=Feature_topoALL{feature}.powspctrm(:)+Feature_topo{subject}{feature}.powspctrm(:);
        LeftATT_topoALL{feature}.powspctrm(:)=LeftATT_topoALL{feature}.powspctrm(:)+Left_eyeAttended_topo{subject}{feature}.powspctrm(:);
        RightATT_topoALL{feature}.powspctrm(:)=RightATT_topoALL{feature}.powspctrm(:)+Right_eyeAttended_topo{subject}{feature}.powspctrm(:);
        
    end
    
    
    
    
end

for feature=1:3,
    Feature_topoALL{feature}.powspctrm= Feature_topoALL{feature}.powspctrm/nSubjects;
    LeftATT_topoALL{feature}.powspctrm(:)=LeftATT_topoALL{feature}.powspctrm(:)/nSubjects;
    RightATT_topoALL{feature}.powspctrm(:)=RightATT_topoALL{feature}.powspctrm(:)/nSubjects;
end
cfg=[];
%

plotcfg.highlight='labels';
plotcfg.highlightchannel={'PO8','PO7','OZ','O2','O1'}
plotcfg.colormap ='jet';
plotcfg.highlightcolor=[1 0 0]
plotcfg.interplimits       = [];
plotcfg.interpolation      = 'v4';
plotcfg.colorbar='East'%'East';
plotcfg.zlim      =  [];%[0 23];
plotcfg.lay=lay;
%
featureLabel={'Left eye', 'Right eye', 'Mouth','Left eye: Left vs Right attended','Right eye:  Right vs Left  attended','Both eyes vs Mouth'};
Feature_topoALL{4}=LeftATT_topoALL{1};
Feature_topoALL{5}=LeftATT_topoALL{1};

Feature_topoALL{4}.powspctrm=LeftATT_topoALL{1}.powspctrm-RightATT_topoALL{1}.powspctrm;
Feature_topoALL{5}.powspctrm=RightATT_topoALL{2}.powspctrm-LeftATT_topoALL{2}.powspctrm;

% Feature_topoALL{4}.powspctrm(:)=Feature_topoALL{1}.powspctrm(:)-Feature_topoALL{2}.powspctrm(:);
% Feature_topoALL{5}.powspctrm(:)=((Feature_topoALL{1}.powspctrm(:)+ Feature_topoALL{2}.powspctrm(:))/2)-Feature_topoALL{3}.powspctrm(:);
% figure, see Morgan et al., 1996 for attention related SSVEP
plotcfg.zlim      =  [];%[0 16];
plotcfg.style =[];%'straight'
%%
plotcfg.colormap ='jet';%[];%'straight'
for feature=5
    %     figure,ft_topoplotER(plotcfg, LeftATT_topoALL{feature}), title(sprintf('%s',featureLabel{feature}));
    %         colormapeditor
    %     figure,ft_topoplotER(plotcfg, RightATT_topoALL{feature}), title(sprintf('%s',featureLabel{feature}));
    % %     colormapeditor
    figure,ft_topoplotER(plotcfg, Feature_topoALL{feature}), title(sprintf('%s',featureLabel{feature}));
    
    figure, ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),(Feature_topoALL{feature}.powspctrm)),colorbar
    colormap('jet')
    % %     colormapeditor
    %     figure, ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),(RightATT_topoALL{feature}.powspctrm)),colorbar
    %     colormap('jet')
    %     colormapeditor
    %     figure, ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),(LeftATT_topoALL{feature}.powspctrm)),colorbar
    %     colormap('jet')
    %     colormapeditor
    %     subplot(3,1,feature)
end
%, title(sprintf('%s',interLabel{feature}))
%  figure, ft_topoplotER(plotcfg, Feature_topoALL{feature}), title(sprintf('%s',featureLabel{feature})),colorbar
% readlocs('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc')
% headplot('setup',elec_montage.elecpos(elec_MoneTwo_ind==0,:,:,:), 'PAUL_SSVEP+epochs.spl');
% figure, ft_plot_topo3d(elec_montage.elecpos(elec_MoneTwo_ind==0,:,:,:),Feature_topoALL{2}.powspctrm-Feature_topoALL{3}.powspctrm,'refine',5),colorbar
%%
plotcfg.colormap ='jet';
plotcfg.zlim      = 'maxmin';
plotcfg.colorbar='East';
figure,ft_topoplotER(plotcfg, Feature_topoALL{4}), title(sprintf('%s',featureLabel{4}));
figure,ft_topoplotER(plotcfg, Feature_topoALL{5}), title(sprintf('%s',featureLabel{6}));
figure, ft_plot_topo3d(elec_montage.elecpos(elec_MoneTwo_ind==0,:,:,:),Feature_topoALL{5}.powspctrm,'refine',5),colorbar
%  save('SNRFeatures_subjectxtrialCond1.mat','Feature_topo','Feature_topoALL')
%% Feature  responses cond 1 vs 2


COND1.showFeat{1}=(COND1.featureResp{1}.powspctrm(:,find(elec_post_ind)));
COND1.showFeat{2}=(COND1.featureResp{2}.powspctrm(:,find(elec_post_ind)));
COND1.showFeat{3}=(COND1.featureResp{3}.powspctrm(:,find(elec_post_ind)));

COND2.showFeat{1}=(COND2.featureResp{1}.powspctrm(:,find(elec_post_ind)));
COND2.showFeat{2}=(COND2.featureResp{2}.powspctrm(:,find(elec_post_ind)));
COND2.showFeat{3}=(COND2.featureResp{3}.powspctrm(:,find(elec_post_ind)));


figure, bar([mean(COND1.showFeat{1}(:)) mean(COND2.showFeat{1}(:)) mean(COND1.showFeat{2}(:)) mean(COND2.showFeat{2}(:)) mean(COND1.showFeat{3}(:)) mean(COND2.showFeat{3}(:))])
set(gca,'XTick',1:6,'XTickLabel',{'LeftE Cond1','LeftE Cond2','RightE Cond1','RightE Cond2','Mouth Cond1','Mouth Cond2'})
%%
% set(gca,'XTick',1:3,'XTickLabel',{'Left hemisphere electrodes','central','Right hemisphere'})
% % which ELECTRODES : elec_post_ind
% for ii=1:3,
% [h,p]=ttest2(COND1.featureResp{ii}.powspctrm(:,find(elec_post_ind)),COND2.featureResp{ii}.powspctrm(:,find(elec_post_ind)))
% end


% INTERFeature  responses cond 1 vs 2 %%%
% for ii=1:7,
% [h,p]=ttest2(COND1.INTERfeatureResp{ii}.powspctrm(:,find(elec_post_ind)),COND2.INTERfeatureResp{ii}.powspctrm(:,find(elec_post_ind)))
% end

COND1.showInterFeat{1}=(COND1.INTERfeatureResp{1}.powspctrm(:,find(elec_post_ind)));%+COND1.INTERfeatureResp{2}.powspctrm(:,find(elec_post_ind)))/2;
COND1.showInterFeat{2}=(COND1.INTERfeatureResp{2}.powspctrm(:,find(elec_post_ind)));%+COND1.INTERfeatureResp{4}.powspctrm(:,find(elec_post_ind)))/2;
COND1.showInterFeat{3}=(COND1.INTERfeatureResp{3}.powspctrm(:,find(elec_post_ind)));%+COND1.INTERfeatureResp{6}.powspctrm(:,find(elec_post_ind)))/2;
COND1.showInterFeat{4}=(COND1.INTERfeatureResp{4}.powspctrm(:,find(elec_post_ind)));
COND1.showInterFeat{5}=(COND1.INTERfeatureResp{5}.powspctrm(:,find(elec_post_ind)));


COND2.showInterFeat{1}=(COND2.INTERfeatureResp{1}.powspctrm(:,find(elec_post_ind)));%+COND2.INTERfeatureResp{2}.powspctrm(:,find(elec_post_ind)))/2;
COND2.showInterFeat{2}=(COND2.INTERfeatureResp{2}.powspctrm(:,find(elec_post_ind)));%+COND2.INTERfeatureResp{4}.powspctrm(:,find(elec_post_ind)))/2;
COND2.showInterFeat{3}=(COND2.INTERfeatureResp{3}.powspctrm(:,find(elec_post_ind)));%+COND2.INTERfeatureResp{6}.powspctrm(:,find(elec_post_ind)))/2;
COND2.showInterFeat{4}=(COND2.INTERfeatureResp{4}.powspctrm(:,find(elec_post_ind)));
COND2.showInterFeat{5}=(COND2.INTERfeatureResp{5}.powspctrm(:,find(elec_post_ind)));

figure, bar([mean(COND1.showInterFeat{1}(:)) mean(COND2.showInterFeat{1}(:)) mean(COND1.showInterFeat{2}(:)) mean(COND2.showInterFeat{2}(:)) mean(COND1.showInterFeat{3}(:)) mean(COND2.showInterFeat{3}(:))])
set(gca,'XTick',1:6,'XTickLabel',{'Both Cond1','Both Cond2','Mouth&Left Cond1','Mouth&Left Cond2','Mouth&Right Cond1','Mouth&Right Cond2'})

figure, bar([mean(COND1.INTERfeatureResp{4}.powspctrm(:)) mean(COND2.INTERfeatureResp{4}.powspctrm(:))])
set(gca,'XTick',1:2,'XTickLabel',{'three Cond1','three Cond2'})


figure, bar([mean(COND1.showInterFeat{1})-(mean(COND1.showInterFeat{2}(:))+mean(COND1.showInterFeat{3}(:)))/2 mean(COND2.showInterFeat{1}(:))-(mean(COND2.showInterFeat{2}(:))+mean(COND2.showInterFeat{3}(:)))/2])
set(gca,'XTick',1:2,'XTickLabel',{'Both Cond1','Both Cond2'})
%%   !!! Quantification of SNR for Features and Intermodulate frequencies !!!!!
%load ElectrodesMatrices.mat
load electrodes_indexes.mat

tagg= [ 6.3158    8.0000    3.5294];
Perm = perms(tagg);
tasks={'Gender','Emotions'};

% clear COND_IM
Subj=1:22%%2:23

Fs=256;
N=5;

nbSub=length(Subj);
% PERM_powspctrm=nan(nbSub,6,64,10139);
PERM_powspctrm=nan(nbSub,6,64,8);
% size(freq{Subj(1)}.powspctrm)
clear PERM_powspctrm ANOVA_COND IM_all PERM
subCounter=0;
%cd '/home/adf/faghelss/Documents/DATA/DATA_EEG/Data_SSVEPBubbles/SSVEP/processedEEG'
path_to_freq='/home/adf/faghelss/Documents/DATA/DATA_EEG/Data_SSVEPBubbles/SSVEP/processedEEG/';
for subject=Subj
    subCounter=subCounter+1
   
  load(strcat(path_to_freq,'frequencySpctrm2_SSVEP_',tasks{mod(subject,2)+1},'_',num2str(subject),'_.mat'),'Freq_sub','freqcfg')
  freq{subject}=Freq_sub;

    counter=0;
    for block=1:Subjects{subject}{5}
        behav{block}=eval(sprintf('load(''/home/adf/faghelss/Documents/DATA/DATA_EEG/Data_SSVEPBubbles/DataBehav/%s_SSVEP%s_%d.mat'')',Subjects{subject}{4},Subjects{subject}{3},block));
        for trial=1:4,counter=counter+1;for tag=1:3,taggs(counter,tag)=behav{block}.data{trial}.freq{tag};end;end
    end
    for trial=1:counter
        for ii=1:6
            if  taggs(trial,:)==Perm(ii,:)
                TrialTagg{subject}(trial)=ii;
            end
        end
    end
  
    
    clear permidx cfg.trials
    for perm=1:6
        PERM{subject}{perm}=[];
        if any(TrialTagg{subject}==perm)
            cfg.trials=find(TrialTagg{subject}==perm);%'all'
            cfg.trials(cfg.trials>min(size(freq{subject}.powspctrm)))=[];
            
        else
            continue
        end
        
        
        PERM{subject}{perm}=cfg.trials;
        if any(PERM{subject}{perm})
        clear datum
        if size(PERM{subject}{perm},2)>1
%         datum=squeeze((freq{subject}.powspctrm(PERM{subject}{perm},:,1:size(freq{Subj(1)}.powspctrm,3)))).*trial_perfo{subject}(PERM{subject}{perm})';
        datum=squeeze((freq{subject}.powspctrm(PERM{subject}{perm},:,1:size(freq{Subj(1)}.powspctrm,3))));
        
        PERM_powspctrm(subCounter,perm,:,:)=squeeze(nanmean(datum));
        else
%         datum=squeeze(freq{subject}.powspctrm(PERM{subject}{perm},:,(1:size(freq{Subj(1)}.powspctrm,3)))).*trial_perfo{subject}(PERM{subject}{perm})';
        datum=squeeze(freq{subject}.powspctrm(PERM{subject}{perm},:,(1:size(freq{Subj(1)}.powspctrm,3))));
        PERM_powspctrm(subCounter,perm,:,:)=squeeze((datum));
        end
        else
        PERM_powspctrm(subCounter,perm,:,:)=nan(squeeze(size(freq{subject}.powspctrm(1,:,1:size(freq{Subj(1)}.powspctrm,3)))));
        end

    end
    %figure, plot(freq{Subj(1)}.freq,(squeeze((PERM_powspctrm(subject,3,electrodesAll,:)))))%,axis([ 2 20 0 .15])
  
end
% *trial_perfo{subject}
electrodesAll=find(elec_post_ind);



%%

clear ALL_GOODindexes
tagg= [ 6.3158    8.0000    3.5294];
Perm = perms(tagg);

hz_intermod=(intermod(tagg',[0.5 50]))
perm=1;
    hz_intermod(8)=abs(2*Perm(perm,1)-2*Perm(perm,2));
    hz_intermod(9)=abs(3*Perm(perm,1)-2*Perm(perm,2));
    hz_intermod(10)=abs(2*Perm(perm,1)-3*Perm(perm,2));
all_taggs=[tagg 6.3158*2    8.0000*2    3.5294*2  hz_intermod hz_intermod*2];

all_taggs((all_taggs>50)&(all_taggs<.9))=[];
for ii=1:length(all_taggs)

            FeatFreq_indx = abs(bsxfun(@minus,all_taggs(ii)',freq{subject}.freq));
            [~,ALL_GOODindexes(ii)] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));

end


noiseIndexes=1:size(FeatFreq_indx,2);
noiseIndexes=find(~ismember(noiseIndexes,ALL_GOODindexes));
IM_Indexes=find(ismember(noiseIndexes,ALL_GOODindexes));
randperm(length(noiseIndexes));
%choosing random noise index for each trial
% noiseIndexes(round(rand(1,1).*2400));
%% Fourier power spectrum
% Subjs_cond=[2:2:22] % or [2:2:22]-1
Emotion_task_sub=[2:2:22];
Gender_task_sub=[1:2:22]
Subjs_cond=[2:22] % or [2:2:22]-1

load electrodes_indexes.mat
elecs=find(elec_post_ind);
elecs(17:18)=[]; % throw out TP7,TP8
elecs(4:5)=[]; % throw out P3,P4
electrodesAll=elecs;

figure,
 for perm=1:6
     hz_intermod=intermod(Perm(perm,:)',[1 20]);
  subplot(3,2,perm),plot(freq{Subj(1)}.freq,(squeeze(nanmean(PERM_powspctrm(Subjs_cond,perm,electrodesAll,:))))),axis([ 2 20 0 .15])
%  subplot(3,2,perm),plot(freq{Subj(1)}.freq,squeeze(nanmean(nanmean(PERM_powspctrm(:,perm,electrodesAll,:)),3)),'-k'),axis([ 2 20 0 .2]),title(sprintf('permutation %d',perm))
         hold on
         bar(Perm(perm,1) ,.015,.4,'r'),
         hold on
         bar(Perm(perm,2),.015,.4,'b')
         hold on
         bar(Perm(perm,3),.015,.4,'g')
         if perm==1
             legend(freq{22}.label(electrodesAll))
         %legend({'uV','Left-eye frequency tag','right-eye tag','mouth tag'})
         end
 end
%%
%

% for perm=1:6

%
% figure, plot(freq{subject}.freq, squeeze(nanmean(PERM_powspctrm(:,perm,find(elec_post_ind),:))));
% xlabel('frequency (Hz)')
% ylabel('phase-locked amplitude (a.u.)')
% hold on
% for ii=1:2 bar(hz_intermod(ii),1,0.1,'m'); end
% for ii=1:2 bar(Perm(perm,ii),1,0.1,'g'); end
% end
%

% DatSize=numel(noiseIndexes);
% DatIdx=randperm(noiseIndexes);
% noiseIndex(DatIdx(round(rand(1).*10000)));

CONDITION=1
interval=200;
buffer=10;
clear SNR_IMeye SNR_noise
subCoutner=0;
for sub=1:subCounter
    subCoutner=subCoutner+1
    for perm=1:6
        clear hz_intermod
                hz_intermod=intermod([Perm(perm,:)],[1,20]);
        
        hz_intermod(1:3)=Perm(perm,:);
        hz_intermod(4:6)=Perm(perm,:).*2;
        hz_intermod(7:9)=Perm(perm,:).*3;
%             hz_intermod(10:12)=Perm(perm,:).*4;
        
%             hz_intermod(8)=abs(2*Perm(perm,1)-2*Perm(perm,2));
%         % % %     hz_intermod(10)=abs(3*Perm(perm,1)-3*Perm(perm,2));
%             hz_intermod(9)=abs(2*Perm(perm,1)-Perm(perm,2));
%             hz_intermod(10)=abs(Perm(perm,1)-2*Perm(perm,2));
%             hz_intermod(11)=abs(2*Perm(perm,1)+Perm(perm,2));
        % %     hz_intermod(11)=abs(2*Perm(perm,1)-3*Perm(perm,2));
        % %     hz_intermod(11)=abs((Perm(perm,1)+Perm(perm,2))-Perm(perm,3));
        % %     hz_intermod(12)=abs((Perm(perm,1)-Perm(perm,2))+Perm(perm,3));
        %     hz_intermod(12)=abs(Perm(perm,1)+2*Perm(perm,2));
        for intermodul=1:9

            FeatFreq_indx = abs(bsxfun(@minus,hz_intermod(intermodul)',freq{subject}.freq));
            [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
            whichFreqSignal=squeeze((PERM_powspctrm(sub,perm,:,index_freqToFeatures)))';
            %             whichFreqNoise=squeeze((PERM_powspctrm(sub,perm,:,noiseIndexes((round(rand(1).*4000))))))';
            if index_freqToFeatures>interval
                whichNegIdx=(index_freqToFeatures-interval):(index_freqToFeatures-buffer);
            end
            whichPosIdx=(index_freqToFeatures+buffer):(index_freqToFeatures+interval);
            if ismember(whichNegIdx,ALL_GOODindexes)
                whichNegIdx(find(ismember(whichNegIdx,ALL_GOODindexes)))=whichNegIdx(find(ismember(whichNegIdx,ALL_GOODindexes)))-1;
            end
            if ismember(whichPosIdx,ALL_GOODindexes)
                whichPosIdx(find(ismember(whichPosIdx,ALL_GOODindexes)))=whichPosIdx(find(ismember(whichPosIdx,ALL_GOODindexes)))+1;
            end
            whichFreqNoiseNeg=squeeze(nanmean(nanmean(PERM_powspctrm(sub,perm,:,whichNegIdx),2),3))';
            whichFreqNoisePos=squeeze(nanmean(nanmean(PERM_powspctrm(sub,perm,:,whichPosIdx),2),3))';
            
            %specific to each electrode
            MwhichFreqNoiseNeg=nanmean(squeeze(nanmean(PERM_powspctrm(sub,perm,:,whichNegIdx),2)),2)';
            MwhichFreqNoisePos=nanmean(squeeze(nanmean(PERM_powspctrm(sub,perm,:,whichPosIdx),2)),2)';
            STDwhichFreqNoise=nanstd([squeeze(nanmean(PERM_powspctrm(sub,perm,:,whichPosIdx),2)) squeeze(nanmean(PERM_powspctrm(sub,perm,:,whichPosIdx),2))],0,2);
           
            
 
            SNR_temp=whichFreqSignal./(nanmean([MwhichFreqNoiseNeg(:) MwhichFreqNoisePos(:)],2))';
            %
            SNR_tempZ=whichFreqSignal-(nanmean([MwhichFreqNoiseNeg(:) MwhichFreqNoisePos(:)],2))';
            SNR_tempZ=SNR_tempZ./STDwhichFreqNoise(:)';
            
            
            %             corrAmp_band{intermodul}(sub,perm,:)=SNR_tempZ'.*squeeze(Beta{CONDITION}(sub,perm,:));
            
            %             SNR_temp_band{intermodul}(sub,perm,:)=squeeze(nanmean(PERM_powspctrm(sub,perm,:,[whichNegIdx whichPosIdx]),4));
            
            %
            SNR_IMeye{intermodul}(sub,perm,:)=SNR_temp;%SNR_temp;%whichFreqSignal;%SNR_temp;%
            
        end
    end
end

% COND_PERM_powspctrm{CONDITION}=SNR_tempZ;
COND_IM{CONDITION}=SNR_IMeye;%SNR_IMeye;%SNR_noise;%SNR_IMeye
% COND_IM{6}=SNR_IMeye;%SNR_IMeye;%SNR_noise;%SNR_IMeye
% COND_IM{CONDITION_attend}=SNR_IMeye;%SNR_IMeye;%SNR_noise;%SNR_IMeye
%  Alpha_corr{CONDITION_attend}=corrAmp_band;%squeeze(nanmean(PERM_powspctrm(:,:,:,ThoseIndexes),4));
%
%% anova for task x feature x electrodes
% electrodes=find(elec_post_ind);
% load ElectrodesMatrices.mat
% 0z 02 PO8 PO6 P6 P8
% 0z 01 PO7 PO5 P5 P7
% 0z
%load('ElectrodesMatrices.mat');
addpath('/home/adf/faghelss/Documents/TOOLBOXES/miscellaneous/shadedErrorBar/')
electrodesAll=1:64;%find(elec_post_ind)-2;
electrodesLeft=find(elec_postLeft_ind);
electrodesLeft([3 8])=[];
electrodesRight=find(elec_postRight_ind);
electrodesRight([2 8])=[];
electrodesMid=find(elec_postMid_ind);
whichelecs=elecs;

whichMeasure=COND_IM

perm=[1:6];
whichSubs=[1:22]

Taggs=Perm(perm,:)
left_eye=[4 7]
COND_all_mean=zeros(3,length(whichSubs),64);
clear ANOVA_COND ANOVA_COND ANOVA_PERM ANOVA_ELEC  ANOVA_HEMIF

 counter=0;
for non_task=left_eye
    counter=counter+1;
COND_all_mean(1,:,:)=squeeze(COND_all_mean(1,:,:))+(squeeze(nanmean(whichMeasure{1}{non_task}(whichSubs,perm,electrodesAll),2)));
end

right_eye= left_eye+1;
 counter=0;
for non_task=right_eye
    counter=counter+1;
    
COND_all_mean(2,:,:)=squeeze(COND_all_mean(2,:,:))+(squeeze(nanmean(whichMeasure{1}{non_task}(whichSubs,perm,electrodesAll),2)));


    
end

mouth=left_eye+2;
 counter=0;
for non_task=mouth
    counter=counter+1;
COND_all_mean(3,:,:)=squeeze(COND_all_mean(3,:,:))+(squeeze(nanmean(whichMeasure{1}{non_task}(whichSubs,perm,electrodesAll),2)));


   
end
COND_all_mean=COND_all_mean./counter;


clear test1 anovaElec_test1 anovaCond_test1 which anovaFeat_test1

which=elec_postLeft_ind+elec_postMid_ind.*2+elec_postRight_ind.*3;
which(which==0)=[];

test1(:,:,1)=nanmean(COND_all_mean(:,:,which==1),3);
test1(:,:,2)=nanmean(COND_all_mean(:,:,which==2),3);
test1(:,:,3)=nanmean(COND_all_mean(:,:,which==3),3);

anovaElec_test1(:,:,1)=ones(size(test1(:,:,1)));
anovaElec_test1(:,:,2)=ones(size(test1(:,:,1))).*2;
anovaElec_test1(:,:,3)=ones(size(test1(:,:,1))).*3;

anovaCond_test1=ones(size(test1));
temp_sub=[2:2:length(whichSubs)];
anovaCond_test1(:,temp_sub,:)=ones(size(test1(1:3,temp_sub,:)))*2;

anovaFeat_test1(1,:,:)=ones(size(test1(1,:,:))).*1;
anovaFeat_test1(2,:,:)=ones(size(test1(1,:,:))).*2;
anovaFeat_test1(3,:,:)=ones(size(test1(1,:,:))).*3;

[P,ANOVATAB,STATS]=anovan(test1(:),{anovaElec_test1(:),anovaFeat_test1(:),anovaCond_test1(:)},'display','on','model','full','varnames',{'electrodes','features','task'})
% 
temp_sub=[1:2:length(whichSubs)];
ys=[nanmean(nanmean(COND_all_mean(1,temp_sub,whichelecs))) nanmean(nanmean(COND_all_mean(2,temp_sub,whichelecs))) nanmean(nanmean(COND_all_mean(3,temp_sub,whichelecs)))]; errorstd=nanstd([(nanmean(COND_all_mean(1,temp_sub,whichelecs),3)); (nanmean(COND_all_mean(2,temp_sub,whichelecs),3)); (nanmean(COND_all_mean(3,temp_sub,whichelecs),3))],0,2)/sqrt(length(temp_sub));
figure, shadedErrorBar(1:size(ys,2),ys,errorstd,'lineprops',{'k-o','markerfacecolor','k'}),title('Feature signal'),legend({'left-eye signal','right-eye signal','mouth signal'})
temp_sub=[2:2:length(whichSubs)];
ys=[nanmean(nanmean(COND_all_mean(1,temp_sub,whichelecs))) nanmean(nanmean(COND_all_mean(2,temp_sub,whichelecs))) nanmean(nanmean(COND_all_mean(3,temp_sub,whichelecs)))]; errorstd=nanstd([(nanmean(COND_all_mean(1,temp_sub,whichelecs),3)); (nanmean(COND_all_mean(2,temp_sub,whichelecs),3)); (nanmean(COND_all_mean(3,temp_sub,whichelecs),3))],0,2)/sqrt(length(temp_sub));
hold on
shadedErrorBar(1:size(ys,2),ys,errorstd,'lineprops',{'r-o','markerfacecolor','r'}),title(sprintf('Feature signal per task (red:Gender,p_interaction=%.4f)',P(6))),legend({'gender task','emotion task'})%,legend({'left-eye signal','right-eye signal','mouth signal'})
if length(perm)>2
set(gca,'XTick',1:3,'XTickLabel',{'left-eye signal','right-eye signal','mouth signal'})
else
set(gca,'XTick',1:3,'XTickLabel',{num2str(Taggs(1)),num2str(Taggs(2)),num2str(Taggs(3))})
end
figure, 
multcompare(STATS,'dimension',[1])


%%

clear ANOVA_COND COND_all
electrodes_right=find(elec_postRight_ind);
electrodes_Left=find(elec_postLeft_ind);
non_task=1
for perm=6:-1:1

ANOVA_COND(1,:,:,:)=ones(size(COND_IM{2}{non_task}(:,perm,electrodesAll)));
ANOVA_COND(2,:,:,:)=ones(size(COND_IM{2}{non_task}(:,perm,electrodesAll))).*2;

COND_all(1,:,:,:)=COND_IM{1}{non_task}(:,perm,electrodesAll)-COND_IM{1}{non_task+1}(:,perm,electrodesAll);
COND_all(2,:,:,:)=COND_IM{2}{non_task}(:,perm,electrodesAll)-COND_IM{2}{non_task+1}(:,perm,electrodesAll);

P=anovan(COND_all(:),{ANOVA_COND(:)},'display','off')
figure, bar([squeeze(nanmean(nanmean(COND_IM{1}{non_task}(:,perm,electrodesAll)))) squeeze(nanmean(nanmean(COND_IM{2}{non_task}(:,perm,electrodesAll))))]), title (sprintf(' modulation at %1.2f Hz Left and %1.2f Hz Right',Perm(perm,1),Perm(perm,2)))
sigstar({[1,2]},[P])
end
%%
for perm=1:6
non_task=5


COND_all(1,:,:,:)=COND_IM{1}{non_task}(:,perm,electrodesAll);
COND_all(2,:,:,:)=COND_IM{2}{non_task}(:,perm,electrodesAll);


P=anovan(COND_all(:),{ANOVA_COND(:)},'display','off')
figure, bar([squeeze(nanmean(nanmean(COND_IM{2}{non_task}(:,perm,electrodesAll)))) squeeze(nanmean(nanmean(COND_IM{1}{non_task}(:,perm,electrodesAll))))]), title (sprintf(' modulation at %1.2f Hz Left and %1.2f Hz Right',Perm(perm,1),Perm(perm,2)))
sigstar({[1,2]},[P])
end
%% 
posterior_elec={'Oz','O1','O2','POz','PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8','P5','P6','T7','T8','P3','P4','P7','P8'}

addpath /home/adf/faghelss/Documents/TOOLBOXES/fieldtrip-20190419/
ft_defaults
load electrodes_indexes.mat
% load('quickcap64.mat') % load 2D layout of electrodes
% load('ANTneuro64_layout.mat') % load 2D layout of electrodes
% load('ElectrodesMatrices.mat');

elcfile = 'standard_waveguard64.elc';
elec = ft_read_sens(elcfile);
plotcfg.elec=elec;

% 
hdmfile = '/home/adf/faghelss/Documents/TOOLBOXES/fieldtrip-20190419/template/headmodel/standard_bem.mat';
vol = ft_headmodel_asa(hdmfile);
% 
% % the electrodes will be projected on the surface
orig = elec;

% project the electrodes on the volume conduction model
[vol, elec] = ft_prepare_vol_sens(vol, elec);
plotcfg.layout='elec1010.lay';
plotcfg.xlim=[ 2.1 2.15];
plotcfg.colorbar='no';%'East';
plotcfg.style='straight';%'no';%'East';
plotcfg.interpolation='linear';

% cfg           = [];
% cfg.method    = 'interactive';
% cfg.elec      = elec_aligned;
% cfg.headshape = vol.bnd(1);
% elec_aligned  = ft_electroderealign(cfg);
load elec_aligned_head_64ASA_WINK.mat
% save('elec_aligned_head_64ASA_WINK.mat','elec_aligned','cfg','elec','vol')

elec_aligned_show=elec;
% elec1010.lay
% EEG1010.lay

elec_aligned_show.cfg.channel=posterior_elec;
elec_aligned_show.chanpos=elec_aligned_show.chanpos(find(elec_post_ind)-2,:);
% elec_aligned_show.chantype=elec_aligned_show.chantype(find(elec_post_ind)-2);
% elec_aligned_show.chanunit=elec_aligned_show.chanunit(find(elec_post_ind)-2);
elec_aligned_show.elecpos=elec_aligned_show.elecpos(find(elec_post_ind)-2,:);
elec_aligned_show.label=elec_aligned_show.cfg.channel;
elec_aligned_show.cfg.elec=[];




plotcfg.channel=freqcfg.channel;
plotcfg.interpolatenan ='yes';

plotcfg.interplimits = 'electrodes';
plotcfg.interpolation ='linear';%'nearest'
plotcfg.colorbar='no';%'East';
plotcfg.style=[];%'straight';
plotcfg.layout='EEG1010.lay';
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';
whichFreq=[Perm(1,1)];%1.695%8%3.531;%6.312;%1.648%3.531% 6.312 8  %  9.8440 %8%hz_intermod(5);
plotcfg.xlim=[whichFreq-0.00015 whichFreq+0.00015];
load ElectrodesMatrices.mat

%% EEG.chanlocs = readasalocs('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc')
clear vector Features_signal
perm=1:6

for ii=1:3,count{ii}=0;end


for COND_task=1:2
Feats_labels={'','left eye','right eye','mouth'};
figure, 
counter=0;
  
if COND_task==1
    subjects_condition=[1:2:22];
    
%     subjects_condition=[1];
elseif COND_task==2
    subjects_condition=[2:2:22];
%     subjects_condition=[2];
end
% subjects_condition(2)=[]; % need to recompute all subject Laetitia data and trials
%  subjects_condition=1:14;
%  subjects_condition(3:4)=[]; % need to recompute all subject Laetitia data and trials
 

for whichInter=1:3%[1 2 3]+6
  
    
    counter=counter+1;
    Features_signal{COND_task}{whichInter}=squeeze(nanmean(COND_IM{1}{whichInter}(subjects_condition,perm,:,:),2));
    switch mod(whichInter,3)
        case 0
            count{3}=count{3}+1;
            Features_signal{COND_task+2}{3}(count{3},:,:)=squeeze(nanmean(COND_IM{1}{whichInter}(subjects_condition,perm,:,:),2));
        case 1
            count{2}=count{2}+1;
            Features_signal{COND_task+2}{2}(count{2},:,:)=squeeze(nanmean(COND_IM{1}{whichInter}(subjects_condition,perm,:,:),2));
        case 2
            count{1}=count{1}+1;
            Features_signal{COND_task+2}{1}(count{1},:,:)=squeeze(nanmean(COND_IM{1}{whichInter}(subjects_condition,perm,:,:),2));
    end
       vector{COND_task}{whichInter}=(squeeze(nanmean(nanmean(COND_IM{1}{whichInter}(subjects_condition,perm,:,:)),2)));%+squeeze(nanmean(nanmean(COND_IM{COND}{whichInter}(subj_2,perm,:,:)),2)))/sqrt(2);%-squeeze(nanmean((COND_IM{COND}{whichInter+3}(subj,perm,:,:)),2));
       Disvector=vector{COND_task}{whichInter};%-vector{2}{1};
       freq_bog{1}=freq{Subj(1)};
       freq_bog{1}.powspctrm=repmat(Disvector',size(freq{Subj(1)}.powspctrm,1),1,size(freq{Subj(1)}.powspctrm,3));
% % ?       mask=squeeze(vector{COND}{whichInter})>critZ{whichInter}(2);
%        plotcfg.highlightchannel=[];%plotcfg.channel{mask};
%        plotcfg.highlightcolor = [1 0 0];
%        plotcfg.highlight          = 'no';
%        cfg.markersymbol
%        show map

       subplot(1,3,counter)
       ft_topoplotER(plotcfg, freq_bog{1}),colormap(jet)
        set(gca, 'CLim', [1.1 4.5]);
% %         set(gca, 'CLim', [-.5 .5]);
%         set(gca, 'CLim', [3 8]);

% head surface (scalp)
% ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
% hold on;
% 
% ft_plot_vol(vol,'facecolor','cortex')
% hold on
% ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
% % % hold on;
% ft_plot_mesh(vol.bnd(3),'facecolor',[0.4 0.6 0.4]);
% % hold on
% % ft_plot_vol(vol, 'facecolor', 'none')
% hold on
% ft_plot_topo3d(elec_aligned.elecpos(elec_MoneTwo_ind==0,:),freq_bog{1}.powspctrm(elec_MoneTwo_ind==0),'refine',1), title(sprintf(' feature %d',whichInter)),colormap('jet')
% ft_plot_topo3d(plotcfg.elec.elecpos(elec_MoneTwo_ind==0,:),freq_bog{1}.powspctrm(elec_MoneTwo_ind==0)), title(sprintf(' feature %d',whichInter)),colormap(plasma)
%         set(gca,'CLim',[-1 2])
%         ft_plot_topo3d(elec.elecpos(elec_post_ind==1,:),vector{COND}{whichInter}(:)),colorbar, title(sprintf('Condition %d : feature %d',COND,whichInter))
% % colormapeditor
%         rotate3d
% ft_plot_vol(vol, 'facecolor', 'none')
%         ft_plot_sens(elec_aligned_show, 'label', 'on');
% hold on
%         ft_plot_sens(elec_aligned,'label','on');
%         figure,topoplot(Disvector(:),(elec.chanpos(elec_MoneTwo_ind==0,:,:,:)))
%         set(gca, 'CLim', [0 0.1]);
%         subplot(1,2,2),ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),mean(Stats_all{feature}.powspctrm((nSubjects+1):end,elec_MoneTwo_ind==0)),'refine',2),title(sprintf('%d',feature));
%         set(gca, 'CLim', [0 0.1]);
% 
%         figure,ft_topoplotER(plotcfg,  Inter_topoALL{feature})%, title(sprintf('%s',interLabel{feature}));
%     %     hold on
%     %     ft_plot_lay(lay)
%     
%     [maxmax,maxIdx]=max(vector{COND}{whichInter}(:))
end
end


%%
critP=.025
counter=0;
figure,

for whichInter=[7:9]
%     if any(size(Features_signal{4}{whichInter})~=size(Features_signal{3}{whichInter}))
%     Features_signal{4}{whichInter}(1:3,:,:)=[];    
%     end
%     
    counter=counter+1;
    clear tmap2 pmap2 stats temp_mat
    h= waitbar(0, 'ttest 1 : 0 % compl�t�');
    for ii= 1:64
        [~,p{ii},~,stats{ii}]=ttest(Features_signal{1}{whichInter}(:,ii),Features_signal{2}{whichInter}(:,ii));
        tmap2(ii) = stats{ii}.tstat;
        pmap2(ii)= p{ii};
        
        waitbar(ii/6, h, sprintf('ttest : %3.2f %% compl�t�', (ii/64)*100));
    end
    delete(h)
    
%     temp_mat{1}{whichInter}=squeeze(sum(Features_signal{3}{whichInter})/sqrt(count{whichInter}));
%     temp_mat{2}{whichInter}=squeeze(sum(Features_signal{4}{whichInter})/sqrt(count{whichInter}));
%         h= waitbar(0, 'ttest 1 : 0 % compl�t�');
%     for ii= 1:64
%         [~,p{ii},~,stats{ii}]=ttest(temp_mat{1}{whichInter}(:,ii),temp_mat{2}{whichInter}(:,ii));
%         tmap2(ii) = stats{ii}.tstat;
%         pmap2(ii)= p{ii};
%         
%         waitbar(ii/6, h, sprintf('ttest : %3.2f %% compl�t�', (ii/64)*100));
%     end
%     delete(h)
%     
    % figure, imagesc(tmap2),colorbar
    % set(gca,'XTick',1:64,'XTickLabel',freq{2}.label)
     if any(pmap2<critP)
       mask=pmap2<critP;
       plotcfg.highlightchannel=[];%plotcfg.channel{mask};
       plotcfg.highlightcolor = [0 0 0];
       plotcfg.highlight          = 'yes';
       plotcfg.highlight=[]
%        cfg.markersymbol
     end
    plotcfg.channel{pmap2<critP};
    vector{COND_task}{whichInter}=tmap2(:);
    Disvector=vector{COND_task}{whichInter};
    freq_bog{1}=freq{2};
    freq_bog{1}.powspctrm=repmat(Disvector',size(freq{2}.powspctrm,1),1,size(freq{2}.powspctrm,3));
    subplot(2,4,counter),title(Feats_labels{mod(whichInter,3)+1})
    ft_topoplotER(plotcfg, freq_bog{1}),colormap(jet),title(sprintf('nb sign %d',sum(pmap2<.01)))
            set(gca, 'CLim', [-2.5 3.5]);
end
%% % find(elec_post_ind)
% DATA1=zeros(size(COND_IM{1}{intermodul}(:,[28:30 61 62])));
% DATA2=zeros(size(COND_IM{1}{intermodul}(:,[28:30 61 62])));
% DATA3=zeros(size(COND_IM{1}{intermodul}(:,[28:30 61 62])));

% DATA1=zeros(size(COND_IM{1}{intermodul}(:,find(elec_post_ind))));
% DATA2=zeros(size(COND_IM{1}{intermodul}(:,find(elec_post_ind))));
% DATA3=zeros(size(COND_IM{1}{intermodul}(:,find(elec_post_ind))));
counterInter=0;
counterCond=1;
clear IM_all ANOVA_COND ANOVA_PERM IM_allINTERMOD ANOVA_COND_INTERMD ANOVA_PERM_INTERMD ANOVA_INTERMOD
for intermodul=[2]
    counterInter=counterInter+1;
    %      figure, bar([nanmean(COND_IM{1}{intermodul}(:,find(elec_post_ind)),2)]-[nanmean(COND_IM{2}{intermodul}(:,find(elec_post_ind)),3)])
    %      hold on
    %      bar([nanmean(COND_IM{5}{intermodul}(:,find(elec_post_ind)),2)])
    %
    %
    %         for ii=4:5
    %         permIdx{ii}=[find(nanmean(squeeze(nanmean(COND_IM{ii}{intermodul}(:,find(elec_post_ind)),3))>1))];
    %         end
    for perm=1:6
        
        
        
        % figure, bar([nanmean(COND_IM{4}{intermodul}(:,find(elec_post_ind)),2)])
        % figure, bar([nanmean(COND_IM{2}{intermodul}(:,find(elec_post_ind)),2)])
        % figure, bar([nanmean(COND_IM{3}{intermodul}(:,find(elec_post_ind)),2)])%+nanmean(SNR_IMeye{2}(:,find(elec_post_ind)),2)]/2)
        % figure, bar([nanmean(SNR_IMeye{intermodul}(:,31),2)])%+nanmean(SNR_IMeye{2}(:,find(elec_post_ind)),2)]/2)
        % figure, bar([nanmean(SNR_IMeye{1}(:,29:31),2)+nanmean(SNR_IMeye{2}(:,29:31),2)]/2
        %figure, bar(mean(SNR_IMeye{1}(:,18),2))
        %
        % DATA1=COND_IM{1}{intermodul}(permIdx{1},find(elec_post_ind));%+COND_IM{1}{2}(permIdx{1},find(elec_post_ind));
        % DATA2=COND_IM{2}{intermodul}(permIdx{2},find(elec_post_ind));%+COND_IM{2}{2}(permIdx{2},find(elec_post_ind));
        % DATA3=COND_IM{3}{intermodul}(permIdx{3},find(elec_post_ind));%+COND_IM{3}{2}(permIdx{3},find(elec_post_ind));
        
        %
        % DATA1=COND_IM{1}{intermodul}(permIdx{1},[29 61 62])+COND_IM{1}{2}(permIdx{1},[29 61 62]);
        % DATA2=COND_IM{2}{intermodul}(permIdx{2},[29 61 62])+COND_IM{2}{2}(permIdx{2},[29 61 62]);
        % DATA3=COND_IM{3}{intermodul}(permIdx{3},[29 61 62])+COND_IM{3}{2}(permIdx{3},[29 61 62]);
        
        % DATA1=DATA1+COND_IM{1}{intermodul}(:,[28:30 61 62]);
        % DATA2=DATA2+COND_IM{2}{intermodul}(:,[28:30 61 62]);
        % DATA3=DATA3+COND_IM{3}{intermodul}(:,[28:30 61 62]);
        
        % DATA1=DATA1+COND_IM{1}{intermodul}(:,find(elec_post_ind));
        % DATA2=DATA2+COND_IM{2}{intermodul}(:,find(elec_post_ind));
        % DATA3=DATA3+COND_IM{3}{intermodul}(:,find(elec_post_ind));
        
        
        for ii=1:2
            
            ANOVA_COND(ii,:,perm,:)=ones(size(COND_IM{ii}{intermodul}(:,perm,find(elec_post_ind)))).*ii;
            IM_all(ii,:,perm,:)=((COND_IM{ii}{intermodul}(:,perm,find(elec_post_ind))));%+(COND_IM{ii}{2}(:,perm,find(elec_post_ind))))/2;
            %             IM_all(ii-3,perm,:)=((COND_IM{ii}{intermodul}(perm,find(elec_post_ind)))+(COND_IM{ii}{2}(perm,find(elec_post_ind))))/2;
            IM_allINTERMOD(counterInter,ii,:,perm,:)=(COND_IM{ii}{intermodul}(:,perm,find(elec_post_ind)));%+(COND_IM{ii}{2}(perm,[28:30 61 62])))/2;
            if perm==6 && ii==4
                counterCond=counterCond+1
            end
        end
        %
        ANOVA_PERM(:,:,perm,:)=ones(2,nbSub,1,size(COND_IM{ii}{intermodul}(:,perm,find(elec_post_ind)),3)).*perm;
        
        
    end
    
    ANOVA_INTERMOD(counterInter,:,:,:,:)=ones(size(IM_allINTERMOD(counterInter,:,:,:,:))).*counterInter;
end
%

ANOVA_COND_INTERMD=repmat(ANOVA_COND,counterInter,1,1,1);
ANOVA_PERM_INTERMD=repmat(ANOVA_PERM,counterInter,1,1,1);
IM_all(ANOVA_COND(:)==0)=[];
ANOVA_COND(ANOVA_COND(:)==0)=[];
% [P,ANOVATAB,STATS]=anovan(IM_all(:),{ANOVA_COND(:), ANOVA_PERM(:)},'model','full','varnames',{'CONDITION','PERMUTATION'})
% [P,ANOVATAB,STATS]=anovan(IM_allINTERMOD(:),{ANOVA_COND_INTERMD(:), ANOVA_PERM_INTERMD(:)},'model','full','varnames',{'CONDITION','PERMUTATION'})
[P,ANOVATAB,STATS]=anova1(IM_all(:),ANOVA_COND(:))
figure,
multcompare(STATS,'dimension',[1])
%%
[h,p1]=ttest2(DATA1(:),DATA3(:))
[h,p2]=ttest2(DATA1(:),DATA2(:))
figure, bar([mean(DATA1(:)) mean(DATA3(:)) mean(DATA2(:))])
sigstar({[1,2],[1,3]},[p1 p2])




% end
%%

CONDITION=3;
switch (CONDITION)
    case 1
        Subj=Subj_cond1((end-7):end);
    case 2
        Subj=Subj_cond2((end-7):end);
    case 3
        Subj=Subj_cond3;
        
end

for subject=Subj
    for perm=1:6
        if any(PERM{subject}{perm})
            for trial=PERM{subject}{perm}
                
            end
        end
    end
end



%% INTERMODULATION FEATURES : SNR
interLabel={'Left && Right', 'Left && Mouth', 'Right && Mouth','All','Blink non-linearity','AVG Mouth-eye pair'};
load('quickcap64.mat') % load 2D layout of electrodes
electrodes = ft_read_sens('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc');
freqcfg.elec=electrodes;
cfg.elecfile=('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc');


Subj_cond1=[ 8 12 14 16 18 20 22 24 26 28 30 33 35 38 41 43 45 48]

Subj_cond2=[ 9 13 17 19 21 23 25 27 29 31 36 39 46 49]

Subj_cond3=[ 32 34 37 40 42 44 47 50]

% all 3 conditions;
Subj_cond1_all=[30   35    38    41    43    45     48];

Subj_cond2_all=[31   36    39    13    15    46    49];

Subj_cond3_all=[32   37    40    42    44    47    50];


% Subj_cond1_all=[30   35    38       45     48];
%
% Subj_cond2_all=[31   36    39       46    49];
%
% Subj_cond3_all=[32   37    40       47    50];

clear Stats_all_band InterFeature_topo_wtvFeat_all Stats_all
%%
% clear freq
CONDITION=2;
%  counter=0;
switch (CONDITION)
    case 1
        %      load SNR_subjectxtrialCond1.mat
        load SNR_subjectxtrialCond1_POSTICA.mat
        %     Subj_cond1=[2 3 4 8 12 14 16 18 20 22 24]
        %     Subj_cond1=[ 12 14 16 18 20 22 24 26 28] % Cond 1 only left out 14
        Subj_cond1=[ 8 12 14 16 18 20 22 24 26 28 30 33 35 38 41 43 45 48]
        
        %     Subj_cond1=[8 12 16 20 22 24 26 28] % 14 not sure of tag, 18 has very low SNR
        %     Subj_cond1=[24] % NEED TO CHECK CHLOE's (20) comp data.
        Subj=Subj_cond1(5:end);
        % Subj=Subj_cond1;
        % Subj=Subj_cond1_all;
    case 2
        %     load SNR_subjectxtrialCond2_POSTICA.mat
        load SNR_subjectxtrialCond2_POSTICA_2.mat
        %     load SNR_subjectxtrialCond2_POSTICA.mat
        %     Subj_cond2=[ 13 15 17 19 21 23 25 27 29] % Cond 2 only left out 15 PROBLEM WITH 27 and 29 subject
        Subj_cond2=[ 9 13 17 19 21 23 25 27 29 31 36 39 46 49]
        %    Subj_cond2=[ 9 13 17 21 23 25 27 29]
        %     Subj_cond2=25; % FSAVOI (19) weird, CATHERINE (9) same pattern as COND1 (8), 15 a bit weird
        Subj=Subj_cond2;%((end-7):end);
        % Subj=Subj_cond2_all;
    case 3
        load SNR_subjectxtrialCond3.mat
        %         Subj_cond3(3)=[];
        %         Subj_cond3(3)=[];
        Subj=Subj_cond3;
        %         Subj=Subj_cond3_all;
        
end
%
BOGUSFREQ=Subj(1);
%
nSubjects=length(Subj);

%


%
rng('default')
clear plotcfg
plotcfg.channel= freq{BOGUSFREQ}.label;
plotcfg.colorbar='East';
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='East';
% load('ANTneuro64_layout.mat') % load 2D layout of electrodes
% plotcfg.elec=electrodes;
% layout=ft_prepare_layout(cfg,freq{BOGUSFREQ})
layout.outline=layout.mask;
plotcfg.layout=lay;
plotcfg.showlabels ='yes';


clear FreqDescriptive

for ii=1:12
    
    InterFeature_topoALL{ii}=AllSub_Features{1};
    InterFeature_topoALL{ii}.label=freq{BOGUSFREQ}.label;
    InterFeature_topoALL{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
    Inter_topoALL{ii}=AllSub_Features{1};
    Inter_topoALL{ii}.label=freq{BOGUSFREQ}.label;
    Inter_topoALL{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
end

for subject=Subj
    
    
    counter=counter+1
    %     subject
    clear data
    
    for ii=1:12
        Inter_topo{subject}{ii}=freq{BOGUSFREQ};
        Inter_topo{subject}{ii}.powspctrm=zeros(1,size(freq{BOGUSFREQ}.powspctrm,2));
        Inter_topo{subject}{ii}.freq=1;
        Inter_difference{subject}{ii}=freq{BOGUSFREQ};
        Inter_difference{subject}{ii}.powspctrm=zeros(size(freq{BOGUSFREQ}.powspctrm,1),size(freq{BOGUSFREQ}.powspctrm,2));
        %         Inter_difference{subject}{ii}.freq=1;
        InterFeature_topo{subject}{ii}=AllSub_Features{1};
        InterFeature_topo{subject}{ii}.label=freq{BOGUSFREQ}.label;
        InterFeature_topo{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
        InterFeature_topo_wtvFeat{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
        Freq_band{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
        Noise{subject}{ii}.powspctrm=zeros(1,length(freq{BOGUSFREQ}.label));
    end
    
    
    
    blockCounter=1;
    trialBlckCounter=0;
    
    nTrials=size(SNR_Elec_powSpectrum{subject},1)
    if subject==39
        nTrials=12;
    end
    for trial=1:nTrials
        
        if subject==22
            eval(sprintf('load(''%s_blck%d_SSVEP_0_1_1.mat'') ',Subjects{subject}{3} ,Subjects{subject}{4}(blockCounter)))
        else
            %             if CONDITION==1
            eval(sprintf('load(''%s_blck%d_SSVEP_0_1_%d.mat'') ',Subjects{subject}{3} ,blockCounter,CONDITION))
            %             else
            %                 eval(sprintf('load(''%s_blck%d_SSVEP_0_1_2.mat'') ',Subjects{subject}{3} ,blockCounter))
            %             end
        end
        
        if  trial== 5 || trial== 9 || trial== 13 || trial== 17
            blockCounter=blockCounter+1;
            
            if subject==22 && CONDITION==1
                eval(sprintf('load(''%s_blck%d_SSVEP_0_1_1.mat'') ',Subjects{subject}{3} ,Subjects{subject}{4}(blockCounter)))
            else
                eval(sprintf('load(''%s_blck%d_SSVEP_0_1_%d.mat'') ',Subjects{subject}{3} ,blockCounter,CONDITION))
            end
            trialBlckCounter=1;
        end
        
        
        if mod(trial,4)~=0
            trial_feature=cell2mat(data{mod(trial,4)}.freq);
            harmon_feat=cell2mat(data{mod(trial,4)}.freq)*2;
        else
            trial_feature=cell2mat(data{4}.freq);
            harmon_feat=cell2mat(data{4}.freq)*2;
        end
        hz_intermod=(intermod([trial_feature]',[1 20]));
        
        hz_intermod(8)=abs(trial_feature(1)-2*trial_feature(2));
        hz_intermod(9)=abs(2*trial_feature(1)-trial_feature(2));
        hz_intermod(10)=abs(2*trial_feature(1)-2*trial_feature(2));
        %         hz_intermod(11)=abs(2*trial_feature(1)+2*trial_feature(2));
        hz_intermod_harmon=(intermod([harmon_feat]',[1 20]));
        
        hz_inter_wtvFeat=intermod([ 6.3150    8.0000    3.5290]',[0.5 20]);
        hz_inter_wtvFeat=sort(hz_inter_wtvFeat(1:8));
        for whichFreq=1:8
            
            FeatFreq_indx = abs(bsxfun(@minus,hz_inter_wtvFeat(whichFreq)',freq{BOGUSFREQ}.freq));
            [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
            InterFeature_topo_wtvFeat{subject}{whichFreq}.powspctrm(:)=InterFeature_topo_wtvFeat{subject}{whichFreq}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
            %             InterFeature_topo_wtvFeat{subject}{whichFreq}.powspctrm(:)=InterFeature_topo_wtvFeat{subject}{whichFreq}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,noiseIndexes(round(rand(1,1).*2400))))';
        end
        for intermodul=1:length(hz_intermod)
            
            if hz_intermod(intermodul)>1
                % finding fundamental frequency attached to trial feature and adding it to feature
                % topographical matrix
                FeatFreq_indx = abs(bsxfun(@minus,hz_intermod(intermodul)',freq{BOGUSFREQ}.freq));
                [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
                all_good_indexesINTER(intermodul)=index_freqToFeatures;
                
                whichFreqSignal=squeeze(freq{subject}.powspctrm(trial,:,index_freqToFeatures))';
                whichFreqNoiseNeg=(squeeze(mean(freq{subject}.powspctrm(trial,:,(index_freqToFeatures-8):(index_freqToFeatures-1)),3)))';
                whichFreqNoisePos=(squeeze(mean(freq{subject}.powspctrm(trial,:,(index_freqToFeatures+1):(index_freqToFeatures+8)),3)))';
                Noise{subject}{intermodul}.powspctrm= Noise{subject}{intermodul}.powspctrm+((whichFreqNoiseNeg+whichFreqNoisePos)/2)';
                SNR_temp=whichFreqSignal'./((whichFreqNoiseNeg+whichFreqNoisePos)/2)';
                %             InterFeature_topo{subject}{intermodul}.powspctrm(:)=InterFeature_topo{subject}{intermodul}.powspctrm(:)+squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
                %             InterFeature_topo{subject}{intermodul}.powspctrm(:)=InterFeature_topo{subject}{intermodul}.powspctrm(:)+squeeze(freq{subject}.powspctrm(trial,:,index_freqToFeatures))';
                InterFeature_topo{subject}{intermodul}.powspctrm(:)=InterFeature_topo{subject}{intermodul}.powspctrm(:)+SNR_temp';
                Feature_pow{intermodul}(counter,trial,:)=SNR_temp';%squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
                Inter_difference{subject}{intermodul}.powspctrm(trial,:)=squeeze(SNR_Elec_powSpectrum{subject}(trial,:,index_freqToFeatures))';
            else
                continue
            end
        end
        
        Band{1}=[2.5 4.5];
        Band{2}=[9.5 14];
        for whichBand=1:2
            
            idx = abs(bsxfun(@minus,Band{whichBand}(1),freq{BOGUSFREQ}.freq));
            [~,Band_indx{1}] = min(idx(:,1:size(idx,2)));
            idx = abs(bsxfun(@minus,Band{whichBand}(2),freq{BOGUSFREQ}.freq));
            [~,Band_indx{2}] = min(idx(:,1:size(idx,2)));
            Freq_band{subject}{whichBand}.powspctrm(:)=Freq_band{subject}{whichBand}.powspctrm(:)+squeeze(mean(freq{subject}.powspctrm(trial,:,Band_indx{1}:Band_indx{2}),3))';
        end
    end
    
    
    Inter_topo{subject}{1}.powspctrm(:)=(InterFeature_topo{subject}{1}.powspctrm(:));%+InterFeature_topo{subject}{2}.powspctrm(:))/(2);
    Inter_topo{subject}{2}.powspctrm(:)=(InterFeature_topo{subject}{3}.powspctrm(:));%+InterFeature_topo{subject}{4}.powspctrm(:))/(2);
    
    Inter_topo{subject}{3}.powspctrm(:)=(InterFeature_topo{subject}{5}.powspctrm(:));%+InterFeature_topo{subject}{6}.powspctrm(:))/(2);
    Inter_topo{subject}{4}.powspctrm(:)= InterFeature_topo{subject}{2}.powspctrm(:);%+InterFeature_topo{subject}{6}.powspctrm(:))/(2*size(freq{subject}{1}.powspctrm,1));
    
    Inter_topo{subject}{5}.powspctrm(:)=(InterFeature_topo{subject}{4}.powspctrm(:));%-(Inter_topo{subject}{2}.powspctrm(:)+Inter_topo{subject}{3}.powspctrm(:))/2;
    Inter_topo{subject}{6}.powspctrm(:)=(InterFeature_topo{subject}{6}.powspctrm(:));%+Inter_topo{subject}{3}.powspctrm(:))/2;
    
    Inter_topo{subject}{7}.powspctrm(:)=((Inter_topo{subject}{7}.powspctrm(:)));%+Inter_topo{subject}{4}.powspctrm(:))/2);
    Inter_topo{subject}{8}.powspctrm(:)=(((InterFeature_topo{subject}{8}.powspctrm(:))));%+Inter_topo{subject}{4}.powspctrm(:))/2);
    Inter_topo{subject}{9}.powspctrm(:)=(((InterFeature_topo{subject}{9}.powspctrm(:))));%+Inter_topo{subject}{4}.powspctrm(:))/2);
    Inter_topo{subject}{10}.powspctrm(:)=(((InterFeature_topo{subject}{10}.powspctrm(:))));%+Inter_topo{subject}{4}.powspctrm(:))/2);
    %     Inter_topo{subject}{11}.powspctrm(:)=(((InterFeature_topo{subject}{11}.powspctrm(:))));%+Inter_topo{subject}{4}.powspctrm(:))/2);
    
    
    % Inter_topoALL{5}.powspctrm(:)=Inter_topoALL{5}.powspctrm(:)+Inter_topo{subject}{5}.powspctrm(:);
    % Inter_topoALL{6}.powspctrm(:)=Inter_topoALL{6}.powspctrm(:)+Inter_topo{subject}{6}.powspctrm(:);
    
    for intermodul=1:length(hz_intermod),
        
        InterFeature_topo{subject}{intermodul}.powspctrm(:)= InterFeature_topo{subject}{intermodul}.powspctrm(:)/nTrials;
        Noise{subject}{intermodul}.powspctrm(:)=Noise{subject}{intermodul}.powspctrm(:)/nTrials;
        
    end
    for whichFreq=1:8
        InterFeature_topo_wtvFeat_all{whichFreq}.powspctrm(counter,:)=InterFeature_topo_wtvFeat{subject}{whichFreq}.powspctrm(:)/nTrials;
    end
    for ii=1:10,
        
        Stats_all{ii}.powspctrm(counter,:)=InterFeature_topo{subject}{ii}.powspctrm(:);
        %         Stats_all{ii}.powspctrm(counter,:)=Noise{subject}{ii}.powspctrm(:);
        Noise_all{ii}.powspctrm(counter,:)=Noise{subject}{ii}.powspctrm(:);
        
    end%Inter_topo{subject}{ii}.powspctrm(:)/nTrials; end
    
    
    for whichBand=1:2
        Stats_all_band{whichBand}.powspctrm(counter,:)=Freq_band{subject}{whichBand}.powspctrm(:)/nTrials;
    end
    
    % Stats_all{7}.powspctrm(counter,:)=Inter_topo{subject}{1}.powspctrm(:)-(Inter_topo{subject}{3}.powspctrm(:)+Inter_topo{subject}{5}.powspctrm(:))/2;
    % Stats_all{8}.powspctrm(counter,:)=Inter_topo{subject}{2}.powspctrm(:)-(Inter_topo{subject}{4}.powspctrm(:)+Inter_topo{subject}{6}.powspctrm(:))/2;
    % Stats_all{9}.powspctrm(counter,:)=(Stats_all{7}.powspctrm(:)+Stats_all{8}.powspctrm(:))/2;
    % figure, ft_topoplotER(plotcfg,Inter_topo{subject}{7})
    %
    % for intermodul=1:6,FreqDescriptive{subject}{intermodul}= ft_freqdescriptives(cfg, Inter_topo{subject}{intermodul});  end
    %
    % Inter_topo{subject}{5}.powspctrm(:)=Inter_topo{subject}{5}.powspctrm(:)/inter_Feature;
end
for subject=Subj
    for intermodul=1:length(hz_intermod),
        
        
        InterFeature_topoALL{intermodul}.powspctrm(:)=InterFeature_topoALL{intermodul}.powspctrm(:)+InterFeature_topo{subject}{intermodul}.powspctrm(:);
        
    end
end

for intermodul=1:length(hz_intermod),
    
    InterFeature_topoALL{intermodul}.powspctrm(:)=InterFeature_topoALL{intermodul}.powspctrm(:)/nSubjects;
end

%%
plotcfg.zlim=[]
for feature=[1];
    
    bogusTopo=Feature_topoALL{feature};
    bogusTopo.powspctrm=TESTANOV(1:11,:)-TESTANOV(12:22,:);%(Stats_all{feature}.powspctrm(1:7,:)+Stats_all{4}.powspctrm(1:7,:))-(Stats_all{feature}.powspctrm(15:21,:)+Stats_all{4}.powspctrm(15:21,:));%Stats_all_band{2}.powspctrm(1:8,:)-Stats_all_band{2}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(1:8,:)-Stats_all{feature}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(17:24,:);%Stats_all{feature}.powspctrm(9:16,:);
    bogusTopo.powspctrm=TESTANOV(23:33,:)-TESTANOV(34:44,:);%
    %  bogusTopo.powspctrm=(Stats_all{feature}.powspctrm(1:7,:))-((Stats_all{feature}.powspctrm(8:14,:)))%+Stats_all{feature}.powspctrm(15:end,:))/2);%Stats_all_band{2}.powspctrm(1:8,:)-Stats_all_band{2}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(1:8,:)-Stats_all{feature}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(17:24,:);%Stats_all{feature}.powspctrm(9:16,:);
    % freq{1}.powspctrm=repmat(mean(Stats_all{feature}.powspctrm(1:9,:)),15,1,2433);
    %     freq{2}.powspctrm=repmat(mean(Stats_all{feature}.powspctrm(10:end,:)),15,1,2433);
    %
    %     figure, ft_topoplotER(plotcfg,Inter_topoALL{1})
    %     figure, ft_topoplotER(plotcfg,InterFeature_topoALL{feature})
    figure, ft_topoplotER(plotcfg,bogusTopo)
    %
end
%%
%   clear StatInterFeatureAllCond Stats_all design
% counter=0;
condition=2;
SubCounter=0;
for ii=1:6,  Stats_all=Inter_topo{subject}{5}; StatInterFeatureAllCond{ii}=Stats_all; end
for subjects=Subj
    SubCounter=SubCounter+1;
    %     for trial=1:size(Inter_topo{subjects}{6}.powspctrm(:,:),1)
    counter=counter+1
    for ii=5
        Stats_all.powspctrm(counter,:)=Inter_topo{subjects}{ii}.powspctrm(:);
        %
        design(1,counter)= SubCounter;
        design(2,counter)= condition;
        StatInterFeatureAllCond{ii}.powspctrm(counter,:)=Inter_topo{subject}{ii}.powspctrm(:);
    end
    %     end
end



%%

cfg.variance      = 'yes'
cfg.jackknife     = 'no'
cfg.keeptrials    = 'yes' % or 'no', estimate single trial power (useful for fourier data) (default = 'no')
% cfg.channel       =
cfg.trials        = 'all' %or a selection given as a 1xN vector (default = 'all')
% cfg.frequency     = %[fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
% cfg.latency       = %[tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
cfg=[];

% [COND1grandavg] = ft_freqgrandaverage(cfg,FreqDescriptive{8}{5},FreqDescriptive{12}{5},FreqDescriptive{14}{5},FreqDescriptive{16}{5},FreqDescriptive{18}{5},FreqDescriptive{22}{5},FreqDescriptive{24}{5})
%
%  [COND1grandavg] = ft_freqgrandaverage(cfg,Inter_topo{8}{5})

Inter_topoALL{1}.powspctrm(:)=(InterFeature_topoALL{1}.powspctrm(:)+InterFeature_topoALL{2}.powspctrm(:))/(2);

% Inter_topoALL{2}.powspctrm(:)=(InterFeature_topoALL{3}.powspctrm(:)+InterFeature_topoALL{4}.powspctrm(:))/(nSubjects*2);
%
% Inter_topoALL{3}.powspctrm(:)=(InterFeature_topoALL{5}.powspctrm(:)+InterFeature_topoALL{6}.powspctrm(:))/(nSubjects*2);
Inter_topoALL{4}.powspctrm(:)=InterFeature_topoALL{7}.powspctrm(:)/(nSubjects);



% Differences between eye integration and mouth eye integration
Inter_topoALL{5}.powspctrm(:)=Inter_topoALL{1}.powspctrm(:)-(Inter_topoALL{2}.powspctrm(:)+Inter_topoALL{3}.powspctrm(:))/2;
Inter_topoALL{6}.powspctrm(:)=(Inter_topoALL{2}.powspctrm(:)+Inter_topoALL{3}.powspctrm(:))/2;

Inter_topoALL{8}.powspctrm=(InterFeature_topoALL{1}.powspctrm(:)-(InterFeature_topoALL{3}.powspctrm(:)+InterFeature_topoALL{5}.powspctrm(:))/2)';
Inter_topoALL{9}.powspctrm=(InterFeature_topoALL{2}.powspctrm(:)-(InterFeature_topoALL{4}.powspctrm(:)+InterFeature_topoALL{6}.powspctrm(:))/2)';
Inter_topoALL{10}.powspctrm=(Inter_topoALL{9}.powspctrm+Inter_topoALL{8}.powspctrm)/2;
% Inter_topoALL{5}.powspctrm(:)=Inter_topoALL{5}.powspctrm(:)/(nSubjects);
% Inter_topoALL{6}.powspctrm(:)=Inter_topoALL{6}.powspctrm(:)/(nSubjects);

%%
% cfg.parameter        = 'powspctrm';%Connect_index_cell{whichConnect};  %'
% % cfg.zlim=[0 1]
plotcfg.highlightchannel=posterior_elec_Left;%{'PO5','P5','TP7','PO7','PO3','T7','O1'};
plotcfg.highlight ='label';
plotcfg.highlightsymbol    = 'o';
plotcfg.highlightcolor= [1 0 0];
plotcfg.highlightsize      = 10;%highlight marker size (default = 6)
plotcfg.highlightfontsize  = 12;%highlight marker size (default = 8)

plotcfg.highlightchannel=posterior_elec;%{'PO6','P6','TP8','PO6','PO4','T8','O2','PO8'};
plotcfg.highlight ='label';
plotcfg.highlightsymbol    = '*';
plotcfg.highlightcolor= [0 0 0];
plotcfg.highlightsize      = 10;%highlight marker size (default = 6)
plotcfg.highlightfontsize  = 12;%highlight marker size (default = 8)

% posterior_elec= {'Oz','Cz','Pz','CPz','O1','O2','POz','PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8','T7','T8'};%'T7','T8''O1', 'O2''Pz','CPz'
% plotcfg.layout = 'butterfly';
plotcfg.channel= 'all';%posterior_elec;
plotcfg.interpolation ='v4'%'nearest'
plotcfg.colormap='jet'
% plotcfg.zlim=[1 1.8]
plotcfg.interplimits='head';
freq{1}=freq{BOGUSFREQ};
freq{2}=freq{1};
plotcfg.style  =[];%'straight'
plotcfg.layout.label(31)={'Oz'}
plotcfg.layout.label(29)={'POz'}
posidx = [1 2 3 33 34 35 36 4 37 5 38 6 39 7 40 8 59 9 41 10 42 11 43 12 60 14 44 15 45 16 46 17 47 18 61 20 48 21 49 22 50 23 62 24 51 25 52 26 53 27 54 28 63 55 56 29 57 58 64 30 31 32 13 19];
for feature=2
    
    bogusTopo=InterFeature_topoALL{feature};
    bogusTopo.powspctrm=(Stats_all{feature}.powspctrm(1:14,:)+Stats_all{2}.powspctrm(1:14,:))-(Stats_all{feature}.powspctrm(15:28,:)+Stats_all{2}.powspctrm(15:28,:));%-(Stats_all{feature}.powspctrm(15:21,:)+Stats_all{4}.powspctrm(15:21,:));%Stats_all_band{2}.powspctrm(1:8,:)-Stats_all_band{2}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(1:8,:)-Stats_all{feature}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(17:24,:);%Stats_all{feature}.powspctrm(9:16,:);
    % bogusTopo.powspctrm=(Stats_all{feature}.powspctrm(9:16,:)+Stats_all{2}.powspctrm(9:16,:))-(Stats_all{feature}.powspctrm(17:24,:)+Stats_all{2}.powspctrm(17:24,:));%-(Stats_all{feature}.powspctrm(15:21,:)+Stats_all{4}.powspctrm(15:21,:));%Stats_all_band{2}.powspctrm(1:8,:)-Stats_all_band{2}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(1:8,:)-Stats_all{feature}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(17:24,:);%Stats_all{feature}.powspctrm(9:16,:);
    % bogusTopo.powspctrm=(Stats_all{feature}.powspctrm(1:8,:)+Stats_all{2}.powspctrm(1:8,:))-(Stats_all{feature}.powspctrm(9:16,:)+Stats_all{2}.powspctrm(9:16,:));%-(Stats_all{feature}.powspctrm(15:21,:)+Stats_all{4}.powspctrm(15:21,:));%Stats_all_band{2}.powspctrm(1:8,:)-Stats_all_band{2}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(1:8,:)-Stats_all{feature}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(17:24,:);%Stats_all{feature}.powspctrm(9:16,:);
    
    %  bogusTopo.powspctrm=(Stats_all{feature}.powspctrm(1:7,:))-((Stats_all{feature}.powspctrm(8:14,:)))%+Stats_all{feature}.powspctrm(15:end,:))/2);%Stats_all_band{2}.powspctrm(1:8,:)-Stats_all_band{2}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(1:8,:)-Stats_all{feature}.powspctrm(9:16,:);%Stats_all{feature}.powspctrm(17:24,:);%Stats_all{feature}.powspctrm(9:16,:);
    %  bogusTopo.powspctrm=(Stats_all{feature}.powspctrm(1:8,:))-((Stats_all{feature}.powspctrm(9:16,:)));
    
    
    % Faut pas mettre 65 et 66 qui sont � COMNT � et � SCALE �. Voici un exemple, avec la fct de bas niveau.
    
    % figure, ft_plot_topo(layout.pos(1:62,1), layout.pos(1:62,2), bogusTopo.powspctrm(posidx), 'mask', layout.mask, 'outline', layout.outline, 'interplim', 'mask', 'interpmethod','v4');
    
    % freq{1}.powspctrm=repmat(mean(Stats_all{feature}.powspctrm(1:9,:)),15,1,2433);
    %     freq{2}.powspctrm=repmat(mean(Stats_all{feature}.powspctrm(10:end,:)),15,1,2433);
    %
    %     figure, ft_topoplotER(plotcfg,Inter_topoALL{1})
    %     figure, ft_topoplotER(plotcfg,InterFeature_topoALL{feature})
    figure, ft_topoplotER(plotcfg,bogusTopo)
    %
end
%%
plotcfg.colorbar='South';
plotcfg.zlim=[];
plotcfg.elec=electrodes;
plotcfg.lay=load('quickcap64.mat') % load 2D layout of electrodes
DATA.elec=[];%electrodes;
DATA.lay=plotcfg.lay;
figure, ft_topoplotER(plotcfg,freq{1})

%%
for feature=1
    if feature==5
        figure,
        subplot(1,2,1),ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),mean(Stats_all{feature}.powspctrm(1:9,:)),'refine',2),title(sprintf('%d',feature));
        set(gca, 'CLim', [0 0.1]);
        subplot(1,2,2),ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),mean(Stats_all{feature}.powspctrm((nSubjects+1):end,:)),'refine',2),title(sprintf('%d',feature));
        set(gca, 'CLim', [0 0.1]);
    else
        figure,
        subplot(1,2,1),ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),mean(Stats_all{feature}.powspctrm(1:9,:)),'refine',2),title(sprintf('%d',feature))
        set(gca, 'CLim', [1 1.15]);
        subplot(1,2,2),ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),mean(Stats_all{feature}.powspctrm((nSubjects+1):end,:)),'refine',2),title(sprintf('%d',feature))
        set(gca, 'CLim', [1 1.15]);
    end
    %     figure,ft_topoplotER(plotcfg,  Inter_topoALL{feature})%, title(sprintf('%s',interLabel{feature}));
    % %     hold on
    % %     ft_plot_lay(lay)
end
%
% InterCOND1=InterFeature_topoALL;

%%

figure,
for feature=5
    
    SumDiff{feature}=mean(Stats_all{feature}.powspctrm(1:nSubjects,:))-mean(Stats_all{feature}.powspctrm((nSubjects+1):end,:));
    colormap('parula')
    ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),(SumDiff{feature}),'refine',1)%, title(sprintf('%s',interLabel{feature}))
    %subplot(2,3,feature),
    % set(gca, 'CLim', [-.15 0.15]);
    % %     hold on
    %     figure,ft_topoplotER(plotcfg,  Inter_topoALL{feature})%, title(sprintf('%s',interLabel{feature}));
    % %     hold on
    % %     ft_plot_lay(lay)
end

%%

figure,
for feature=1:5
    
    SumDiff{feature}=mean(Stats_all{feature}.powspctrm(1:nSubjects,:))-mean(Stats_all{feature}.powspctrm((nSubjects+1):end,:));
    colormap('parula')
    subplot(2,3,feature), ft_plot_topo(electrodes.elecpos(elec_MoneTwo_ind==0,1),electrodes.elecpos(elec_MoneTwo_ind==0,2),(SumDiff{feature}),'refine',1)%,'mask',lay.mask,'outline',lay.outline)%, title(sprintf('%s',interLabel{feature}))
    set(gca, 'CLim', [-.15 0.15]);
    % %     hold on
    %     figure,ft_topoplotER(plotcfg,  Inter_topoALL{feature})%, title(sprintf('%s',interLabel{feature}));
    % %     hold on
    % %     ft_plot_lay(lay)
end
%  InterCOND2=InterFeature_topoALL;
% figure,ft_topoplotER(plotcfg, Feature_topoALL{4}), title(sprintf('%s',featureLabel{4}));
%  figure,ft_topoplotER(plotcfg, InterFeature_topoALL{5}), title(sprintf('%s',featureLabel{1}));
%

% figure, ft
%  InterFeature_topoALL


%save('Intermodul_featuresCOND1.mat','Inter_topoALL','InterFeature_topoALL')

%%

interLabel={'Left && Right', 'Left && Mouth', 'Right && Mouth','All','Blink non-linearity','AVG Mouth-eye pair'};
load('quickcap64.mat') % load 2D layout of electrodes
freqcfg.layout=lay;


%  load Intermodul_featuresCOND1.mat
%    load Intermodul_featuresCOND2.mat
% % % %
%   COND2_topo=Inter_topoALL;
for ii=1:12
    DIFF_topo{ii}=COND1_topo{1};
end

plotcfg.colorbar='South';
plotcfg.zlim=[]%[-2 6];
secndplotcfg=plotcfg;
secndplotcfg.zlim=[-3 4];
figure,
for ii=5
    
    DIFF_topo{ii}.powspctrm=COND1_topo{ii}.powspctrm-COND2_topo{ii}.powspctrm;
    
    subplot(1,3,1), ft_topoplotER(secndplotcfg,  COND1_topo{ii}), title(sprintf('Condition wink/blink',interLabel{ii}));
    subplot(1,3,2), ft_topoplotER(secndplotcfg,  COND2_topo{ii}), title(sprintf('Condition one eye closing detection',interLabel{ii}));
    
    subplot(1,3,3), ft_topoplotER(plotcfg,  DIFF_topo{ii}), title(sprintf('%s',interLabel{ii}));
    
end


figure,  ft_topoplotER(plotcfg,  DIFF_topo{ii}), title(sprintf('%s',interLabel{ii}));
hold on
ft_plot_lay(lay)
% DIFF_topo{8}.powspctrm=(COND1_topo{3}.powspctrm+COND2_topo{3}.powspctrm)/2;
% DIFF_topo{9}.powspctrm=(COND1_topo{2}.powspctrm+COND2_topo{2}.powspctrm)/2;
% DIFF_topo{10}.powspctrm=(COND1_topo{1}.powspctrm+COND2_topo{1}.powspctrm)/2;
% DIFF_topo{11}.powspctrm=DIFF_topo{10}.powspctrm-(DIFF_topo{9}.powspctrm+DIFF_topo{8}.powspctrm)/2;

% DIFF_topo{12}.powspctrm=DIFF_topo{8}.powspctrm-DIFF_topo{9}.powspctrm;
% figure, ft_topoplotER(plotcfg,  DIFF_topo{8}), title(sprintf('%s',interLabel{ii}));
% figure, ft_topoplotER(plotcfg,  DIFF_topo{9}), title(sprintf('%s',interLabel{ii}));
% figure, ft_topoplotER(plotcfg,  DIFF_topo{10}), title(sprintf('%s',interLabel{ii}));
% figure, ft_topoplotER(plotcfg,  DIFF_topo{12}), title(sprintf('%s',interLabel{ii}));

% figure,
% for ii=5
%
% DIFF_topo{ii}.powspctrm=COND1_topo{ii}.powspctrm-COND2_topo{ii}.powspctrm;
% subplot(1,3,1), ft_topoplotER(secndplotcfg,  COND1_topo{ii}), title(sprintf('Condition wink/blink',interLabel{ii}));
% subplot(1,3,2), ft_topoplotER(secndplotcfg,  COND2_topo{ii}), title(sprintf('Condition one eye closing detection',interLabel{ii}));
% subplot(1,3,3), ft_topoplotER(plotcfg,  DIFF_topo{ii}), title(sprintf('%s',interLabel{ii}));
% end
% figure, ft_topoplotER(plotcfg,  DIFF_topo{ii})

electrodes = ft_read_sens('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc')

counter=0;
clear indx_tempToAsa
for electrodes=1:size(electrodes.label,1)
    
    
    
    if any(find(strcmp(electrodes.label(electrodes),DIFF_topo{1}.label)))
        counter=counter+1;
        indx_tempToAsa(counter)=find(strcmp(electrodes.label(electrodes),DIFF_topo{1}.label));
    end
    
end
% figure,    [pos, tri] = icosahedron162;
%     mesh.pos = pos;
%     mesh.tri = tri;
%     ft_plot_mesh(mesh, 'facecolor', 'skin', 'edgecolor', 'none')
%     camlight
%%     hold on
figure, ft_plot_topo3d(electrodes.chanpos(indx_tempToAsa,:) ,DIFF_topo{1}.powspctrm,'contourstyle','black','refine',3)


%%
motus_topo{4}=DIFF_topo{1};
motus_topo{4}.powspctrm=DIFF_topo{1}.powspctrm-(DIFF_topo{3}.powspctrm+DIFF_topo{2}.powspctrm)/2;

figure,ft_topoplotER(plotcfg,  motus_topo{4}), title(sprintf('%s',interLabel{4}));


%% DESCRIPTIVE
cfg.variance      = 'yes'
cfg.jackknife     = 'no'
cfg.keeptrials    = 'yes' % or 'no', estimate single trial power (useful for fourier data) (default = 'no')
% cfg.channel       =
cfg.trials        = 'all' %or a selection given as a 1xN vector (default = 'all')
% cfg.frequency     = %[fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
% cfg.latency       = %[tmin tmax] or 'all', to specify a subset of latencies (default = 'all')

FreqDescriptiveCOND2= ft_freqdescriptives(cfg, COND2_topo{5})

% [grandavg] = ft_freqgrandaverage(cfg, freq{2}{frequency}, freq{3}{frequency},freq{8}{frequency})


figure, plot(grandavg.freq, abs(grandavg.powspctrm));
xlabel('frequency (Hz)')
ylabel('phase-locked amplitude (a.u.)')
legend(freqcfg.channel)
% save ('Intermodulary_frequenciesSubjects.mat','freq')


%% COND 1 VERSUS 2

% WinkNonLinearity=Inter_difference{1};
% load CathWinkTask_1to3.mat
% load CathNOWinkTask_4to7.mat

DiffWink=Winktask;
DiffWink.powspctrm=Winktask.powspctrm(:)' - NOWinkTask.powspctrm(:)';

figure,ft_topoplotER(plotcfg, DiffWink), title(sprintf('%s',interLabel{6}))


% save('CathWinkTask_1to3.mat','WinkNonLinearity','Inter_topo')
%  save('CathNOWinkTask_4to7.mat','WinkNonLinearity','Inter_topo')

% save('Catherine_NONLINEARDIFF.mat','DiffWink')


%% STATISTICS

% In MEEG experiments, there are two types of UO: subjects and trials
% (i.e., epochs of time within a subject). Concerning the concept of a between-
% or a within-UO design, there are two schemes according to which the UOs can be
% assigned to the experimental conditions: (1) every UO is assigned to only one
% of a number of experimental conditions (between UO-design; independent samples),
% or (2) every UO is assigned to multiple experimental conditions in a particular
% order (within UO-design; dependent samples).


% power spectrum for each Subjects, intermodulatory tagg, trial and electrode.
% Inter_difference{subject}{intermodul}.powspctrm(trial,:)

% design(1,design(1,:)==0)=design(1,design(1,:)~=0);
% 1:76= cond1

% load Features_frequenciesSubjects.mat
%  load 'SNRFeatures_subjectxtrialCond1.mat'
load('quickcap64.mat') % load 2D layout of electrodes
% cfg.hdr                 = ft_read_header(cfg.dataFile);
clear statcfg
statcfg.avgovertime = 'no';
cfg.avgoverfreq = 'yes';
statcfg.channel          = 1:64;
statcfg.statistic = 'depsamplesT';%'mean'%'ft_statfun_depsamplesT';%'ft_statfun_indepsamplesT';% 'ft_statfun_actvsblT';
statcfg.correctm         = 'cluster'; % montecarlo : need to put 'max'
statcfg.clusteralpha     = 0.05;
statcfg.clusterstatistic = 'maxsum';
statcfg.minnbchan        = 2;
statcfg.tail             = 0;
statcfg.clustertail      = 0;
statcfg.alpha            = 0.050;
statcfg.numrandomization = 1000;
% statcfg.frequency=freq_F{3}{1}.freq;

% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'distance';
cfg_neighb.layout      =lay;
statcfg.neighbours       = ft_prepare_neighbours(cfg_neighb,Stats_all);

% subj = 7;
% design =[];% zeros(2,2*subj);
% for i = Subj
%     for trial=1:size(Inter_difference{i}{7}.powspctrm(:,:),1)
%   design(1,ii*trial) = i;
%     end
% end
% for i = 1:subj
%   design(1,subj+i) = i;
% end
% design(2,1:subj)        = 1;
% design(2,subj+1:2*subj) = 2;


% design(1,1:111)=1;
% design(1,112:222)=2;
% design(2,1:111)=1:111;
% design(2,112:222)=1:111;
statcfg.resampling ='bootstrap'
statcfg.design   = design';%design;
statcfg.ivar     = 2;
statcfg.uvar     = 1;

statcfg.method='montecarlo'%5'stats' %'analytic';
statcfg.correctm         = 'max';

stat=ft_freqstatistics(statcfg,Stats_all)

cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.zlim   = [];
cfg.layout = lay;
figure,ft_clusterplot(cfg, stat);

%% ANOVA SNR : OCCIPIT x condition
electrodes = ft_read_sens('/home/adf/faghelss/Documents/DATA/DATA_EEG/standard_waveguard64.elc');

ThreeDCap=electrodes.label;
parietal_elec={'Cz','Pz','CPz','P7','P3', 'P4','P1', 'P2', 'P6','P8','P5','CP5','CP6'};
parietal_elec_Left ={'P1','P3','P5','P7','CP5'};
parietal_elec_Right={'P2','P4','P6','P8','CP6'};
parietal_elec_mid={'Cz','Pz','CPz'};
posterior_elec={'Oz','O1','O2','POz','PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8','P5','P6','T7','T8','P3','P4','P7','P8'};%'T7','T8''O1', 'O2''Pz','CPz'
% posterior_elec_Left = {'PO5','P5','TP7','PO7','PO3','O1','T7','P3','P7'};%'T7'
posterior_elec_Left = {'PO7','PO3','PO5','P7','O1','T7','P3','P7','P5','TP7'};%'T8'
% posterior_elec_Right = {'PO6','P6','TP8','PO8','PO4','O2','T8','P4','P8'};%'T8'
posterior_elec_Right = {'PO8','PO4','PO6','P8','O2','T8','P4','P6','P8','TP8'};%'T8'
posterior_elec_mid = {'Oz','POz'};%'POz', 'POz','O2','O1' Pz','CPz'
MoneTwo={'M1','M2'};%'POz', 'POz','O2','O1'
Coh_elec_all={'T7','T8', 'POz', 'O1', 'Oz', 'O2','PO3','PO4','PO5', 'PO6', 'TP7', 'TP8','PO7','PO8','CP5','CP6','FT7','FT8','AF3','AF4'};

% posterior_elec= {'Oz','O1','O2','POz','PO5', 'PO3', 'PO4', 'PO6', 'PO7','PO8','P5','P6','P7','P8'};%'T7','T8''O1', 'O2''Pz','CPz'
% posterior_elec_Left = {'PO5','P5','PO7','PO3','T7','TP7','P7','O1'};%'T7',
% posterior_elec_Right = {'PO6','P6','PO8','PO4','T8','TP8','P8','O2'};%'T8'
% posterior_elec_mid = {'Oz','POz','Pz'};%'POz', 'POz','O2','O1' Pz','CPz'


% posterior_elec= {'PO3','PO5','PO7','P5','O1'};
% posterior_elec= {'PO5'};%'T7','T8''O1', 'O2''Pz','CPz'
% posterior_elec= {'PO2','PO4','PO6','P6'};
% posterior_elec_Left = {'PO3','PO5','PO7','P5'};%'T7',
% % posterior_elec_Left = {'PO2','PO4','PO8','P6'};%'T7',
% posterior_elec_Right = {'O2','PO4','PO6','P6'};%
% posterior_elec_mid = {'Oz','O1'};%'POz', 'POz','O2','O1' Pz','CPz'

% posterior_elec_Temp = {'T7','T8','TP7','TP8'};
% posterior_elec_Left = {'O1','P7','P3','P5', 'P1','PO5', 'PO3','PO7'};
% posterior_elec_Right = {'O2','P8','P4','P6', 'P2','PO6', 'PO4','PO8'};
%  Stats=data_postica{12}{1};
%  Stats.powspctrm=Stats_all{4}.powspctrm;%(Stats_all{1}.powspctrm+Stats_all{2}.powspctrm)/2;
% Stats.powspctrm=(Stats_all{1}.powspctrm+Stats_all{2}.powspctrm)/2;
%     Stats.powspctrm=(Stats_all{3}.powspctrm)./((Stats_all{3}.powspctrm+Stats_all{2}.powspctrm)+(Stats_all{1}.powspctrm));
%     Stats.powspctrm=(Stats_all{8}.powspctrm)./((Stats_all{3}.powspctrm+Stats_all{2}.powspctrm)+(Stats_all{1}.powspctrm));%-Stats_all{1}.powspctrm;%+Stats_all{3}.powspctrm)/2);%-(Stats_all{1}.powspctrm);
%    Stats.cohspctrm=squeeze(mean((CohFeature{1}.cohspctrm(:,31,:)),2));%-(CohInterFeature{3}.cohspctrm(:,:,:));%zscore((Stats_all{3}.powspctrm))%-zscore(STATS_FEATURES{3}.powspctrm);%-zscore((STATS_FEATURES{1}.powspctrm+STATS_FEATURES{2}.powspctrm)/2);
% Stats.powspctrm=InterFeature_topo_wtvFeat_all{1}.powspctrm;
% Stats.powspctrm=(Stats_all{10}.powspctrm)%;+Stats_all{2}.powspctrm)/2;%+Stats_all{5}.powspctrm;%POWattendedFeatures.powspctrm(:,:);%((Stats_all{2}.powspctrm))%+(Stats_all{2}.powspctrm))/2;%+(Stats_all{4}.powspctrm))%-((Stats_all{3}.powspctrm+Stats_all{2}.powspctrm))/2);

% Stats.powspctrm=((Stats_all{4}.powspctrm))-((Stats_all{5}.powspctrm+Stats_all{6}.powspctrm))/2);%+Stats_all{2}.powspctrm+Stats_all{10}.powspctrm;%./((Stats_all{1}.powspctrm+Stats_all{2}.powspctrm)/2);
%      Stats.powspctrm=Stats_all_band{2}.powspctrm%InterFeature_topo_wtvFeat_all{8}.powspctrm%((Stats_all{2}.powspctrm));%+Stats_all{2}.powspctrm+Stats_all{3}.powspctrm))/3;%SNR_eyesAllcond;%

% % % % Stats.powspctrm=(Stats_all{1}.powspctrm+Stats_all{2}.powspctrm)-(Stats_all{3}.powspctrm);
% clear BrainSign_all
% Subj=Subj_cond2;
% whichSub=((17+1):(17+length(Subj)));
% % whichSub(2)=[];
% % whichSub(11)=[];
% subCounter=0;
% for subject=whichSub
%     subCounter=subCounter+1;
% BrainSign_all(subCounter)=mean(mean(Stats.powspctrm(subject,find(elec_post_ind))));
% end
% %
% [r,p]=corrcoef(RT_allcond{2}(:),BrainSign_all)
% % CORRELATE NON_LINEARUTY TO INDIVIDUAL PERFORMANCE (RT)
% figure, scatter(RT_allcond{2}(:),BrainSign_all),lsline, title(sprintf('corr between brain signal and RT : %.2f, p=%.2f',r(2),p(2)))
% % % % % %
%
% data_postica=raw_data;
%  Stats.label=elec_montage.label;%data_postica.label;
Stats.label=freq{22}.label%freq{9}.label;
%
% mean((mean(Stats.powspctrm(1:nSubjects,find(elec_postLeft_ind)))))
nSubjects=14;
%

clear elec_post_ind elec_postLeft_ind elec_postRight_ind elec_postMid_ind
for elec=1:size(Stats.label,1)
    
    if any(find(strcmp(Stats.label(elec),posterior_elec)))
        elec_post_ind(elec)= 1;
    else
        elec_post_ind(elec)= 0;
    end
    
    if any(find(strcmp(Stats.label(elec),posterior_elec_Left)))
        elec_postLeft_ind(elec)= 1;
    else
        elec_postLeft_ind(elec)= 0;
    end
    
    if any(find(strcmp(Stats.label(elec),posterior_elec_Right)))
        elec_postRight_ind(elec)= 1;
    else
        elec_postRight_ind(elec)= 0;
    end
    if any(find(strcmp(Stats.label(elec),posterior_elec_mid)))
        elec_postMid_ind(elec)= 1;
    else
        elec_postMid_ind(elec)= 0;
    end
    if any(find(strcmp(Stats.label(elec),parietal_elec)))
        elec_parietal_ind(elec)= 1;
    else
        elec_parietal_ind(elec)= 0;
    end
    if any(find(strcmp(Stats.label(elec),parietal_elec_Left)))
        elec_parietalLeft_ind(elec)= 1;
    else
        elec_parietalLeft_ind(elec)= 0;
    end
    if any(find(strcmp(Stats.label(elec),parietal_elec_Right)))
        elec_parietalRight_ind(elec)= 1;
    else
        elec_parietalRight_ind(elec)= 0;
    end
    if any(find(strcmp(Stats.label(elec),parietal_elec_mid)))
        elec_parietalMid_ind(elec)= 1;
    else
        elec_parietalMid_ind(elec)= 0;
    end
    
    if any(find(strcmp(ThreeDCap(elec),MoneTwo)))
        elec_MoneTwo_ind(elec)= 1;
    else
        elec_MoneTwo_ind(elec)= 0;
    end
    
    if any(find(strcmp(Stats.label(elec),Coh_elec_all)))
        Coh_elec_all_ind(elec)= 1;
    else
        Coh_elec_all_ind(elec)= 0;
    end
    
end
elec_MoneTwo_ind(13)=1;
elec_MoneTwo_ind(19)=1;
elec_MoneTwo_ind(60:64)=0;


save('electrodes_indexes.mat','elec_post_ind','elec_postLeft_ind','elec_postRight_ind','elec_postMid_ind','elec_parietal_ind','elec_parietalLeft_ind','elec_parietalLeft_ind','elec_parietalMid_ind')
%%
counter1=0;

clear Elec_matrix
for elec=1:size(Stats.label,1)
    
    if (elec_post_ind(elec)+elec_postLeft_ind(elec))==2
        counter1=counter1+1;
        Elec_matrix(counter1)= 1;
    end
    if (elec_post_ind(elec)+elec_postRight_ind(elec))==2
        counter1=counter1+1;
        Elec_matrix(counter1)= 3;
    end
    if (elec_post_ind(elec)+elec_postMid_ind(elec))==2
        counter1=counter1+1;
        Elec_matrix(counter1)= 2;
    end
    %
    %     if (elec_parietal_ind(elec)+elec_parietalLeft_ind(elec))==2
    %         counter1=counter1+1;
    %         Elec_matrix(counter1)= 1;
    %     end
    %     if ( elec_parietal_ind(elec)+elec_parietalRight_ind(elec))==2
    %         counter1=counter1+1;
    %         Elec_matrix(counter1)= 3;
    %     end
    %     if (elec_parietal_ind(elec)+elec_parietalMid_ind(elec))==2
    %         counter1=counter1+1;
    %         Elec_matrix(counter1)= 2;
    %     end
    
end
%
%find(elec_postMid_ind)
%  Stats.cohspctrm=[];
%  Stats.cohspctrm(1:nSubjects,:)=squeeze(mean(abs(CohInterFeature{1}{4}.cohspctrm(:,29:31,:)),2));
%  Stats.cohspctrm((nSubjects+1):(2*nSubjects),:)=squeeze(mean(abs(CohInterFeature{2}{4}.cohspctrm(:,29:31,:)),2));
%
%  Stats.cohspctrm(1:nSubjects,:)=squeeze((Stats.powspctrm(:,:)));
%   Stats.cohspctrm((nSubjects+1):(2*nSubjects),:)=Stats_all{2}.powspctrm;%squeeze((Stats.cohspctrm(:,18,:)));
%   Stats.cohspctrm((nSubjects+10):(3*nSubjects),:)=Stats_all{3}.powspctrm;%squeeze((Stats.cohspctrm(:,18,:)));
%
Stats.cohspctrm=Stats.powspctrm;
%
clear Subjects_Matrix Electrode_Matrix Condition_Matrix vd
%[p, table1] = anova_rm({Stats_all.powspctrm(1:8,find(elec_post_ind)),Stats_all.powspctrm(9:16,find(elec_post_ind))});% Between two groups and 3 first, non-induction blocks

vd=Stats.powspctrm(:,find(elec_post_ind));
% vd=Stats.cohspctrm(:,find(elec_post_ind));


Condition_Matrix=zeros(size(vd));
% Condition_Matrix(1:18,:)=1;
% Condition_Matrix((18+1):(32),:)=2;
% Condition_Matrix((33):40,:)=3;
Condition_Matrix(1:nSubjects,:)=1;
Condition_Matrix(((nSubjects+1):(2*nSubjects)),:)=2;
% Condition_Matrix(((2*nSubjects)+1):(3*nSubjects),:)=3;
% Condition_Matrix(((3*nSubjects)+1):(4*nSubjects),:)=2;
% Condition_Matrix(((4*nSubjects)+1):(5*nSubjects),:)=3;
% Condition_Matrix(((5*nSubjects)+1):(6*nSubjects),:)=3;

% Condition_Matrix(1:18,:)=1;
% Condition_Matrix((19):(32),:)=2;
% Condition_Matrix((33):(40),:)=3;
%
for ii=1:nSubjects*2
    %    Subjects_Matrix(ii+nSubjects,:)=ones(1,size(vd,2)).*ii;
    %    Subjects_Matrix(ii,:)=ones(1,size(vd,2)).*ii;
    %     Electrode_Matrix(ii,:)=Elec_matrix;
    %     Electrode_Matrix(ii+nSubjects,:)=Elec_matrix;%repmat(Elec_matrix,size(vd,1),1);
    %     Electrode_Matrix(ii+nSubjects*2,:)=Elec_matrix;
    
    
    Electrode_Matrix(ii,:)=Elec_matrix;
    %    RT_Matrix(
end

%%anovan(vd(:),{Condition_Matrix(:), Subjects_Matrix(:),Electrode_Matrix(:)},'model','full','varnames',{'Condition','Subjects','Electrodes'})
% [P,T,STATS,TERMS]= anovan(vd(:),{Electrode_Matrix(:)},'model','full','varnames',{'Electrodes'},'display','on')
[P,T,STATS,TERMS]= anovan(vd(:),{Electrode_Matrix(:), Condition_Matrix(:)},'model','full','varnames',{'Electrodes','Condition'},'display','on')
%
% %
% %
figure,

multcompare(STATS,'Dimension',[2],'Display','on')
%
%%
% [h,p1{1}]=ttest2(vd(Electrode_Matrix==1),vd(Electrode_Matrix==2))
% [h,p2{1}]=ttest2(vd(Electrode_Matrix==1),vd(Electrode_Matrix==3))
% [h,p3{1}]=ttest2(vd(Electrode_Matrix==2),vd(Electrode_Matrix==3))
%
% figure, bar([nanmean(vd(Electrode_Matrix==1)) nanmean(vd(Electrode_Matrix==2)) nanmean(vd(Electrode_Matrix==3))])%,axis([ 0 4 1 1.8])
% sigstar({[1,2] [1,3] [2,3]},[p1{1} p2{1} p3{1}])
% set(gca,'XTick',1:3,'XTickLabel',{'Left hemisphere electrodes','central','Right hemisphere'})

[h,p1{1}]=ttest2(vd(Condition_Matrix==1),vd(Condition_Matrix==2))
[h,p2{1}]=ttest2(vd(Condition_Matrix==1),vd(Condition_Matrix==3))
[h,p3{1}]=ttest2(vd(Condition_Matrix==2),vd(Condition_Matrix==3))
% [h,p4{1}]=ttest2(vd(Condition_Matrix==2),vd(Condition_Matrix==3))
% [h,p5{1}]=ttest2(vd(Condition_Matrix==2),vd(Condition_Matrix==4))
% [h,p6{1}]=ttest2(vd(Condition_Matrix==3),vd(Condition_Matrix==4))

figure, bar([nanmean(vd(Condition_Matrix==1)) nanmean(vd(Condition_Matrix==2)) nanmean(vd(Condition_Matrix==3)) ]), title(sprintf('Eye Intermodulate Response between conditions,p=%.3f',P(2)))%, axis([ 0 4 1 1.2])
% figure, bar([nanmean(vd(Condition_Matrix==1)) nanmean(vd(Condition_Matrix==2)) nanmean(vd(Condition_Matrix==3))  ]), axis([ 0 4 3 5])% title(sprintf('Eye Intermodulate Response between conditions,p=%.3f',P(2))), axis([ 0 3 7 11])
% sigstar({[1,2] [1,3] [1,4] [2,3] [3,4]},[p1{1} p2{1} p3{1} p4{1} p6{1}])
% sigstar({[1,2] [1,3],[2,3]},[p1{1} p2{1} p3{1}])
sigstar({[1,2] [1,3] [2,3]},[p1{1} p2{1} p3{1}])
% set(gca,'XTick',1:3,'XTickLabel',{'Left eye','Right eye','Mouth'}),ylabel(' Signal to Noise Ration (SNR)')
set(gca,'XTick',1:2,'XTickLabel',{'Non-linear (Wink/Blink) Condition','Linear (Blink only) condition'}),ylabel(' Signal to Noise Ration (SNR,%)')
% nanmean(vd(Condition_Matrix==3))

% figure, boxplot(mean(Stats.powspctrm(1:nSubjects,find(elec_post_ind))))
% hold on
% boxplot(mean(Stats.powspctrm((nSubjects+1):(2*nSubjects),find(elec_post_ind))),'--r')
% hold on
% boxplot(nanmean(Stats.powspctrm(((nSubjects*2)+1):(3*nSubjects),find(elec_post_ind))),'-k')

%%
figure, boxplot([(Stats.powspctrm(Condition_Matrix==1))  (Stats.powspctrm(Condition_Matrix==2)) ],'notch' ,'on')
% hold on
% boxplot(mean(Stats.powspctrm((nSubjects+1):(2*nSubjects),find(elec_post_ind))),'--r')
% hold on
% boxplot(nanmean(Stats.powspctrm(((nSubjects*2)+1):(3*nSubjects),find(elec_post_ind))),'-k')
%
%
%
figure,

multcompare(STATS,'Dimension',[1 2],'Display','on')


%%
anova1(vd(:),Condition_Matrix(:))
%
anova1(vd(:),Electrode_Matrix(:))
%  anovan(vd(:),{Condition_Matrix(:)},'model','full','varnames',{'Condition'})

%
% condition1=Stats;
% condition2=Stats;
% condition1.powspctrm=Stats_all{5}.powspctrm(1:nSubjects,:);
% condition2.powspctrm=Stats_all{5}.powspctrm((nSubjects+1):(2*nSubjects),:);

% colormap='parula';
% dataVector=((Stats_all{1}.powspctrm((nSubjects+1):(2*nSubjects),:)));
% dataVector
% eloc = readlocs('BioSemi64.loc');
% figure;topoplot(dataVector,eloc,'style','blank','colormap',colormap('default'))
% set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
%
% figure, plot(squeeze(mean(mean(vd(1:nSubjects,Elec_matrix==1,:)))),'--r')
% hold on
% plot(squeeze(mean(mean(vd((nSubjects+1):(2*nSubjects),Elec_matrix==1,:)))),'--b')
% hold on
% plot(squeeze(mean(vd(1:nSubjects,1,:)))-squeeze(mean(vd((nSubjects+1):(2*nSubjects),1,:))),'g')
% set(gca,'XTick',1:23,'XTickLabel',posterior_elec)

figure, topoplot(mean(Stats.cohspctrm(1:8,:))-mean(Stats.cohspctrm(9:16,:)),'BioSemi64.loc'),colorbar
%%
figure
for row=1:4
    for col=1:4
        subplot(4,4,(row-1)*4+col),bar(squeeze(mean(vd(1:8,row,col))))
        ylim([0 1])
    end
end


% figure, plot(mean(vd(1:nSubjects,:)),'--r')
% hold on
% plot(mean(vd((nSubjects+1):(2*nSubjects),:)),'--b')
% hold on
% plot(mean(vd(1:nSubjects,:))-mean(vd((nSubjects+1):(2*nSubjects),:)),'g')
% set(gca,'XTick',1:23,'XTickLabel',posterior_elec)
% figure,  ft_topoplotER(plotcfg,  Stats)
%     hold on
%     ft_plot_lay(lay)

% mean((mean(Stats.powspctrm(1:nSubjects,find(elec_postRight_ind)))))
% EEGloc=load('BioSemi64.loc')
% figure; topoplot([],'BioSemi64.loc','style','blank','electrodes','labelpoint')
% figure, mesh(mean((Stats_all{1}.powspctrm((nSubjects+1):(2*nSubjects),:))))
%% MULTIPLE COMPARISONS BETWEEN ELECTRODES
clear tmap1 stats pmap1 Hmap1 p H
% tmap1=zeros(size(CohInterFeature{1}{1}.wplispctrm,2));

% tmap2=zeros(1,size(CohInterFeature{1}.cohspctrm(:,1,1:64),2));
n1=0;
n2=0;
feat=2;
interm=1;
perm=[1:6]
% COND1Coh{5}=abs(COND1.iCohInter{1})-(abs(COND1.iCohInter{3})+abs(COND1.iCohInter{5}))/2;
% COND2Coh{5}=abs(COND2.iCohInter{1})-(abs(COND2.iCohInter{3})+abs(COND2.iCohInter{5}))/2;
% COND1Coh{5}=squeeze(nanmean(COND1Coh{5},2));
% COND2Coh{5}=squeeze(nanmean(COND2Coh{5},2));
% COND1Coh{feat}=((squeeze(nanmean(abs(COND1.iCohInter{feat})+abs(COND1.iCohInter{feat+1}),2))));%-(((squeeze(nanmean(abs(COND1.iCohInter{3})+abs(COND1.iCohInter{4}),2))))+((squeeze(nanmean(abs(COND1.iCohInter{5})+abs(COND1.iCohInter{6}),2))))/2);
% COND2Coh{feat}=(squeeze(nanmean(abs(COND2.iCohInter{feat})+abs(COND2.iCohInter{feat+1}),2)));%-(((squeeze(nanmean(abs(COND2.iCohInter{3})+abs(COND2.iCohInter{4}),2))))+((squeeze(nanmean(abs(COND2.iCohInter{5})+abs(COND2.iCohInter{6}),2))))/2);
% COND1Coh{feat}=((squeeze(nanmean(abs(COND1.iCohInter{1})+abs(COND1.iCohInter{2}),2))))-(((squeeze(nanmean(abs(COND1.iCohInter{3})+abs(COND1.iCohInter{4}),2))))+((squeeze(nanmean(abs(COND1.iCohInter{5})+abs(COND1.iCohInter{6}),2))))/2);
% COND2Coh{feat}=(squeeze(nanmean(abs(COND2.iCohInter{1})+abs(COND2.iCohInter{2}),2)))-(((squeeze(nanmean(abs(COND2.iCohInter{3})+abs(COND2.iCohInter{4}),2))))+((squeeze(nanmean(abs(COND2.iCohInter{5})+abs(COND2.iCohInter{6}),2))))/2);


% COND1cohFeat=squeeze(nanmean(abs(COND1.iCohFeat{feat})))+squeeze(nanmean(abs(COND1.iCohFeat{feat+1})));
% COND2cohFeat=squeeze(nanmean(abs(COND2.iCohFeat{feat})))+squeeze(nanmean(abs(COND2.iCohFeat{feat+1})));
% COND1cohFeat=COND1cohFeat-squeeze(nanmean(abs(COND1.iCohFeat{feat+2})));
% COND2cohFeat=COND2cohFeat-squeeze(nanmean(abs(COND2.iCohFeat{feat+2})));
%
% EyesVsOthersCond1=iCohInterCond1All{1}-(iCohInterCond1All{2}+iCohInterCond1All{3})/2;
% EyesVsOthersCond2=iCohInterCond2All{1}-(iCohInterCond2All{2}+iCohInterCond2All{3})/2;
nbElec=find(elec_post_ind);
h= waitbar(0, 'ttest 1 : 0 % compl�t�');
jj=1;
% for ii= nbElec%find(channIdx)%find(elec_post_ind)%find(channIdx)
%     n1=n1+1;
for ii=nbElec %find(channIdx)%find(elec_post_ind)%
    
    
    
    [H{ii,jj},p{ii,jj},~,stats{ii,jj}]=ttest(nanmean(COND_IM{4}{interm}(:,perm,ii),2),nanmean(COND_IM{2}{interm}(:,perm,ii),2));
    tmap1(ii,jj) = stats{ii,jj}.tstat;
    pmap1(ii,jj) =p{ii,jj};
    Hmap1(ii,jj) = H{ii,jj};
end
%     waitbar((ii/nbElec), h, sprintf('ttest : %3.2f %% compl�t�', (ii/14)));
% end
delete(h)
pmap1=(pmap1);
tmap1=(tmap1);
[pthr,pcor,padj]=fdr(pmap1(pmap1>0));

pmap1(pmap1>0)=pcor;
% pmap1(pmap1>.05)=nan;
% pmap1(pmap1<=.05&(tmap1)>0)=1;
% pmap1(pmap1<=.05&(tmap1)<0)=-1;

%  figure, imagesc(squeeze(tmap1)),colorbar

figure, subplot(1,2,1),imagesc(squeeze(tmap1)),colorbar
% figure, imagesc(squeeze(pmap1)),colorbar
% set(gca,'XTick',1:nbElec,'XTickLabel',Coh_elec_all)
% set(gca,'YTick',1:nbElec,'YTickLabel',Coh_elec_all)
%  subplot(1,2,2), imagesc((squeeze(pmap1)<=(.05/(length(pmap1(pmap1>0)))))),colorbar %&squeeze(tmap1)>0
subplot(1,2,2), imagesc((squeeze(pmap1)<=.05&squeeze(tmap1)>0)) %&squeeze(tmap1)>0
% subplot(1,2,2), imagesc(pmap1),colorbar %&squeeze(tmap1)>0

% figure, imagesc(squeeze(pmap1)),colorbar
% set(gca,'XTick',1:nbElec,'XTickLabel',Coh_elec_all)
% set(gca,'YTick',1:nbElec,'YTickLabel',Coh_elec_all)
%  contour(pmap1,[.00001 .05],'r')
% hold on
% contour(Hmap1==1,'r')
%
% figure, imagesc(Hmap1),colorbar
% hold on
%% SHOW COHERENCE ON

FAKETOPO=zeros(1,64);
FAKETOPO(Coh_elec_all_ind==1)=squeeze(nanmean(COND1Coh{feat}(:,1,:)-COND2Coh{feat}(:,1,:)));
squeeze(nanmean(COND1Coh{feat}(:,1,:)-COND2Coh{feat}(:,1,:)))

plotcfg.colorbar='South';
plotcfg.zlim=[-.4 .25];
plotcfg.elec=electrodes;

for feature=1
    
    figure,
    ft_plot_topo3d(electrodes.elecpos(:,:,:,:),FAKETOPO),title(sprintf('%d',feature)),colorbar;
    %  subplot(1,2,1),
    %  set(gca, 'CLim', [0 0.1]);
    %  subplot(1,2,2),ft_plot_topo3d(electrodes.elecpos(elec_MoneTwo_ind==0,:,:,:),mean(Stats_all{feature}.powspctrm((nSubjects+1):end,:)),'refine',2),title(sprintf('%d',feature));
    %  set(gca, 'CLim', [0 0.1]);
end
%     figure,ft_topoplotER(plotcfg,  Inter_topoALL{feature})%, title(sprintf('%s',interLabel{feature}));
% %     hold on
% %     ft_plot_lay(lay)
% end
%


%%
clear propMax maxIcohChannel allPropMax Y
for feat=[1]
    for subject=1:8
        propMax{feat}(subject,:,:)=zeros(14);
        [Y,maxIcohChannel{feat}(subject)]=max(COND1.iCohInter{feat}(subject,:));
        propMax{feat}(subject,maxIcohChannel{feat}(subject))=1;
    end
    allPropMax{feat}=squeeze(sum(propMax{feat}));
    figure, imagesc(allPropMax{feat}),colorbar
    set(gca,'XTick',1:14,'XTickLabel',toChann),xlabel('To channel')
    set(gca,'YTick',1:14,'YTickLabel',posterior_elec),ylabel('From channel')
end


figure, subplot(1,2,1),imagesc(squeeze(nanmean(iCohInterCond1All{1})),[0 .6])
set(gca,'XTick',1:14,'XTickLabel',toChann),xlabel('To channel')
set(gca,'YTick',1:14,'YTickLabel',posterior_elec),ylabel('From channel')
subplot(1,2,2),imagesc(squeeze(nanmean(iCohInterCond2All{1})),[0 .6])
set(gca,'XTick',1:14,'XTickLabel',toChann),xlabel('To channel')
set(gca,'YTick',1:14,'YTickLabel',posterior_elec),ylabel('From channel')
%
% h= waitbar(0, 'ttest 1 : 0 % compl�t�');
% for ii= 1:62
%
%         [~,~,~,stats{ii}]=ttest2(EEGtimeALL(2).data(ii,jj,:),Stats_all{4}.powspctrm(9:16,ii,1));
%         tmap2(ii) = stats{ii}.tstat;
%
%     waitbar(ii/6, h, sprintf('ttest : %3.2f %% compl�t�', (ii/64)*100));
% end
% delete(h)
%
% figure, imagesc(tmap2),colorbar

%%

% for ii=1:66
%
%     if any(isnan(CohInterFeature{1}{1}.wplispctrm(1:8,ii,1))) %|| ii~= find(elec_post_ind)
%         channIdx(ii)=0;
%     else
%         channIdx(ii)=1;
% %     elseif
%     end
% end
clear tmap2 pmap2 stats
% nbSubject=8;
counter=0;
h= waitbar(0, 'ttest 1 : 0 % compl�t�');
for ii= 1:find(elec_post_ind)
    counter=counter+1;
    [~,p{ii},~,stats{ii}]=ttest2(Stats_all{1}.powspctrm(:,ii),Stats_all{3}.powspctrm(:,ii));
    tmap2(ii) = stats{ii}.tstat;
    pmap2(ii)= p{ii};
    
    waitbar(ii/6, h, sprintf('ttest : %3.2f %% compl�t�', (ii/64)*100));
end
delete(h)

figure, imagesc(tmap2),colorbar
set(gca,'XTick',1:17,'XTickLabel',posterior_elec)
figure, imagesc(pmap2<.05),colorbar


%%
figure, imagesc(tmap2),colorbar
clear plotcfg

plotcfg.channel=cfg.hdr.label(1:64)%freqcfg.channel;
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='East';
plotcfg.layout=lay;
% plotcfg.foilim=[;
plotcfg.showlabels ='yes';
plotcfg.zlim = []%[0 0.5];


% make a plot
cfg = [];
cfg.highlightsymbolseries = ['*','*','.','.','.'];
cfg.layout = lay%'CTF151.lay';
cfg.contournum = 0;
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.zlim = [-5 5];
figure, ft_clusterplot(cfg,stat);


%% Mutli channel plot of timex frequency analysis

timecfg = [];

timecfg.showlabels   = 'yes';
timecfg.layout       = lay;
timecfg.channel      = 'Oz';
figure,
ft_singleplotTFR(timecfg, freq);




%% Granger
Connectcfg.foi= 8;
Connectcfg.method='granger';%'coh', 'granger'

stat = ft_connectivityanalysis(Connectcfg,freq{subject}{1})


figure
for row=1:4
    for col=1:4
        subplot(4,4,(row-1)*4+col);
        plot(stat.freq, squeeze(stat.grangerspctrm(row,col,:)))
        ylim([0 1])
    end
end

newcfg.parameter = 'grangerspctrm';
newcfg.zlim      = [0 1];
figure, ft_connectivityplot(newcfg, stat);

%% Coherence

% causal outflow hub? coordinating other large-scale networks
% (Sridharan et al., 2008; Menon and Uddin, 2010).
% In networks, hubs allow for increased levels of information flow between
% distant nodes, acting as integrators (Honey et al., 2007).
%
% A causal outflow hub is a region that has a high number of causal outflow
% connections and low number of causal inflow connections as assessed using
% Granger causality analyses.


% Nolte et al., 2004  : In general, if the imaginary part of C(x,y) is positive,
% then x and y are interacting and x is earlier than y; indicating
% that information is flowing from x to y: At specific
% frequencies, however, ?earlier? and ?later? are ambigious;
% e.g. at 10 Hz 10 ms earlier is the same as 90 ms later. For the
% present interpretation we assumed that the smaller delay in
% absolute value is the more probable explanation; e.g. in the
% above example we would favor ?10 ms earlier? over ?90 ms
% later?. Note, that we can make this interpretation just from
% the sign of the imaginary part of coherency without actually
% calculating a delay for which we would need a reliable real
% part of coherency

%%%
% Dear Maxime,
%
% simply use the second calling format as described in the help message
% of newcrossf "help newcrossf". To compare the coherence between
% channels 10 and 11 between datasets 1 and 2 type in:
%
% newcrossf( { ALLEEG(1).data(10,:,:) ALLEEG(2).data(10,:,:) },
% {ALLEEG(1).data(11,:,:) ALLEEG(2).data(11,:,:)}, ALLEEG(1).pnts,
% [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ALLEEG(1).srate);

% find(elec_post_ind)
%  newcrossf({ EEGtimeALL(1).data(30,:,:)  EEGtimeALL(2).data(31,:,:)},...
% 17152,[100 67000],256,0)
%
{ EEGtimeALL(1).data(58,:,:)  EEGtimeALL(2).data(58,:,:)}
newcrossf({ EEGtimeALL(1).data(32,:,:)  EEGtimeALL(2).data(32,:,:)},...
    { EEGtimeALL(1).data(14,:,:)  EEGtimeALL(2).data(14,:,:)},...
    1152,[0 4500],256,0,'maxfreq',50)%,'baseline',[10 300])
% newcrossf({ mean(EEGtimeALL(1).data(29:31,:,:),1)  mean(EEGtimeALL(2).data(29:31,:,:),1)},...
% { mean(EEGtimeALL(1).data(find(elec_postRight_ind),:,:))  mean(EEGtimeALL(2).data(find(elec_postRight_ind),:,:))},256,[0 1000],256,0,'maxfreq',50,'condboot','abs','boottype','randall','type','amp')%,'base 'padratio' 'elocs','BioSemi64.loc'
%



% { EEGtimeALL(1).data(58,:,:)  EEGtimeALL(2).data(58,:,:)}
newcrossf({ EEGtimeALL(1).data(32,:,:)  EEGtimeALL(2).data(32,:,:)},...
    { EEGtimeALL(1).data(14,:,:)  EEGtimeALL(2).data(14,:,:)},...
    1024,[0 2048],1024,0,'maxfreq',20)%,'baseline',[10 300])


n=0;
for subject=Subj_cond1
    
    COND='BLINK';
    eval(sprintf('load(''TrialedEPOCHED%s_data%s_%s'',''raw_data'')',COND,Subjects{subject}{1},Subjects{subject}{5}))
    
    trials=size(raw_data.trial,2);
    for trial=1:trials
        n=n+1;
        EEGtimeALL(2).data(:,:,n)= raw_data.trial{trial}(:,1:3072);  % or 17152
    end
    
    
    
end
% %30,18
%
% newcrossf({ (EEGtimeALL(1).data(29,:,:))  (EEGtimeALL(2).data(29,:,:))},...
% { (EEGtimeALL(1).data([18],:,:))  (EEGtimeALL(2).data([18],:,:))},256,[0 1000],256,0,'maxfreq',25,'condboot','abs','boottype','randall','alpha',.05)%,'base 'padratio' 'elocs','BioSemi64.loc'
%
%
% newcrossf({ (EEGtimeALL(1).data(29,:,:))  (EEGtimeALL(2).data(29,:,:))},256,[0 1024],256,0,'maxfreq',25,'condboot','abs','boottype','randall','alpha',.05)%,'base 'padratio' 'elocs','BioSemi64.loc'


% newcrossf({ EEGtimeALL(1).data(52,:,:)  EEGtimeALL(2).data(52,:,:)},...
% { EEGtimeALL(1).data(54,:,:)  EEGtimeALL(2).data(54,:,:)},1152,[0 4500],256,0,'maxfreq',20,'condboot','abs','boottype','shuffle','alpha',.05,'elocs','BioSemi64')%,'base 'padratio'
% % load Features_frequenciesSubjects.mat
fs=1024;
load('quickcap64.mat') % load 2D layout of electrodes
% TrialTime=(size(raw_data.time,1)/1024);
% figure, mscohere(freq.powspctrm(5,:),freq.powspctrm(2,:),[],[],[],fs,'onesided')
% figure,
% for row=1:8%size(freq.powspctrm,1)
% for col=1:8%size(freq.powspctrm,1)
%   subplot(8,8,(row-1)*8+col);
%   mscohere(freq.powspctrm(row,:),freq.powspctrm(col,:),hanning(1024),[],[],fs,'onesided')
%   ylim([0 1])
% end
% end
% % % % % %
% load Intermodul_featuresCOND1.mat
%  cfg.hdr.label(1:64)


%%%%%%%%% CHECK OUT pac.m  AND crossf AND timefreq!!!!%%%%%%
% Subj_cond1=[ 8 12 16 18 20 22 24 26 28] % Cond 1 only left out 14
% Subj_cond2=[ 9 13 17 19 21 23 25 27 29] % Cond 2 only left out 15
%%
clear freqcfg freq cff cohcfg

cohcfg            = [];
cohcfg.output     = 'fourier';%;'powandcsd';
cohcfg.method     = 'mtmfft';
cohcfg.taper ='dpss';
cohcfg.tapsmofrq  = 0.015;
cohcfg.foilim     = [1 20];

cfg.pad='nextpow2';

cohcfg.tapsmofrq  = 0.015;
cohcfg.keeptrials='yes';


cfg.channel    =[1:64];%freqcfg.c

%%

CONDITION=1;
% Subj=[Subj_cond2(end-1)]
% Subj=Subj_cond1((end-7):end);/
Subj=2:3
clear fd freqPowNcsd AllcohChan
nFrequencies=length(cohcfg.foilim(1):cohcfg.tapsmofrq:cohcfg.foilim(2));
n=0;
cohAll=[];
feat= [3.5290    8.0000    6.3150];
Perm = perms(feat);
nSub=0;
% clear iCohFeat iCohFeatBOOT iCohInter iCohInterBOOT

for subject=Subj
    nSub=nSub+1
    %     eval(sprintf('load(''Trialed_PostIcadata%s_%s.mat'')',Subjects{subject}{1},Subjects{subject}{5}))
    %     eval(sprintf('load(''Trialed_PostIcadataSMOOTHRej%d_%s'')',subject,Subjects{subject}{5}))
    
    if subject<18
         eval(sprintf('load(''Trialed70_avg_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
    else
        eval(sprintf('load(''Trialed_avg_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
    end
%     eval(sprintf('load(''Trialed_data%s_%s.mat'')',Subjects{subject}{1},Subjects{subject}{5}))
    data_postica=raw_data;
    %     raw_data=data_postica;
    
    
    
    clear wav wavChanAll
    Fs=256;
    startTrial=1*Fs;
    endTrial=66*Fs;
    
    nbTrial=max(size(data_postica.trial));
    for trial=1:nbTrial
        % wav(:,trial,1)=raw_data.trial{trial}(18,:);
        %         wav(:,trial,2)=data_postica.trial{trial}(31,:);
        for channel=1:64
            wavChanAll(:,trial,channel)=data_postica.trial{trial}(channel,:);%startTrial:endTrial
        end
    end
    N=size(wavChanAll(:,trial,channel),1);
    
    % Frequency specifications:
    dF = Fs/N;                      % hertz
    f = -Fs/2:dF:Fs/2-dF;
    
    
    %     %
    %
    %        freqPowNcsd{subject}        = ft_freqanalysis(cohcfg, data_postica);
    
    
    
    %     cfg.method     = 'pdc';%'wpli';%'coh';%'wpli_debiased';%'coh';   %'csd'%;'corr'%;'granger'%'powcorr'%'plv'%'coh';%'wpli_debiased';%'powcorr_ortho';%'plv'%'coh';'powcrr'
    % %      cfg.complex    = 'real';
    %     cfg.removemean  = 'no';
    %     cfg.bandwidth = 0.45;
    clear permidx
    for perm=1:6
        fd{subject}{perm}=[];
        if any(TrialTagg{subject}==perm)
            cfg.trials=find(TrialTagg{subject}==perm);%'all'
            cfg.trials(cfg.trials>max(size(data_postica.trial)))=[];
        else
            continue
        end
        
        if length(cfg.trials)>1
            permidx(perm)=1;
            nT=0;
            for trial=cfg.trials
                nC1=0;
                nC2=0;
                nT=nT+1;
                clear psd X1 X2
                for channelFrom=1:64%find(Coh_elec_all_ind)%[find(elec_post_ind) 22:23 59 60]%[31]
                    nC1=nC1+1;
                    nC2=0;
                    for channelTo=1:64%find(Coh_elec_all_ind)%[find(elec_post_ind) 22:23 59 60] %find(elec_parietalMid_ind)%59:60%find(elec_parietalMid_ind)%[26]
                        nC2=nC2+1;
                        wave1=wavChanAll(:,trial,channelFrom);
                        wave2=wavChanAll(:,trial,channelTo);
                        
                        X1 = fftshift(fft(wave1));
                        X2 = fftshift(fft(wave2));
                        
%                         FourierPermChannel(subject,perm,channelTo,:)=X1;
                        
                        numsmp = length(wave1);
                        psd(:,1) = 2.*abs(X1).^2./(numsmp.^2);
                        psd(:,2) = 2.*abs(X2).^2./(numsmp.^2);
                        csd = 2.*(X1.*conj(X2)./(numsmp.^2));
                        
                        fd{subject}{perm}.csd(nT,:,nC1,nC2) = csd;
                        fd{subject}{perm}.psd(nT,:,:,nC1,nC2) = psd;
                    end
                end
                
            end
            
            %             fd{subject}{perm}             = ft_connectivityanalysis(cfg, freqPowNcsd{subject});
        else
            permidx(perm)=0;
            continue
        end
    end
    
    %     if nSub==1
    %         iCoh{feat}=nan(length(Subj),6,nC1,nC2);
    %     end
    %
    for perm=find(permidx)
        
        Features1=Perm(perm,:);
        Features2=[Features1 Features1.*2];
        hz_intermod=intermod(Features1',[0.5 20]);
        
        for ii=1:nC1
            for jj=1:nC2
                sumpsd=squeeze(sum(fd{subject}{perm}.psd(:,:,:,ii,jj)));
                sumcsd=squeeze(sum(fd{subject}{perm}.csd(:,:,ii,jj)));
                
                
                coh = imag(sumcsd' ./ sqrt(sumpsd(:,1) .* sumpsd(:,2)));
                fd{subject}{perm}.cohChan(:,ii,jj)=coh;
                AllcohChan(nSub,perm,:,ii,jj)=coh;
                
                for feat=1:6
                    FeatFreq_indx = abs(bsxfun(@minus,Features2(feat),f));
                    [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
                    whichFreqFeat= [(index_freqToFeatures-1):(index_freqToFeatures+1)];
                    
                    
                    iCohFeat{feat}(nSub,perm,ii,jj)=AllcohChan(nSub,perm,whichFreqFeat(2),ii,jj);
                    iCohFeatBOOT{feat}(nSub,perm,ii,jj)=AllcohChan(nSub,perm,ceil(16000*rand(1)),ii,jj);
                end
                
                for interfeat=1:7
                    FeatFreq_indx = abs(bsxfun(@minus,hz_intermod(interfeat),f));
                    [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
                    whichFreqFeat= [(index_freqToFeatures-1):(index_freqToFeatures+1)];
                    
                    iCohInter{interfeat}(nSub,perm,ii,jj)=AllcohChan(nSub,perm,whichFreqFeat(2),ii,jj);
                    iCohInterBOOT{interfeat}(nSub,perm,ii,jj)=AllcohChan(nSub,perm,ceil(16000*rand(1)),ii,jj);
                end
                
                
                
            end
        end
        
    end
    %     sumcsd=squeeze(sum(fd{subject}{3}.csd));
    %     figure, plot(f,sumcsd)
    
    
    
    
    % % for perm=1
    % trials=find(TrialTagg{subject}==perm);%'all'
    % trials(trials>max(size(data_postica.trial)))=[];
    % end
    %
    % for trial=trials
    %     n=n+1;
    %     EEGtimeALL(2).data(:,:,n)= raw_data.trial{trial}(:,1:15360);  % or 17152
    % end
    %           for trial=1:max(size(raw_data.trial))
    %            n=n+1;
    %           EEGtimeALL(1).data(:,:,n)= raw_data.trial{trial}(:,1:340);  % or 17152
    %           end
    
    
    %     end
    
end

%%
% % % % % % load('iCoh_Conditions.mat') not the latest
% COND2=load('iCohCond2.mat');
% COND1=load('iCohCond1.mat');
% TO SAVE : which electrodes : Coh_elec_all_ind
Coh_elec_all={'T7','T8', 'POz', 'O1', 'Oz', 'O2','PO5', 'PO6', 'TP7', 'TP8','PO7','PO8','CP5','CP6','FT7','FT8','AF3','AF4'};
posterior_elec= Coh_elec_all;%{'T7','T8', 'POz', 'O1', 'Oz', 'O2','PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8','CP5','CP6','FT7','FT8'};
%  toChann={'CP5','CP1','CP2','CP6','T7','T8','FT7','FT8'};
toChann=posterior_elec;
legLabels=posterior_elec;
% for to=5
counter=0;
for feat=[1:2]
    counter=counter+1;
    %     iCohFeatCond2{feat}(isnan(iCohFeatCond2{feat}))=0;
    
    % CohMatSubjects{feat}=squeeze(nanmean((abs(iCohFeat{feat}(:,:,:,:))),2));
    % CohMatSubjectsCond2{feat}=squeeze(nanmean((abs(iCohFeatCond2{feat}(:,:,:,:))),2));
    %
    %CohMatAll{feat}=tril(squeeze(nanmean((abs(CohMatSubjectsCond2{feat})),1)));
    % % cohMat=squeeze(nanmean(squeeze(nanmean((abs(iCohInter{feat}(:,:,:,:))),1))))-squeeze(nanmean(squeeze(nanmean((abs(iCohInterCond2{feat}(:,:,:,:))),1))));
    % %  Diff2=squeeze(nanmean(squeeze(nanmean(((iCohFeat{feat}(:,:,:,:))),1))))'-squeeze(nanmean(squeeze(nanmean(((iCohFeat{2}(:,:,:,:))),1))))';
    % % subplot(1,2,counter),imagesc(cohMat'),colorbar
    %  subplot(1,3,counter),imagesc(CohMatAll{feat}),colorbar
    % subplot(1,3,feat),imagesc(Diff2,[-.25 .25]),colorbar
    Tz=(abs(COND1.iCohInter{feat})-mean(abs(COND1.iCohInterBOOT{feat}(:))))/std(abs(COND1.iCohInterBOOT{feat}(:)));
    % figure, bar(nanmean(squeeze(nanmean((abs(iCohFeat{feat}(:,:,:,to))),1)))),title(sprintf('iCoherence of occipital elec with %s',toChann{to})),axis([0 15 0 .5])
    % figure, imagesc(squeeze(nanmean(squeeze(nanmean((abs(iCohFeat{feat}(:,:,:,:))),1))))'),colorbar
    % figure, imagesc((Diff2+Diff1)/sqrt(2)>0.10),colorbar
    % figure, imagesc(Diff2),colorbar
    figure, imagesc(tril((squeeze(nanmean(nanmean(abs(Tz))))))),colorbar
    set(gca,'XTick',1:18,'XTickLabel',posterior_elec),xlabel('From channel'), title(sprintf('feature %d',feat'))
    set(gca,'YTick',1:18,'YTickLabel',posterior_elec),ylabel('To channel')
    % hold on
    % imagesc(squeeze(nanmean(squeeze(nanmea22n((abs(iCohFeat{feat}(:,:,:,:))),1))))'>.3)
    
end
% end
%%
for interfeat=1:2
    iCohInter{interfeat}(iCohInter{interfeat}==0)=nan;
    figure, bar(nanmean(squeeze(nanmean((iCohInter{interfeat}(:,:,:,2)),1)))),title(sprintf('iCoherence of occipital elec with %s',data_postica.label{18}))
    set(gca,'XTick',1:18,'XTickLabel',posterior_elec)
end


%%
% CHANLOCS = readeetraklocs('/Applications/fieldtrip-20151020/template/electrode/standard_1020.elc');
elec_montage = ft_read_sens('/Users/SimonFaghel-Soubeyrand/Desktop/DATA_EEG/standard_waveguard64.elc');
for to=1
    FtsevenCohLeft=abs(iCohFeatCond2{1}(:,:,:,1:(end-4)));
    FtsevenCohRight=abs(iCohFeatCond2{2}(:,:,:,1:(end-4)));
    FtsevenCohMouth=abs(iCohFeatCond2{3}(:,:,:,1:(end-4)));
    [h,p1{to}]=ttest2(FtsevenCohLeft(:),FtsevenCohRight(:))
    [h,p2{to}]=ttest2(FtsevenCohMouth(:),FtsevenCohLeft(:))
    [h,p3{to}]=ttest2(FtsevenCohMouth(:),FtsevenCohRight(:))
    
    figure, bar([nanmean(FtsevenCohLeft(:)) nanmean(FtsevenCohRight(:)) nanmean(FtsevenCohMouth(:))]),title(sprintf('left eye, right eye, mouth, occ to %s', 'occ')),axis([0 4 0 .5])
    sigstar({[1,2] [1,3] [2,3]},[p1{to} p2{to} p3{to}])
    
end

%%
%
for to=1:18
    FtsevenCohLeft=(abs(iCohInterCond1{1}(:,:,:,to)+abs(iCohInterCond1{2}(:,:,:,to))))/2;%+abs(iCohInter{2}(:,:,:,:));
    FtsevenCohRight=(abs(iCohInterCond1{3}(:,:,:,to)+abs(iCohInterCond1{4}(:,:,:,to))))/2;
    FtsevenCohMouth=(abs(iCohInterCond1{5}(:,:,:,to)+abs(iCohInterCond1{6}(:,:,:,to))))/2;
    [h,p1{to}]=ttest2(FtsevenCohLeft(:),FtsevenCohRight(:))
    [h,p2{to}]=ttest2(FtsevenCohMouth(:),FtsevenCohLeft(:))
    [h,p3{to}]=ttest2(FtsevenCohMouth(:),FtsevenCohRight(:))
    
    figure, bar([nanmean(FtsevenCohLeft(:)) nanmean(FtsevenCohRight(:)) nanmean(FtsevenCohMouth(:))]),title(sprintf('left eye, right eye, mouth, occ to %s', posterior_elec{to})),axis([0 4 0 .4])%,axis([0 4 -0.1 .1])
    sigstar({[1,2] [1,3] [2,3]},[p1{to} p2{to} p3{to}])
    
end

% ds.chanPairs* (required) - N x 2 matrix, with N being the number of
% connected channel pairs. For example, ds.chanPairs = [7, 12; 13 20];
% specifies two connected channel pairs (7, 12) and (13, 20).
displayStruc.channPairs=[];
for jj=1:18
    for ii=1:18,
        ds.channPairs=[displayStruc.channPairs;ii, jj];
    end
end
% * *ds.connectStrength* (optional) - N x 1 matrix, a vector specifying
% connection strengths. If unspecified, then the connections will be
% rendered in a color at the center of the current colormap.
ds.connectStrength=squeeze(nanmean(nanmean(abs(iCohInterCond1{1}(:,:,:,:)))));


figure, ft_plot_topo(elec_montage.chanpos(Coh_elec_all_ind==1,:,:,:),squeeze(nanmean(nanmean((iCohInter{1}(:,:,:,1)),1)))),colorbar
% ft_topoplotCC(

test=data_postica;
test.cohspctrm=ds.connectStrength;
test.avg.cohspctrm=mean(ds.connectStrength(:));
cfg.channel=Coh_elec_all;
cfg.layout=elec_montage;
cfg.elec=elec_montage;
figure, ft_multiplotCC(cfg,test)

figure,ft_plot_topo(squeeze(elec_montage.chanpos(Coh_elec_all_ind==1,1)),squeeze(elec_montage.chanpos(Coh_elec_all_ind==1,2)),squeeze(nanmean(nanmean((iCohInter{1}(:,:,:,1)),1))))%,'mask',lay.mask)

% topoplot_connect(ds,elec_montage.chanpos)

%% COMPARE CONDITIONS

% posterior_elec= {'T7','T8', 'POz', 'O1', 'Oz', 'O2','PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'};
% toChann={'CP5','CP1','CP2','CP6','T7','T8','FT7','FT8'};
% toChann=posterior_elec;
% legLabels=posterior_elec;

% figure,
for feat=[1:6]
    iCohInterCond2{feat}((iCohInterCond2{feat})==0)=NaN;
    % iCohInterCond2{feat}((iCohInterCond2{feat})==0)=NaN;
    
    % % ShowCond2{feat}=squeeze(nanmean(squeeze(nanmean(abs(iCohInterCond2{feat}),1))))';
    % ShowCond1{feat}=squeeze(nanmean(squeeze(nanmean(abs(iCohInterCond1{feat}),1))))'-squeeze(nanmean(squeeze(nanmean(abs(iCohInterCond1{3}),1))))';
    % % Diff{feat}=(iCohInterCond1{feat}(:,:,:,:))-(iCohInterCond2{feat}(:,:,:,:));
    %
    % % Diff{feat}=squeeze(nanmean(squeeze(nanmean(Diff{feat},1))))';
    %
    %
    % subplot(1,2,feat),imagesc(ShowCond1{feat}),colorbar
    % set(gca,'XTick',1:14,'XTickLabel',posterior_elec),xlabel('From channel')
    % set(gca,'YTick',1:14,'YTickLabel',toChann),ylabel('To channel')
    % subplot(1,2,2),imagesc(ShowCond2{feat},[.3 .45]),colorbar
    % set(gca,'XTick',1:14,'XTickLabel',posterior_elec),xlabel('From channel')
    % set(gca,'YTick',1:8,'YTickLabel',toChann),ylabel('To channel')
    
    % figure, imagesc(Diff{feat},[0 .2]),colorbar
    % set(gca,'XTick',1:14,'XTickLabel',posterior_elec),xlabel('From channel')
    % set(gca,'YTick',1:8,'YTickLabel',toChann),ylabel('To channel')
    
    % hold on
    % imagesc(squeeze(nanmean(squeeze(nanmean((abs(iCohFeat{feat}(:,:,:,:))),1))))'>.3)
    
end
% end

%%

iCohInterCond1All{1}=squeeze(nanmean(abs(iCohInterCond1{1})+abs(iCohInterCond1{2})/2,2));
iCohInterCond1All{2}=squeeze(nanmean(abs(iCohInterCond1{3})+abs(iCohInterCond1{4})/2,2));
iCohInterCond1All{3}=squeeze(nanmean(abs(iCohInterCond1{5})+abs(iCohInterCond1{6})/2,2));
iCohInterCond1All{4}=squeeze(nanmean(abs(iCohInterCond1{7}),2));

iCohInterCond2All{1}=squeeze(nanmean(abs(iCohInterCond2{1})+abs(iCohInterCond2{2})/2,2));
iCohInterCond2All{2}=squeeze(nanmean(abs(iCohInterCond2{3})+abs(iCohInterCond2{4})/2,2));
iCohInterCond2All{3}=squeeze(nanmean(abs(iCohInterCond2{5})+abs(iCohInterCond2{6})/2,2));
iCohInterCond2All{4}=squeeze(nanmean(abs(iCohInterCond2{7}),2));
% iCohInterCond2All{3}


%%
for subject=1:8
    figure,
    for ii=1
        iCohInterCond1All{ii}(isnan(iCohInterCond1All{ii}))=0;
        iCohInterCond2All{ii}(isnan(iCohInterCond2All{ii}))=0;
        ShowCond1{ii}=squeeze(nanmean((iCohInterCond1All{ii}),1))';
        ShowCond2{ii}=squeeze(nanmean((iCohInterCond2All{ii}),1))';
        Diff=tril((squeeze(iCohInterCond1All{ii}(subject,:,:))-(squeeze(iCohInterCond2All{ii}(subject,:,:)))));%squeeze(iCohInterCond1All{ii}(subject,:,:))
        Diff=Diff;%-tril((squeeze(iCohInterCond1All{2}(subject,:,:))+squeeze(iCohInterCond2All{2}(subject,:,:)))/2);
        subplot(1,1,ii),imagesc(Diff),colorbar, title(sprintf('mean iCoh between condition : %.2f',mean(Diff(:))))
        
        set(gca,'XTick',1:14,'XTickLabel',posterior_elec),xlabel('From channel')
        set(gca,'YTick',1:14,'YTickLabel',toChann),ylabel('To channel')
        % subplot(1,2,2),imagesc(squeeze(iCohInterCond2All{ii}(subject,:,:)),[0 .70]),colorbar
        set(gca,'XTick',1:14,'XTickLabel',posterior_elec),xlabel('From channel')
        set(gca,'YTick',1:14,'YTickLabel',toChann),ylabel('To channel')
        
    end
end
% iCohRightEMouthCond2All(isnan(iCohRightEMouthCond2All))=0;
% iCohLeftEMouthCond2All(isnan(iCohLeftEMouthCond2All))=0;
% iCohBothEyesCond2All(isnan(iCohBothEyesCond2All))=0;









%%
figure, cfg.zlim      = [];cfg.parameter='psispctrm';
cfg.channel='all'%{'Oz','O2','PO6','CP6','FT7'};
cfg.refchannel={'Cz'};

cfg.xlim=[hz_intermod(2) hz_intermod(2)+0.15]
% ft_connectivityplot(cfg,fd{subject}{1})
ft_topoplotER(cfg,fd{subject}{2})

Features1=Perm(1,:);
Features2=[Features1 Features1.*2];
hz_intermod=intermod(Features1',[0.5 20]);


FeatFreq_indx = abs(bsxfun(@minus,Features2(feat),fd{subject}{1}.freq));
[~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
whichFreqFeat= [(index_freqToFeatures-1):(index_freqToFeatures+1)]

%   save('wpliStructure.mat','fd','TrialTagg','-v7.3' )
%%
% newcrossf({ mean(EEGtimeALL(1).data(29:31,:,:),1)  mean(EEGtimeALL(2).data(29:31,:,:),1)},...
% {(EEGtimeALL(1).data(26,:,:))  (EEGtimeALL(2).data(26,:,:))},340,[-500 1000],256,0,'maxfreq',20,'type','phasecoher')%,'base 'padratio' 'elocs','BioSemi64.loc'
% %%
% newcrossf({ (EEGtimeALL(1).data(18,:,:))  (EEGtimeALL(2).data(18,:,:))},...
% {(EEGtimeALL(1).data(14,:,:))  (EEGtimeALL(2).data(14,:,:))},128,[0 500],256,0,'maxfreq',100)%,'base 'padratio' 'elocs','BioSemi64.loc'
%

% segmentLength=round((numel(mean(EEGtimeALL(2).data(18,:,:),3))/20))
% figure, spectrogram(mean(EEGtimeALL(2).data(18,:,:),3),segmentLength,[],[],256,'yaxis')
%
% figure, mscohere(mean(EEGtimeALL(2).data(18,:,:),3),mean(EEGtimeALL(2).data(30,:,:),3))
%  ft_connectivityplot(cfg, fd{subject}{perm});

%%
freqPowNcsd{8}{perm}.fourierspctrm=freqPowNcsd{8}{perm}.fourierspctrm(1,:,1:1283);
freqPowNcsd{8}{perm}.freq=freqPowNcsd{8}{perm}.freq(1,1:1283);

cfg.parameter='fourierspctrm';
APPFREQALL{perm} =ft_appendfreq(cfg,freqPowNcsd{8}{perm},freqPowNcsd{12}{perm},freqPowNcsd{14}{perm})
grandcfg.parameter='cohspctrm';


cohcfg.keeptrials='no';
for subject=Subj_cond1
    [FreqDescr{subject}] = ft_freqdescriptives(cohcfg,append{subject})
end

[grandavg] = ft_freqgrandaverage(cohcfg,FreqDescr{14}, FreqDescr{12},FreqDescr{16},FreqDescr{18},FreqDescr{20}...
    ,FreqDescr{22},FreqDescr{24})


cfg.keeptrials='no';
grandavgWPLI=ft_freqgrandaverage(cfg,append{subject});

% n=0;
% fd_all=fd{subject};
% fd_all.cohspctrm=zeros(size(fd{subject}.cohspctrm));
% for subject=Subj_cond1
%     n=n+1;
% if~n==1
% if ~size(fd{subject}.cohspctrm,3)>size(fd{Subj_cond1(n-1)}.cohspctrm,3)
% fd_all.cohspctrm=fd_all.cohspctrm+fd{subject}.cohspctrm;
% else
%   fd_all.cohspctrm=fd_all.cohspctrm+fd{subject}.cohspctrm(:,:,size(fd{Subj_cond1(n-1)}.cohspctrm,3));
%  end
% end
%
% end
% fd_all.cohspctrm=fd_all.cohspctrm;

cfg.layout=lay;
% eloc = readlocs('BioSemi64.loc');
%%

load 'coherenceStructure.mat'

whichConnect=2;
Connect_index_cell={'plvspctrm','cohspctrm','wplispctrm','powcorrspctrm','crsspctrm','psispctrm'}
plotCFG_index={'plvspctrm','cohspctrm','powspctrm','crsspctrm'}

cfg_index= {'plv','coh','powcorr','powcorr_ortho','csd'};

Subj_cond1=[8 12 16 18 20 22 24 26 28] % Cond 1 only
Subj_cond2=[9 13 17 19 21 23 25 27 29] % Cond 2 only
load('quickcap64.mat') % load 2D layout of electrodes
freqcfg.channel=ft_channelselection({'all','-EOG','-M1','-M2', '-Fp1','-Fp2','-Fpz'},fd{12}{1}.label)
% Perm =
%
%     6.3150    8.0000    3.5290
%     6.3150    3.5290    8.0000
%     8.0000    6.3150    3.5290
%     8.0000    3.5290    6.3150
%     3.5290    8.0000    6.3150
%     3.5290    6.3150    8.0000

feat= [3.5290    8.0000    6.3150];
Perm = perms(feat);

Perm_Taggs{1}={'Mouth','Left eye','Right eye'};
Perm_Taggs{2}={'Right eye','Left eye','Mouth'};
Perm_Taggs{3}={'Mouth','Right eye','Left eye'};
Perm_Taggs{4}={'Right eye','Mouth','Left eye'};
Perm_Taggs{5}={'Left eye','Mouth','Right eye'};
Perm_Taggs{6}={'Left eye','Right eye','Mouth'};

clear Feature_coherence
clear Feature_Structure


cfg.layout           = lay;
cfg.channel          = freqcfg.channel;
plotcfg.parameter= Connect_index_cell{whichConnect};
clear InterFeature_coherence Feature_Structure Feature_coherence InterFeature_coherence_Structure InterFeature_coherence_Structure
cond=1;

nSub=0;
for subject=[Subj_cond2(1:end)]
    nSub=nSub+1;
    
    for hh=1:6
        idxPerm(hh)=~isempty(fd{subject}{hh});
    end
    
    for feat=1:6,Feature_coherence{nSub}{feat}=zeros(66,66);end
    for interfeat=1:7,InterFeature_coherence{nSub}{interfeat}=zeros(66,66);end
    
    
    
    for Permutation=find(idxPerm)
        
        Features1=Perm(Permutation,:);
        Features2=[Features1 Features1.*2];
        hz_intermod=intermod(Features1',[0.5 20]);
        
        
        for feat=1:6
            
            FeatFreq_indx = abs(bsxfun(@minus,Features2(feat),fd{subject}{Permutation}.freq));
            [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
            whichFreqFeat= [(index_freqToFeatures-1):(index_freqToFeatures+1)];
            
            Feature_coherence{nSub}{feat}(:,:)=Feature_coherence{nSub}{feat}+eval(sprintf('mean((fd{subject}{Permutation}.%s(:,:,whichFreqFeat(2))),3)',Connect_index_cell{whichConnect}));
        end
        
        for interfeat=1:7
            FeatFreq_indx = abs(bsxfun(@minus,hz_intermod(interfeat),fd{subject}{Permutation}.freq));
            [~,index_freqToFeatures] = min(FeatFreq_indx(:,1:size(FeatFreq_indx,2)));
            whichFreqInterFeat= [(index_freqToFeatures-1):(index_freqToFeatures+1)];
            
            
            InterFeature_coherence{nSub}{interfeat}(:,:)=InterFeature_coherence{nSub}{interfeat}+eval(sprintf('mean((fd{subject}{Permutation}.%s(:,:,whichFreqInterFeat(2))),3)',Connect_index_cell{whichConnect}));
        end
    end
    for feat=1:6
        Feature_coherence{nSub}{feat}(:,:)=Feature_coherence{nSub}{feat}(:,:)./sum(idxPerm);
        if ~isempty(fd{subject}{1})
            Feature_Structure{nSub}{feat}=fd{subject}{1};
        else
            Feature_Structure{nSub}{feat}=fd{subject}{2};
        end
        eval(sprintf('SS=max(size(Feature_Structure{nSub}{feat}.%s(:,:,1:end)));',Connect_index_cell{whichConnect}));
        eval(sprintf('Feature_Structure{nSub}{feat}.%s(:,:,1:SS)=repmat(Feature_coherence{nSub}{feat}(:,:),1,1,SS);',Connect_index_cell{whichConnect}));
    end
    
    for interfeat=1:7
        InterFeature_coherence{nSub}{interfeat}(:,:)= InterFeature_coherence{nSub}{interfeat}(:,:)./sum(idxPerm);
        if ~isempty(fd{subject}{1})
            InterFeature_coherence_Structure{nSub}{interfeat}=fd{subject}{1};
        else
            InterFeature_coherence_Structure{nSub}{interfeat}=fd{subject}{2};
        end
        eval(sprintf('SS=max(size(InterFeature_coherence_Structure{nSub}{interfeat}.%s(:,:,1:end)));',Connect_index_cell{whichConnect}));
        eval(sprintf('InterFeature_coherence_Structure{nSub}{interfeat}.%s(:,:,1:end)=repmat(InterFeature_coherence{nSub}{interfeat}(:,:),1,1,SS);',Connect_index_cell{whichConnect}));
    end
    % fd{subject}{Permutation}
    
    
    for ii=1:7,eval(sprintf('CohInterFeature{cond}{ii}.%s(nSub,:,:)=InterFeature_coherence{nSub}{ii}(:,:);',Connect_index_cell{whichConnect})); end
    
    for ii=1:6,eval(sprintf('CohFeature{cond}{ii}.%s(nSub,:,:)=Feature_coherence{nSub}{ii}(:,:);',Connect_index_cell{whichConnect})); end
    % cfg.alphaparam=0.6;
    % cfg.arrowhead        ='none';
    cfg.channel          =freqcfg.channel;%{'T8','P8','P4','P6', 'P2','PO6','TP8','PO8','T7','P7','P3','P5', 'P1','PO5','TP7','PO7'}%posterior_elec;%'all' ;
    %
    %  cfg.foi=Features(feat);
    % for feat=1:3
    %      ft_topoplotCC(cfg,Feature_Structure{nSub}{feat}), title(sprintf('%s',Perm_Taggs{6}{feat})),colorbar
    % end
end
%%
% Out-degree: Number of causal outflow connections from a node in the network to any other node,
% In-degree: Number of causal in-flow connections to a node in the network from any other node, and
cfg.channel          =freqcfg.channel;
SS=3308;
SS=200;
% Feature Coherence
clear Coherence_FeatureALL

for feat=1:3
    Coherence_FeatureALL{feat}=fd{9}{1};
    eval(sprintf('Coherence_FeatureALL{feat}.%s=zeros(size(Feature_Structure{1}{feat}.%s(:,:,1:SS)));',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}))
    Coherence_FeatureALL{feat}.freq(:)=[];
    Coherence_FeatureALL{feat}.freq=fd{12}{1}.freq(1:SS);
end

whichSub=[1:9]; % maybe throw sub 3?

clear cond1 cond2 zDiff
for subject=whichSub
    % figure,
    for feat=1:3
        eval(sprintf('cond1=(squeeze(CohFeature{1}{feat}.%s(subject,:,:)));',Connect_index_cell{whichConnect}))
        bootCoher= cond1;
        bootIndex=ceil(size(bootCoher,2)*rand(size(bootCoher)));
        bootStrap=cond1(bootIndex);
        zCond1=(cond1-mean(bootStrap(~isnan(bootStrap))))/std(bootStrap(~isnan(bootStrap)));
        
        eval(sprintf('cond2=(squeeze(CohFeature{2}{feat}.%s(subject,:,:)));',Connect_index_cell{whichConnect}))
        bootCoher= cond2;
        bootIndex=ceil(size(bootCoher,2)*rand(size(bootCoher)));
        bootStrap=cond2(bootIndex);
        zCond2=(cond2-mean(bootStrap(~isnan(bootStrap))))/std(bootStrap(~isnan(bootStrap)));
        
        zDiff=zCond2-zCond1/sqrt(2);
        zSum=zCond2+zCond1/sqrt(2);
        zCond1(isnan(zCond1))=mean(bootStrap(~isnan(bootStrap)));
        %         pCond1 = (1-normcdf(zCond1,mean(zCond1(:)),std(zCond1(:))));
        
        eval(sprintf('Coherence_FeatureALL{feat}.%s(:,:,1:SS)=Coherence_FeatureALL{feat}.%s+repmat(cond2,1,1,SS);',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}));% InterFeature_coherence_Structure{subject}{interfeat}.cohspctrm(:,:,1:SS);
        cfg.xlim=[1 1.1];
        
    end
end

for feat=1:3, eval(sprintf('Coherence_FeatureALL{feat}.%s(:,:,1:end)=Coherence_FeatureALL{feat}.%s(:,:,1:end)./(length(whichSub));',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect})); end
%
cfg.refchannel       ={'Oz'}%{'P5','P6','CP5','CP6'}%{'PO6','PO8','PO4'}%{'Oz','O1','O2'}%{'PO6','PO8','PO4'}%{'Oz','O1','O2'};%{'Cz','FCz'}%{'P5','P6','CP5','CP6'}%{'PO6','PO8'}%{'FT7','FT8'}%{'CPz'}%{'T8','PO8'}; % T8 PO8
cfg.colorbar         ='South';
cfg.directionality='inflow';
%  cfg.zlim=[0 .4];
cfg.zlim='maxmin';
cfg.xlim=[3.52 3.53]

cfg.parameter        = Connect_index_cell{whichConnect};%Connect_index_cell{whichConnect};  %'
%
figure,
for feat=1:3
    subplot(1,3,feat),ft_topoplotER(cfg,Coherence_FeatureALL{feat}), title(sprintf('%s',Perm_Taggs{6}{feat}))
end

%
% for feat=1:3, reshape(find(Coherence_FeatureALL{feat}.cohspctrm(:,:,:)<-3.5),[66 66 200]);end
% cfg.maskparameter=
figure,
cfg.parameter        = 'pvalues';%Connect_index_cell{whichConnect};  %'
cfg.zlim=[0 1]
cfg.highlightchannel={'PO6','TP8','CP6','Oz','TP7','CP5'};
cfg.highlight ='label';
cfg.highlightsymbol    = 'o';
cfg.highlightcolor= [1 0 0];
cfg.highlightsize      = 10;%highlight marker size (default = 6)
cfg.highlightfontsize  = 12;%highlight marker size (default = 8)

cfg.zlim=[0 .1]
for feat=1:3
    Zmat=squeeze(Coherence_FeatureALL{feat}.cohspctrm(:,:,1));
    Zmat(isnan(Zmat))=mean(Zmat(~isnan(Zmat)));
    pMat= (1-normcdf(Zmat,mean(Zmat(:)),std(Zmat(:))));
    
    Coherence_FeatureALL{feat}.pvalues(:,:,1:SS)=repmat(pMat,1,1,SS);
    subplot(1,3,feat),ft_topoplotER(cfg,Coherence_FeatureALL{feat}), title(sprintf('%s',Perm_Taggs{6}{feat}))
end

%%
figure
for row=14:18
    for col=14:18
        subplot(4,4,(row-1)*4+col);
        plot(fd{8}{5}.freq, squeeze(fd{8}{5}.plvspctrm(row,col,:)))
        ylim([0 1])
    end
    
end


%% INTERMODULATE FEATURE COHERENCE

cfg.channel          =freqcfg.channel;
SS=3308;
SS=200;
% Feature Coherence
clear Coherence_FeatureALL

for interfeat=1:7
    Coherence_FeatureALL{interfeat}=fd{9}{1};
    eval(sprintf('Coherence_FeatureALL{interfeat}.%s=zeros(size(Feature_Structure{1}{1}.%s(:,:,1:SS)));',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}))
    Coherence_FeatureALL{interfeat}.freq(:)=[];
    Coherence_FeatureALL{interfeat}.freq=fd{12}{1}.freq(1:SS);
end

whichSub=[1:9];
clear cond1 cond2 zDiff
for subject=whichSub
    % figure,
    for interfeat=1:7
        eval(sprintf('cond1=squeeze(CohInterFeature{1}{%d}.%s(%d,:,:));',interfeat,Connect_index_cell{whichConnect},subject))
        bootCoher= cond1;
        bootIndex=ceil(size(bootCoher,2)*rand(size(bootCoher)));
        bootStrap=cond1(bootIndex);
        zCond1=(cond1-mean(bootStrap(~isnan(bootStrap))))/std(bootStrap(~isnan(bootStrap)));
        
        eval(sprintf('cond2=squeeze(CohInterFeature{2}{%d}.%s(%d,:,:));',interfeat,Connect_index_cell{whichConnect},subject))
        bootCoher= cond2;
        bootIndex=ceil(size(bootCoher,2)*rand(size(bootCoher)));
        bootStrap=cond2(bootIndex);
        zCond2=(cond2-mean(bootStrap(~isnan(bootStrap))))/std(bootStrap(~isnan(bootStrap)));
        
        zDiff=zCond2-zCond1/sqrt(2);
        zSum=zCond2+zCond1/sqrt(2);
        
        eval(sprintf('Coherence_FeatureALL{%d}.%s(:,:,1:SS)=Coherence_FeatureALL{%d}.%s+repmat(cond1-cond2,1,1,SS);',interfeat,Connect_index_cell{whichConnect},interfeat,Connect_index_cell{whichConnect}));% InterFeature_coherence_Structure{subject}{interfeat}.cohspctrm(:,:,1:SS);
        cfg.xlim=[1 1.1];
        
    end
end

for interfeat=1:7, eval(sprintf('Coherence_FeatureALL{interfeat}.%s(:,:,1:end)=Coherence_FeatureALL{interfeat}.%s(:,:,1:end)./(length(whichSub));',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}));end


%
eval(sprintf('Coherence_FeatureALL{1}.%s(:,:,1:end)=(Coherence_FeatureALL{1}.%s(:,:,1:end)+Coherence_FeatureALL{2}.%s(:,:,1:end))./(2);',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}));
eval(sprintf('Coherence_FeatureALL{2}.%s(:,:,1:end)=(Coherence_FeatureALL{3}.%s(:,:,1:end)+Coherence_FeatureALL{4}.%s(:,:,1:end))./(2);',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}));
eval(sprintf('Coherence_FeatureALL{3}.%s(:,:,1:end)=(Coherence_FeatureALL{5}.%s(:,:,1:end)+Coherence_FeatureALL{6}.%s(:,:,1:end))./(2);',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}));
eval(sprintf('Coherence_FeatureALL{4}.%s(:,:,1:end)=Coherence_FeatureALL{7}.%s(:,:,1:end);',Connect_index_cell{whichConnect},Connect_index_cell{whichConnect}));
%
% clear Diff
% for ii=1:4
%  Diff{ii}=Coherence_FeatureALL{ii};
%  Diff{ii}.cohspctrm=(Coherence_FeatureALLCOND1{ii}.cohspctrm(:,:,:))-(Coherence_FeatureALLCOND2{ii}.cohspctrm(:,:,:));
% end
%

% figure, imagesc(squeeze(mean(fd{12}{5}.plvspctrm(find(elec_post_ind),find(elec_post_ind),:),1))),colorbar
figure,
cfg.layout=lay;
cfg.parameter=Connect_index_cell{whichConnect};
cfg.refchannel       ={'CP5','CP6'}%{'Oz','O1','O2'}%{'CP5','CP6'}%{'POz','Oz','O1','O2'}%{'P5','P6','CP5','CP6'}%{'Oz','O1','O2'}%{'P5','P6','PO8','PO6','PO4','PO7','PO5','PO3','POz'}%{'Oz'};%{'PO8','PO6','PO4'}%{'O1','Oz','O2'}%{'P5','P6','CP5','CP6'}%{'Oz'}%{'CPz'}%{'T8','PO8'}; % T8 PO8
cfg.colorbar         ='South';
cfg.directionality='inflow';
%  cfg.zlim=[-4 4];
cfg.zlim='maxmin';
for interfeat=1:4
    
    subplot(2,2,interfeat),ft_topoplotER(cfg,Coherence_FeatureALL{interfeat})
    % subplot(2,2,interfeat),ft_topoplotER(cfg,Diff{interfeat})
end

%%
% for feat=1:3, reshape(find(Coherence_FeatureALL{feat}.cohspctrm(:,:,:)<-3.5),[66 66 200]);end
% cfg.maskparameter=
figure,
cfg.parameter        = 'pvalues';%Connect_index_cell{whichConnect};  %'
cfg.zlim=[0 .1]
cfg.highlight ='label';
cfg.highlightsymbol    = 'o';
cfg.highlightcolor= [1 0 0];
cfg.highlightsize      = 8;%highlight marker size (default = 6)
cfg.highlightfontsize  = 12;%highlight marker size (default = 8)
cfg.colorbar         ='no';
for feat=1:7
    Zmat=squeeze(Coherence_FeatureALL{feat}.wplispctrm(:,:,1));
    Zmat(isnan(Zmat))=mean(Zmat(~isnan(Zmat)));
    pMat= (1-normcdf(Zmat,mean(Zmat(:)),std(Zmat(:))));
    
    Coherence_FeatureALL{feat}.pvalues(:,:,1:SS)=repmat(pMat,1,1,SS);
    subplot(2,2,feat),ft_topoplotER(cfg,Coherence_FeatureALL{feat})
end


%%
cfg.refchannel       ={'P5','P6','PO8','PO6','PO4','PO7','PO5','PO3','POz'}%{'P5','P6','CP5','CP6'}%{'PO6','PO8'}%{'P5','P6','CP5','CP6'}%{'PO6','PO8'}%%{'CP5','CP6'}%'FT7','FT8',{'CP5','CP6'}%{'CPz'}%{'T8','PO8'}; % T8 PO8 Real good IM with CP6-CP5 : STS??, Superior temporal G?GSM or FT7,FT8 from L. Koessler et al. / NeuroImage 46 (2009) 64?72
cfg.colorbar         ='no';%'South';
cfg.directionality='inflow';
cfg.zlim='maxmin';
figure
for interfeat=1:7
    
    subplot(2,4,interfeat),ft_topoplotER(cfg,Coherence_FeatureALL{interfeat})
    % subplot(2,2,interfeat),ft_topoplotER(cfg,Diff{interfeat})
end

%%

cfg.widthparam='12';
cfg.alphaparam=1;
cfg.arrowhead        ='none';%{'T8','P8','P4','P6', 'P2','PO6','TP8','PO8'}
cfg.channel          = Coh_elec_all(1:12);%{'T8','P8','P4','P6', 'P2','PO6','TP8','PO8','T7','P7','P3','P5', 'P1','PO5','TP7','PO7','Oz','POz','O1','O2'}%posterior_elec;%'all' ;
cfg.foi=1.1;%Features(feat);
cfg.directionality=[];
cfg.arrowoffset = 1;
for feat=1:3
    % ft_topoplotCC(cfg,Coherence_FeatureALL{feat}), title(sprintf('%s',Perm_Taggs{6}{feat})),colorbar
    ft_topoplotCC(cfg,Coherence_FeatureALL{feat}),colorbar%, title(sprintf('%s',Perm_Taggs{6}{feat})),colorbar
    
end
%%



% figure,
% for feat=1
% % Coherence_FeatureALL{feat}.cohspctrm(:,:,1:100)=Coherence_FeatureALL{feat}.cohspctrm(:,:,1:100)./nSub;
% ft_connectivityplot(cfg,Coherence_FeatureALL{feat})%, title(sprintf('%s',Perm_Taggs{6}{feat}))
% end
%
% figure,
% for feat=1
% % Coherence_FeatureALL{feat}.cohspctrm(:,:,1:100)=Coherence_FeatureALL{feat}.cohspctrm(:,:,1:100)./nSub;
% ft_multiplotCC(cfg, )%, title(sprintf('%s',Perm_Taggs{6}{feat}))
% end


% diff=Coherence_FeatureALL{1};
% diff.cohspctrm=Coherence_FeatureALL{3}.cohspctrm-(Coherence_FeatureALL{1}.cohspctrm+Coherence_FeatureALL{2}.cohspctrm)/2;
% ft_topoplotCC(cfg,diff), title('diff')
%%
for subject=Subj_cond1(1)
    cfg                  = [];
    cfg.parameter        = 'cohspctrm'%'powcorrspctrm';%plvspctrm'% 'cohspctrm';%'wpli_debiasedspctrm'%'powcorrspctrm'%'plvspctrm'%'cohspctrm';
    cfg.xlim             = [1 20];
    % cfg.foi            = [6.32];
    
    cfg.channel          = [1:64]%freqcfg.channel;
    cfg.refchannel       = {'PO6'}; % T8 PO8
    cfg.layout           = lay;
    cfg.showlabels       = 'yes';
    cfg.directionality = 'inflow';
    figure; ft_multiplotER(cfg, fd{subject})
    
end
%%
% elec = ft_read_sens('/Users/SimonFaghel-Soubeyrand/Documents/MATLAB/fieldtrip-master/template/electrode/standard_1020.elc')
subject=8
% posterior_elec_Left = {'T7','P7','P3','P5', 'P1','PO5','TP7','PO7'};
% posterior_elec_Right = {'T8','P8','P4','P6', 'P2','PO6','TP8','PO8'};
% posterior_elec_mid = {'Pz','POz','Oz','O2','O1','PO4','PO3'};
cfg.refchannel       = [];%{'Oz'}; % T8 PO8
cfg.layout           = lay;
cfg.channel          ={'T8','Oz','T7'} ;
cfg.arrowhead        ='none';
cfg.foi=6.32;
figure;ft_topoplotCC(cfg,fd{subject}{1}),colorbar



%% Same with intermodulate frequencies
hz_intermod=[4.4710   11.5290    1.6850   14.3150    2.7860    9.8440   17.8440];


for subject=Subj_cond1
    figure;
    for ii=4
        cfg.refchannel       = 'Oz'; % T8 PO8
        cfg.layout           = lay;
        switch (ii)
            case 1
                cfg.xlim             = [hz_intermod(1)-0.02 hz_intermod(1)+.02];
            case 2
                cfg.xlim             = [hz_intermod(3)-.02 hz_intermod(3)+.02];
            case 3
                cfg.xlim             = [hz_intermod(5)-.02 hz_intermod(5)+.02];
            case 4
                cfg.xlim             = [hz_intermod(end)-.02 hz_intermod(end)+.02];
        end
        cfg.channel          =[1:64];
        cfg.parameter        = 'cohspctrm';
        cfg.arrowhead        ='both';
        cfg.zlim             = 'maxmin';%[0 0.5];
        cfg.colorbar         ='South';
        cfg.foi=3;
        % subplot(1,1,ii),
        ft_topoplotER(cfg,fd{subject})
    end
end



%%

plotcfg.channel          = [1:64]%freqcfg.channel;
plotcfg.refchannel       = {'O2'}%{'O1','O2','Oz'}; % T8 PO8
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='East';
plotcfg.layout=lay;
plotcfg.parameter= 'cohspctrm';
plotcfg.xlim           = [6.3 6.32]%[3.50 3.53];%[6.3 6.32];8]  [3.5290]  [6.3150]
plotcfg.showlabels ='yes';

figure, ft_topoplotER(plotcfg,fd)

%%

plotcfg.channel          = [1:64]%freqcfg.channel;
plotcfg.refchannel       = {'Oz'}%{'F8','T8'}%{'O1','O2','Oz'};%{'O2'}% % T8 PO8
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='South';
plotcfg.layout=lay;
plotcfg.parameter= 'cohspctrm'%'plvspctrm'%'cohspctrm';
plotcfg.showlabels ='yes';



freq={3.5290,6.315, 8,4.47}

figure,
for frequency=1:4
    plotcfg.xlim           = [freq{frequency}-0.015 freq{frequency}+0.015];%[6.3 6.32];8]  [3.5290]  [6.3150]
    subplot(1,length(freq),frequency), ft_topoplotER(plotcfg,fd_all)
end
%%
intermodulaires=intermod([3.5290 6.3158 8],[0.5 20])
freq=intermodulaires




figure,
for frequency=1:4
    plotcfg.xlim           = [freq(frequency)-0.015 freq(frequency)+0.015];%[6.3 6.32];8]  [3.5290]  [6.3150]
    subplot(1,4,frequency), ft_topoplotER(plotcfg,fd_all)
end


% bar(Left_eye(1),0.01,0.1, 'r')
% bar(Left_eye(2),0.1,0.1, 'r')
% bar(Right_eye(1),0.1,0.1, 'b')
% bar(Right_eye(2),0.1,0.1, 'b')
% bar(Right_eye(1)*4,0.1,0.1, 'b')
% bar(Mouth(1)+Right_eye(1),1,0.1, 'c')
% bar(Mouth(1)-Right_eye(1),1,0.1, 'c')
% bar(Mouth(1),1,0.1, 'g')
% bar(Mouth(2),1,0.1, 'g')
% % bar(Mouth(1)+Left_eye(1)+Right_eye(1),1,0.1, 'k')
% for ii=1:length(hz_intermod)  bar(hz_intermod(ii),1,0.1,'m'); end
% for ii=27:35
%     fd.label{ii}
%
% end
%
% % elec_18_T8
% figure,
% plot(fd.freq, squeeze(fd.cohspctrm(18,30,:)))
% title(sprintf('connectivity between %s and %s', fd.label{18}, fd.label{30}));
% xlabel('freq (Hz)')
% ylabel('coherence')

%% Feature Coherence : Is it ok with single trial analysis?..

% which electrode for reference?
% which connectivity index ?

whichConnect=2;
Connect_index_cell={'plvspctrm', 'cohspctrm','powcorrspctrm','crsspctrm'}
plotCFG_index={'plvspctrm','cohspctrm','powspctrm','crsspctrm'}

cfg_index= {'plv','coh','powcorr','powcorr_ortho','csd'};


clear plotcfg fd Feature_topo
% rng('default')
featureLabel={'Left eye', 'Right eye', 'Mouth','Right versus Left Eye','Right eye versus Mouth'};
plotcfg.channel          = [1:64]%freqcfg.channel;
plotcfg.refchannel       = {'O1','O2','Oz'}; % T8 PO8
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='South';
plotcfg.layout=lay;
plotcfg.parameter= Connect_index_cell{whichConnect};%'cohspctrm';
plotcfg.showlabels ='yes';

clear Feature_topo

nChannels=size(freqPowNcsd.fourierspctrm,2);
nFreqs=size(freqPowNcsd.freq,2);
for ii=1:3
    %     Feature_topo{ii}=fd;
    %     Feature_topo{ii}.cohspctrm=zeros(nChannels,nChannels,nFreqs);
    Feature_topo{ii}.spctrm=zeros(nChannels,nChannels);
end

fd_all.coh=zeros(nChannels,nChannels,nFreqs);

blockCounter=1;
trialBlckCounter=0;
harmon=0;
for trial=1:16%size(freq_F{1}.powspctrm,1)
    
    cfg            = [];
    cfg.method     = cfg_index{whichConnect}%'coh';%'wpli_debiased';%'powcorr_ortho';%'plv'%'coh';'powcrr'
    cfg.complex    = 'imag'
    cfg.removemean  = 'yes';
    cfg.keeptrials='yes';
    cfg.trials     = [trial];
    cfg.channel    =[1:64]%freqcfg.channel%posterior_elec;
    fd             = ft_connectivityanalysis(cfg, freqPowNcsd );
    
    
    fd_all.coh=fd_all.coh+eval(sprintf('fd.%s',Connect_index_cell{whichConnect}));
    trialBlckCounter=trialBlckCounter+1;
    eval(sprintf('load(''MARIEEVE_SSVEPWINK%d_SSVEP_0_1.mat'') ', blockCounter))
    
    
    if  trial== 5 || trial== 9 || trial== 13
        eval(sprintf('load(''MARIEEVE_SSVEPWINK%d_SSVEP_0_1.mat'') ', blockCounter))
        %          load MARIEEVE_SSVEPINDUCTION3_SSVEP_0_0_0.mat
        
        blockCounter=blockCounter+1;
        if blockCounter == 2
            
        else
            
        end
        trialBlckCounter=1;
    end
    
    
    trial_features=cell2mat(data{trialBlckCounter}.freq);
    harmon_feat=cell2mat(data{trialBlckCounter}.freq)*2;
    taggs=[trial_features harmon_feat];
    for frequency=1:6
        
        C = abs(bsxfun(@minus,taggs(frequency)',fd.freq));
        [~,index_freq] = min(C(:,1:size(C,2)));
        if frequency ==3 || frequency ==6
            Feature_topo{3}.spctrm(:,:)=Feature_topo{3}.spctrm(:,:)+eval(sprintf('(fd.%s(:,:,index_freq))'' ',Connect_index_cell{whichConnect}));
        else
            Feature_topo{mod(frequency,3)}.spctrm(:,:)=Feature_topo{mod(frequency,3)}.spctrm(:,:)+eval(sprintf('(fd.%s(:,:,index_freq))'' ',Connect_index_cell{whichConnect}));
        end
    end
end
% % % % % % %
for feature=1:3, Feature_topo{feature}.spctrm(:,:)=Feature_topo{feature}.spctrm(:,:)/(16); end
fd_all.coh=fd_all.coh/16;
%%
plotCFG_index={'plvspctrm','cohspctrm','powcorrspctrm'}

clear Feature_topoCoh
plotcfg.parameter=plotCFG_index{whichConnect};% Connect_index_cell{whichConnect};%'cohspctrm';
plotcfg.refchannel       = {'Oz'};%{'O1'};%{'FT7','F7','T7'}%{'FT8','F8','T8'}%{'O1','O2','Oz'};% {'Oz',}%% T8 PO8
figure,
for feature=1:3
    
    
    Feature_topoCoh{feature}=fd;
    eval(sprintf('Feature_topoCoh{feature}.%s(:,:,:)=repmat(Feature_topo{feature}.spctrm(:,:),1,1,nFreqs);',Connect_index_cell{whichConnect}));
    
    
    plotcfg.xlim           = [2.3 2.4];%[6.3 6.32];8]  [3.5290]  [6.3150]
    subplot(1,3,feature), ft_topoplotER(plotcfg,fd)
    ft_topoplotER(plotcfg, Feature_topoCoh{feature})%, title(sprintf('%s',featureLabel{feature}));
    
end
%%
plotcfg.refchannel       ={'Oz','O1','O2'}%{'P2','P1','Pz'}%{'FT8','F8','T8'}% {'Oz','O1','O2'};%{'O
plotcfg.colorbar='East';
figure,
Feature_topoMean=fd;
Feature_topoMean.spctrm=zeros(size(Feature_topo{feature}.spctrm,1),size(Feature_topo{feature}.spctrm,1))
for feature= 1:3, Feature_topoMean.spctrm(:,:)=Feature_topoMean.spctrm(:,:)+Feature_topo{feature}.spctrm(:,:); end
Feature_topoMean.spctrm(:,:)=Feature_topoMean.spctrm(:,:)/3;
eval(sprintf('Feature_topoMean.%s(:,:,:)=repmat(Feature_topoMean.spctrm(:,:),1,1,nFreqs);',Connect_index_cell{whichConnect}));
ft_topoplotER(plotcfg, Feature_topoMean)%, title(sprintf('%s',featureLabel{feature}));

% eval(sprintf('Feature_topoMean.%s(:,:,:)=fd_all.coh;',Connect_index_cell{whichConnect}));







%% INTERMOD CONNECTIVITY


whichConnect=2;
Connect_index_cell={'plvspctrm', 'cohspctrm','powcorrspctrm','crsspctrm'}
plotCFG_index={'plvspctrm','cohspctrm','powspctrm','crsspctrm'}

cfg_index= {'plv','coh','powcorr','powcorr_ortho','csd'};


clear plotcfg fd Feature_topo
% rng('default')
featureLabel={'Left eye', 'Right eye', 'Mouth','Right versus Left Eye','Right eye versus Mouth'};
plotcfg.channel          = [1:64]%freqcfg.channel;
plotcfg.refchannel       = {'O1','O2','Oz'}; % T8 PO8
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='South';
plotcfg.layout=lay;
plotcfg.parameter= Connect_index_cell{whichConnect};%'cohspctrm';
plotcfg.showlabels ='yes';


clear Feature_topo

nChannels=size(freqPowNcsd.fourierspctrm,2);
nFreqs=size(freqPowNcsd.freq,2);
for ii=1:3
    %     Feature_topo{ii}=fd;
    %     Feature_topo{ii}.cohspctrm=zeros(nChannels,nChannels,nFreqs);
    Feature_topo{ii}.cohspctrm=zeros(nChannels,nChannels);
    
end

fd_all.coh=zeros(nChannels,nChannels,nFreqs);

blockCounter=1;
trialBlckCounter=0;
harmon=0;
for trial=1:8%size(freq_F{1}.powspctrm,1)
    
    cfg            = [];
    cfg.method     = 'coh';%'wpli_debiased';%'powcorr_ortho';%'plv'%'coh';'powcrr'
    cfg.complex    = 'imag'
    cfg.removemean  = 'yes';
    cfg.keeptrials='yes';
    cfg.trials     = [trial];
    cfg.channel    =[1:64]%freqcfg.channel%posterior_elec;
    fd             = ft_connectivityanalysis(cfg, freqPowNcsd );
    
    
    fd_all.coh=fd_all.coh+fd.cohspctrm;
    trialBlckCounter=trialBlckCounter+1;
    eval(sprintf('load(''MARIEEVE_SSVEPWINK%d_SSVEP_0_1.mat'') ', blockCounter))
    
    
    if  trial== 5 || trial== 9 || trial== 13
        eval(sprintf('load(''MARIEEVE_SSVEPWINK%d_SSVEP_0_1.mat'') ', blockCounter))
        %          load MARIEEVE_SSVEPINDUCTION3_SSVEP_0_0_0.mat
        
        blockCounter=blockCounter+1;
        if blockCounter == 2
            
        else
            
        end
        trialBlckCounter=1;
    end
    
    
    trial_features=cell2mat(data{trialBlckCounter}.freq);
    harmon_feat=cell2mat(data{trialBlckCounter}.freq)*2;
    taggs=[trial_features harmon_feat];
    intermodulaires=intermod([trial_features harmon_feat],[1.5 20]);
    for frequency=1:6
        
        C = abs(bsxfun(@minus,intermodulaires(frequency)',fd.freq));
        [~,index_freq] = min(C(:,1:size(C,2)));
        if frequency ==3 || frequency ==6
            Feature_topo{3}.cohspctrm(:,:)=Feature_topo{3}.cohspctrm(:,:)+(fd.cohspctrm(:,:,index_freq))';
        else
            Feature_topo{mod(frequency,3)}.cohspctrm(:,:)=Feature_topo{mod(frequency,3)}.cohspctrm(:,:)+(fd.cohspctrm(:,:,index_freq))';
        end
    end
end
% % % % % % %
for feature=1:3, Feature_topo{feature}.cohspctrm(:,:)=Feature_topo{feature}.cohspctrm(:,:)/sqrt(16); end
%%

plotcfg.refchannel       = {'Oz'};%{'FT7','F7','T7'}%{'FT8','F8','T8'}%{'O1','O2','Oz'};% {'Oz',}%% T8 PO8
figure,
for feature=1:3
    
    
    Feature_topoCoh{feature}=fd;
    Feature_topoCoh{feature}.cohspctrm(:,:,:)=repmat(Feature_topo{feature}.cohspctrm(:,:),1,1,nFreqs);
    
    
    plotcfg.xlim           = [2.3 2.4];%[6.3 6.32];8]  [3.5290]  [6.3150]
    subplot(1,3,feature), ft_topoplotER(plotcfg,fd)
    ft_topoplotER(plotcfg, Feature_topoCoh{feature})%, title(sprintf('%s',featureLabel{feature}));
    
end

%%


trial_features=cell2mat(data{trialBlckCounter}.freq);
harmon_feat=cell2mat(data{trialBlckCounter}.freq)*2;
taggs=[trial_features harmon_feat];
intermodulaires=intermod([trial_features harmon_feat],[1.5 20])

%%
plotcfg.channel          = [1:64]%freqcfg.channel;
plotcfg.refchannel       = {'Oz','O1','O2'}%{'F8','T8'}%{'O1','O2','Oz'};%{'O2'}% % T8 PO8
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='South';
plotcfg.layout=lay;
plotcfg.parameter= 'cohspctrm'%'plvspctrm'%'cohspctrm';
plotcfg.showlabels ='yes';



C = abs(bsxfun(@minus,intermodulaires(frequency)',fd.freq));
[~,index_freq] = min(C(:,1:size(C,2)));

meanSpctr=eval(sprintf('(fd.%s(:,:,index_freq))'' ',Connect_index_cell{whichConnect}));

freq=mat2cell(intermodulaires,1);
figure,
for frequency=1:21
    
    C = abs(bsxfun(@minus,intermodulaires(frequency)',fd.freq));
    [~,index_freq] = min(C(:,1:size(C,2)));
    
    meanSpctr.spctrm=meanSpctr.spctrm+eval(sprintf('(fd.%s(:,:,index_freq))'' ',Connect_index_cell{whichConnect}));
    
    %
    % plotcfg.xlim           = [freq{1}(frequency)-0.01 freq{1}(frequency)+0.01];%[6.3 6.32];8]  [3.5290]  [6.3150]
    %  subplot(5,5,frequency), ft_topoplotER(plotcfg,fd)
end

% bar(Left_eye(1),0.01,0.1, 'r')
% bar(Left_eye(2),0.1




%% TIME FREQUENCY ANALYSIS

posterior_elec = {'T7','T8','P7','P3','Pz', 'P4', 'P8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'};
%       posterior_elec = {'O1', 'Oz', 'O2'};
%     posterior_elec = {'Pz','Cz','CPz'};
%     posterior_elec ={'PO8','PO6','TP8','T8'};
% TrialTime=raw_data.time{1}(end);
fs=1;% 1/TrialTime;%2.5;


cfg              = [];
cfg.keeptrials   ='no';
cfg.output       = 'pow';
cfg.baseline     = [];%[.52 .6];
cfg.baselinetype = [];%'relchange'; %'db';% relchange chck out laurent used it
cfg.channel      = posterior_elec;%[1:64]%freq.channel%{'Oz','Pz','T8','PO8'};
cfg.method       = 'mtmconvol';
cfg.pad          = 'nextpow2';
cfg.taper        = 'hanning';
cfg.foi          = 1:fs:20;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*.5;   % length of time window = 0.5 sec 7./cfg.foi;
cfg.toi          = 0.5:0.05:1.25;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.trials='all';%find(TrialTagg{subject}==3);%'all'

clear TFR
for perm=1:6
    SubCounter=0;  
for subject=2:15
    SubCounter=SubCounter+1;
    cfg.trials=find(InterResponsesTrial{subject}==perm)
    eval(sprintf('load(''/Volumes/Simon Backup 2TB/DATA_EEG/TrialedResponse_avg_data_SSVEP_%s_%s'',''raw_data'')',Subjects{subject}{1},Subjects{subject}{3}))
%     freqcfg.trials='all';
    if ~isempty(cfg.trials)
    TFR{perm}{SubCounter} = ft_freqanalysis(cfg, raw_data);
    else
    TFR{perm}{SubCounter}=[];
    end
end
end
% for trial=trials
%
%
% cfg.trials=[trial];
% TFRhann{trial} = ft_freqanalysis(cfg, raw_data);
% end
%%
%     figure,
% for trial=1%:5
% ft_singleplotTFR(cfg,TFRhann{trial});
% end
% figure,
% % ft_singleplotTFR(cfg,TFRhann);
% clear DIFF
% subCounter=0;
% for subject=Subj_cond1(1:5)
%     subCounter=subCounter+1;
%     
%     DIFF{subCounter}=TFRhann{subCounter}{perm};
%     DIFF{subCounter}.powspctrm=(zscore(TFRhann{subCounter}{perm}.powspctrm(:,:,:))-zscore(TFRhann{subCounter+30}{perm}.powspctrm(:,:,:)))/2;
%     DIFF_all(subCounter,:,:,:)=(zscore(TFRhann{subCounter}{perm}.powspctrm(:,:,:))-zscore(TFRhann{subCounter+30}{perm}.powspctrm(:,:,:)))/2;
%     
% end
for perm=1:6
switch perm
    case 5
        [grandavg{perm}] = ft_freqgrandaverage(cfg, TFR{perm}{1},  TFR{perm}{2},  TFR{perm}{3},TFR{perm}{4} ,TFR{perm}{5},TFR{perm}{6},TFR{perm}{7},TFR{perm}{8},TFR{perm}{10},TFR{perm}{11},TFR{perm}{12},TFR{perm}{13},TFR{perm}{14})
    case 1
        [grandavg{perm}] = ft_freqgrandaverage(cfg, TFR{perm}{1},  TFR{perm}{3},  TFR{perm}{4},TFR{perm}{5} ,TFR{perm}{6},TFR{perm}{7},TFR{perm}{8},TFR{perm}{10},TFR{perm}{11},TFR{perm}{12},TFR{perm}{13},TFR{perm}{14})
    otherwise
    [grandavg{perm}] = ft_freqgrandaverage(cfg, TFR{perm}{1},TFR{perm}{2},  TFR{perm}{3},  TFR{perm}{4},TFR{perm}{5} ,TFR{perm}{6},TFR{perm}{7},TFR{perm}{8},TFR{perm}{10},TFR{perm}{11},TFR{perm}{12},TFR{perm}{13},TFR{perm}{14})
end

% [grandavg] = ft_freqgrandaverage(cfg, DIFF{1}, DIFF{2},DIFF{3}, DIFF{4}, DIFF{5});%DIFF{6},DIFF{7},DIFF{8},DIFF{9},DIFF{10},DIFF{11},DIFF{12},DIFF{13},DIFF{14},DIFF{15},DIFF{16},DIFF{17},DIFF{18});
% [grandavg{perm}] = ft_freqgrandaverage(cfg, TFR{1}{perm},  TFR{2}{perm},  TFR{3}{perm},TFR{4}{perm} ,TFR{5}{perm},TFR{6}{perm},TFR{7}{perm},TFR{8}{perm},TFR{9}{perm},TFR{10}{perm},TFR{11}{perm},TFR{12}{perm},TFR{13}{perm},TFR{14}{perm})
% [grandavg] = ft_freqgrandaverage(cfg, TFRhann{1+30}{perm},  TFRhann{2+30}{perm},  TFRhann{3+30}{perm},TFRhann{4+30}{perm} ,TFRhann{5+30}{perm});

% [grandavg] = ft_freqgrandaverage(cfg, TFRhann{1}{perm},TFRhann{3}{perm},TFRhann{5}{perm})
% [grandavg] = ft_freqgrandaverage(cfg, TFRhann{12}, TFRhann{14}, TFRhann{16}, TFRhann{18}, TFRhann{20}, TFRhann{22}, TFRhann{24}, TFRhann{26},TFRhann{28})

% DIFF.powspctrm=(TFRhann{subject}.powspctrm(:,:,:)-TFRhann{subject+1}.powspctrm(:,:,:))/2;

% cfg.layout=lay;%electrodes;
%  [grandavg1] = ft_freqgrandaverage(cfg, DIFF)
plotcfg=[];
plotcfg.channel={'O1', 'Oz', 'O2'};%{'PO8','PO6','P6','O1', 'Oz', 'O2'};%%{'PO7','PO5','P5', 'Oz'};%%{'PO8','PO6','P6','O1', 'Oz', 'O2'};%{'O1', 'Oz', 'O2'}; {'PO8','PO6','TP8','T8'};

% figure, 
% ft_singleplotTFR(plotcfg, grandavg{perm}),line([0 1], [4 4]),line([0 1], [6 6]),line([0 1], [8 8])
% figure, imagesc(squeeze(nanmean(grandavg{perm}.powspctrm(:,:,:)))),line([0 21], [4 4]),line([0 21], [6 6]),,line([0 21], [8 8])
figure, imagesc(squeeze(nanmean(grandavg{perm}.powspctrm(:,:,:)))-squeeze(nanmean(grandavg{1}.powspctrm(:,:,:)))),line([0 21], [4 4]),line([0 21], [6 6]),line([0 21], [8 8])
set(gca,'XTick',1:length(grandavg{perm}.time),'XTickLabels',grandavg{perm}.time)
set(gca,'YTick',1:length(grandavg{perm}.freq),'YTickLabels',grandavg{perm}.freq')
end
% % % % save('DIFF_TimeFreqOcc_WINKvsBLINK.mat','TFRhann','DIFF','cfg')
%%
cfg = [];
load('quickcap64.mat') % load 2D layout of electrodes
cfg.elec=electrodes;
cfg.lay=[];%lay;
cfg.baseline     = [-0.5 -0.25];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.1 1];
cfg.zlim         = [];
cfg.ylim         = [4.22 4.5];
cfg.ylim         = [11.22 11.52];
cfg.marker       = 'on';
figure
ft_topoplotTFR(cfg, grandavg);

%  [grandavg2] = ft_freqgrandaverage(cfg, TFRhann{subject+1})
% figure,
% ft_singleplotTFR(cfg,TFRhann{subject+1});
%%
C = abs(bsxfun(@minus,6.529,TFR.freq));
[~,index_freqToFeatures] = min(C(:,1:size(C,2)));
cfg = [];
cfg.channel      = {'O1','PO1'}%{'T7','T8','P7','P3','Pz', 'P4', 'P8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'}%{'Oz','POz','O1','O2'}%{'Oz'};%{'Oz','POz','O1','O2'}%{'Pz','Oz'}%{'PO8'};
cfg.baseline     = [-0.5 -0.1];
% cfg.baseline      = [-inf 0];
cfg.baselinetype = 'relchange';%'absolute'; %'db';%
cfg.zlim         = 'maxmin'
cfg.showlabels   = 'yes';
cfg.marker       = 'on';
cfg.layout       = lay;
figure
%ft_multiplotTFR(cfg, TFRhann);
ft_singleplotTFR(cfg,TFR);
%ft_topoplotTFR(cfg,TFRhann);


figure, plot(mean(squeeze(TFR.powspctrm(:,170,:))))
hold on
plot(mean(squeeze(TFR.powspctrm(:,370,:))))
hold on
plot(mean(squeeze(TFR.powspctrm(:,467,:))))
% plot(squeeze(mean(mean(squeeze(TFRhann.powspctrm(:,601:734,:)))))-mean(squeeze(TFRhann.powspctrm(:,170,:)))')


%%

% Highlighting artifacts... 1) EOG ---> 2) Muscle---> 3) Jump

% EOG ARTIFACT HIGHLIGHTING
cfg.preproc.bpfilter    = 'yes';
cfg.preproc.bpfilttype  = 'but';
cfg.preproc.bpfreq      = [1 15];
cfg.preproc.bpfiltord   = 5.3;
cfg.preproc.rectify     = 'yes';



cfg.artfctdef.eog.cutoff       = 4;
cfg.artfctdef.eog.trlpadding   = 0.5;
cfg.artfctdef.eog.fltpadding   = 0.1;
cfg.artfctdef.eog.artpadding   = 0.1;

cfg.artfctdef.eog.channel      = cfg.hdr.label;



[cfg, artifact_EOG] = ft_artifact_eog(cfg);
cfg.artfctdef.eog.artifact=artifact_EOG;

%
% % MUSCLE ARTIFACTS Parameter

cfg.preproc.bpfilter    = 'yes'
cfg.preproc.bpfreq      = [110 140]
cfg.preproc.bpfiltord   =  8
cfg.preproc.bpfilttype  = 'but'
cfg.preproc.rectify     = 'yes'

cfg.artfctdef.muscle.channel     = cfg.hdr.label;
%cfg.artfctdef.muscle.cutoff      = z-value at which to threshold (default = 4)
cfg.artfctdef.muscle.trlpadding  = 0.1
cfg.artfctdef.muscle.fltpadding  = 0.1
cfg.artfctdef.muscle.artpadding  = 0.1

[cfg, artifact_muscle] = ft_artifact_muscle(cfg);
cfg.artfctdef.muscle.artifact =artifact_muscle;
%
% JUMP ARTIFACT

cfg.artfctdef.jump.channel       = cfg.hdr.label;
cfg.artfctdef.jump.medianfilter  = 'yes';
cfg.artfctdef.jump.medianfiltord = 9;
cfg.artfctdef.jump.absdiff       = 'yes';
[cfg, artifact_jump] = ft_artifact_jump(cfg);
cfg.artfctdef.minaccepttim = 5;
cfg.artfctdef.feedback ='yes';
cfg.artfctdef.reject='partial';%'partial';
Cleanercfg = ft_rejectartifact(cfg);


Cleaner_data=ft_preprocessing(Cleanercfg,raw_data);

%%

% load 'Features_frequenciesSubjects_raw.mat'
% load Features_frequenciesSubjects_SmoothRejection.mat





% plotcfg.zlim = [0 0.5];


CONDITION=1;

% switch CONDITION
%     case 1
%         Subj_cond1=[8 12 14 16 18]
% load Features_frequenciesSubjects.mat
%     case 2
% %         load 'Features_frequenciesSubjectsCOND2.mat'
%         Subj_cond2=[9 13 15 17 19]
% end

%   Subj_cond1=[8 12 14 16]
%   Subj_cond2=[13 15 17]

Subj_cond2=[8]
%  Subj_cond2=[13 15 17]

%%
CONDITION=1;
rng('default')
clear plotcfg
featureLabel={'Left eye', 'Right eye', 'Mouth','Right versus Left Eye','Right eye versus Mouth','Both eyes vs Mouth'};
plotcfg.channel= freqcfg.channel;
plotcfg.colorbar='East';
plotcfg.layout=lay;
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';
plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='East';
load('quickcap64.mat') % load 2D layout of electrodes
plotcfg.layout=lay;
plotcfg.showlabels ='yes';

Subj_cond1(Subj_cond1==35)=[];

for ii=1:3,
    AllSub_Features{ii}=freq_F{8}{1};
    AllSub_Features{ii}.powspctrm=zeros(1,size(freq_F{8}{1}.powspctrm,2));
end

for subject=Subj_cond1(1:end)
    clear Feature_topo
    for ii=1:3
        
        Feature_topo{ii}=freq_F{subject}{1};
        Feature_topo{ii}.powspctrm=zeros(1,size(freq_F{subject}{1}.powspctrm,2));
        
    end
    
    blockCounter=1;
    trialBlckCounter=0;
    harmon=0;
    % load CATHERINESSVEPBLINKCOND2Right_blck1_SSVEP_0_1_2.mat
    for trial=1:size(freq_F{subject}{1}.powspctrm,1)
        trialBlckCounter=trialBlckCounter+1;
        
        if CONDITION==1
            eval(sprintf('load(''%s_blck%d_SSVEP_0_1_1.mat'') ',Subjects{subject}{3} ,blockCounter))
        else
            eval(sprintf('load(''%s_blck%d_SSVEP_0_1_2.mat'') ',Subjects{subject}{3} ,blockCounter))
        end
        
        if  trial== 5 || trial== 9 || trial== 13
            blockCounter=blockCounter+1;
            
            if CONDITION==1
                eval(sprintf('load(''%s_blck%d_SSVEP_0_1_1.mat'') ',Subjects{subject}{3} ,blockCounter))
            else
                eval(sprintf('load(''%s_blck%d_SSVEP_0_1_2.mat'') ',Subjects{subject}{3} ,blockCounter))
            end
            
            
            trialBlckCounter=1;
        end
        
        
        trial_features=cell2mat(data{trialBlckCounter}.freq);
        harmon_feat=cell2mat(data{trialBlckCounter}.freq)*2;
        for frequency=1:6
            C = abs(bsxfun(@minus,trial_features',freq_F{subject}{frequency}.freq));
            [~,index_freqToFeatures] = min(C(:,1:size(C,2)));
            if ~any(index_freqToFeatures)
                C = abs(bsxfun(@minus,harmon_feat',freq_F{subject}{frequency}.freq));
                [~,index_freqToFeatures] = min(C(:,1:size(C,2)))
                harmon=harmon+1; %  this/3 should be equal to number of trials
            end
            Feature_topo{index_freqToFeatures}.powspctrm(:)=Feature_topo{index_freqToFeatures}.powspctrm(:)+freq_F{subject}{frequency}.powspctrm(trial,:)';
            
        end
    end
    %%%
    
    for feature=1:3
        
        Feature_topo{feature}.powspctrm(:)=Feature_topo{feature}.powspctrm(:)/((size(freq_F{subject}{feature}.powspctrm,1)));
        AllSub_Features{feature}.powspctrm(:)=AllSub_Features{feature}.powspctrm(:)+Feature_topo{feature}.powspctrm(:);
        
        
        %         figure,ft_topoplotER(plotcfg, Feature_topo{feature}), title(sprintf('%s',featureLabel{feature}));
        
    end
    
    %
    Feature_Diff{1}=Feature_topo{2};
    Feature_Diff{2}=Feature_Diff{1};
    Feature_Diff{3}=Feature_Diff{1};
    %
    
    Feature_Diff{1}.powspctrm(:)=Feature_topo{2}.powspctrm(:)-Feature_topo{1}.powspctrm(:);
    Feature_Diff{2}.powspctrm(:)=Feature_topo{2}.powspctrm(:)-Feature_topo{3}.powspctrm(:);
    Feature_Diff{3}.powspctrm(:)=((Feature_topo{2}.powspctrm(:)+Feature_topo{1}.powspctrm(:))/2)-Feature_topo{3}.powspctrm(:);
    
    LEFT_INDUCTION=Feature_Diff{1};
    
    % figure,ft_topoplotER(plotcfg, Feature_Diff{1}), title(sprintf('%s',featureLabel{4}));
    % figure,ft_topoplotER(plotcfg, Feature_Diff{2}), title(sprintf('%s',featureLabel{5}));
    % figure,ft_topoplotER(plotcfg, Feature_Diff{3}), title(sprintf('%s',featureLabel{6}));
    %
    
    
    
end

for feature=1:3
    AllSub_Features{feature}.powspctrm(:)=AllSub_Features{feature}.powspctrm(:)/length(Subj_cond1);
    figure,ft_topoplotER(plotcfg, AllSub_Features{feature}), title(sprintf('%s',featureLabel{feature}));
end
% AllSub_Features{feature}=Inter_difference{subject}{1};
% All_cfgDiff.powspctrm=MDiff.powspctrm;
% save('topofreq.mat','AllSub_Features')
%  figure,ft_topoplotER(plotcfg, All_cfgDiff), title(sprintf('%s',interLabel{6}))
% % save('inductionMarie.mat','RIGHT_INDUCTION','LEFT_INDUCTION')



%% INTERMODULATORY FREQUENCIES :
clear freqcfg

plotcfg.interpolatenan ='yes';
plotcfg.interplimits = 'head';
plotcfg.interpolation ='v4';%'nearest'
plotcfg.colorbar='East';
plotcfg.layout=lay;
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';

hz_intermod=intermod([3.529 8 6.315],[1.5 20]);

TrialTime=50;%size(cell2mat(raw_data.time(1)),2)/1024;

clear freq freqcfg

posterior_elec = {'T7','T8','P7','P3','Pz', 'P4', 'P8', 'POz', 'O1', 'Oz', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'TP7', 'TP8','PO7','PO8'};
% [1:64]
freqcfg.channel=ft_channelselection({'all','-EOG','-M1','-M2','-AF7','-AF5','-Fpz','-AF3','-Fp2','-AF4','-Fp1','-F3','-F5','-F7','-F1','-Fp1','-EOG'},cfg.hdr.label(all_electrodes))
load('quickcap64.mat') % load 2D layout of electrodes
plotcfg.layout=lay;
plotcfg.colorbar='East';
plotcfg.parameter= 'powspctrm';
plotcfg.showlabels ='yes';
%     load('quickcap64.mat') % load 2D layout of electrodes



freqcfg.layout=lay;

cfg.hdr.label(1:64)
Subj_cond1=[2 3 4 8] % Cond 1 only
Subj_cond2=[6 9] % Cond 2 only

clear freqcfg freq
for subject=Subj_cond2
    %     eval(sprintf('load(''Trialed_PostIcadata%s_%s.mat'')',Subjects{subject}{1},Subjects{subject}{5}))
    eval(sprintf('load(''Trialed_PostIcadataSMOOTHRej%s_%s'')',Subjects{subject}{1},Subjects{subject}{5}))
    
    for INTER=1:size(hz_intermod,2)
        
        %     freqcfg.foi   = hz_intermod(INTER);
        freqcfg.foilim  = [hz_intermod(INTER)-0.015 hz_intermod(INTER)+0.015];
        
        freqcfg.taper ='dpss';
        freqcfg.tapsmofrq  = 1/TrialTime;
        freqcfg.pad='nextpow2'
        
        freqcfg.method= 'mtmfft';%'wavelet';
        
        freqcfg.keeptrials='yes';
        [freq{subject}{INTER}] = ft_freqanalysis(freqcfg, data_postica)
        
        %     figure,ft_topoplotER(plotcfg, freq{INTER}),% title(sprintf('%s',interLabel{INTER}))
        
    end
    
end

% save('Intermodulary_frequenciesSubjectsCOND2.mat','freq','freqcfg')
%% SOURCE ANALYSIS TEST

% read in the template MRI
addpath /Users/SimonFaghel-Soubeyrand/Documents/MATLAB/fieldtrip-master
       mri =ft_read_mri('/Users/SimonFaghel-Soubeyrand/Documents/MATLAB/fieldtrip-master/external/spm8/templates/T1.nii');
load('/Volumes/Simon Backup 2TB/DATA_EEG/Anatomical Sources/segmentedmri.mat');
% segment the MRI
%         cfg = [];cfg.output = {'brain' 'skull' 'scalp'};   segmentedmri = ft_volumesegment(cfg, mri);
%         cfg = [];cfg.output = {'brain'};   segmentedmri = ft_volumesegment(cfg, mri);
% 
% cfg = [];
% cfg.write      = 'no';
% [segmentedmri] = ft_volumesegment(cfg, mri);
% create the headmodel (BEM)
        cfg = [];
        cfg.method ='openmeeg';%'openmeeg'; % TODO FIXME download openmeeg
%             cfg.method ='dipoli'; % dipoli only works under linux

        hdm = ft_prepare_headmodel(cfg, segmentedmri);
%         hdm = ft_prepare_headmodel(cfg, mri);
%         >>
% >>         elec = ft_read_sens('standard_1010.elc');
% >>         hdm = ft_convert_units(hdm, elec.unit);
% >>
% >>         cfg = [];
% >>         cfg.grid.xgrid = -125:8:125;
% >>         cfg.grid.ygrid = -125:8:125;
% >>         cfg.grid.zgrid = -125:8:125;
% >>         cfg.grid.tight = 'yes';
% >>         cfg.grid.unit = hdm.unit;
% >>         cfg.inwardshift = -1.5;
% >>         cfg.vol = hdm;
% >>         grid = ft_prepare_sourcemodel(cfg)
% >>         grid = ft_convert_units(grid, elec.unit);
% >>
% >>         figure;
% >>         hold on;
% >>         ft_plot_mesh(hdm.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha',
% >>         0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
% >>         ft_plot_mesh(grid.pos(grid.inside, :));
% >>
% >>         % this step is not necessary, cause you will see that
% >>         everything is
% >>         already aligned
% >>         cfg = [];
% >>         cfg.method = 'interactive';
% >>         cfg.elec = elec;
% >>         cfg.headshape = hdm.bnd(1);
% >>         tmp = ft_electroderealign(cfg);
% >>         elec = tmp; % I had a bug here that I couldn't assign elec
% >>         directly
% >>
% >>         %% verify location of occipital electrodes
% >>
% >>         occ_elec = elec;
% >>         occ_chan = ft_channelselection({'O*', 'PO*', 'Cz*', 'Fz*'},
% >>         elec.label);
% >>         occ_idx = match_str(elec.label, occ_chan);
% >>         occ_elec.chanpos = occ_elec.chanpos(occ_idx, :);
% >>         occ_elec.elecpos = occ_elec.elecpos(occ_idx, :);
% >>         occ_elec.label = occ_elec.label(occ_idx, :);
% >>         figure;
% >>         ft_plot_sens(occ_elec)
% >>         hold on;
% >>         ft_plot_vol(ft_convert_units(hdm, elec.unit))/
% 
% 
% 
% 
% 
% 
