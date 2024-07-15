% Simulation program for WTC-GLM and CTW-GLM
close all

%% ADDPATH
% Please add path to ...
%  Homer2 (hmrIntensity2OD.m, hmrOD2Conc.m, addSimHRF.m) from https://homer-fnirs.org/
%  NIRS-SPM (detrend_wavelet_MDL.m) from https://www.nitrc.org/projects/nirs_spm/
%  huppertt-nirs-toolbox (ar_fit.m) from https://github.com/huppertt/nirs-toolbox
%   to fix the AR order, modify the function

% Set folder name of resting data
root_dir_restdata = '\Resting\'; % change

% Simulation settings
L = 300;% 5 min

N_sh = 100;
N_sim = 500; %

Fs = 10; % 10 Hz
delay_true = 0; % no-delay
prob_gens = [(0:0.0025:0.01)']; % GLM

norm_prob = @(x) x/sum(x(:));
pick1ele = @(A) A(randi(numel(A)));
pick1 = @(A) A(1);
extm = @(x)(x(:));

Pmed=100;

% Load HRFs (50 Hz)
hrf_100 = load([root_dir_restdata,'code\hrf_100.mat'],'hrf');

dir_NIRS_home1 = [root_dir_restdata,'resting_state_1\'];
dir_NIRS_each1 = dir([dir_NIRS_home1,'Subj*']);
dir_NIRS_each1 = dir_NIRS_each1(1:end-1); % exclude excel

dir_NIRS_home2 = [root_dir_restdata,'resting_state_2\'];
dir_NIRS_each2 = dir([dir_NIRS_home2,'Subj*']);
dir_NIRS_each2 = dir_NIRS_each2(1:end-1); % exclude excel

dir_NIRS_each = cat(1,dir_NIRS_each1,dir_NIRS_each2);

nirs_org_50Hz = cell(size(dir_NIRS_each ,1),1);

hrf_100.hrf.hrf_conc = resample(hrf_100.hrf.hrf_conc ,1,5); % downsampling 50 Hz -> 10 Hz
ratio_hrf_Oxy = max(hrf_100.hrf.hrf_conc(:,1))/max(hrf_2sec);
ratio_hrf_Deoxy = -max(-hrf_100.hrf.hrf_conc(:,2))/max(hrf_2sec);

hrf_100.hrf.hrf_d = resample(hrf_100.hrf.hrf_d ,1,5);

hrf_s = hrf_100.hrf.hrf_conc(:,1);
L_hrf = size(hrf_s,1); %176

%% Load data for baseline
dc_org_10Hz_30mm = cell(size(dir_NIRS_each ,1),1);
for k_each = 1:size(dir_NIRS_each ,1)
    dir_NIRS = fullfile(dir_NIRS_each(k_each).folder,dir_NIRS_each(k_each).name);
    filename = [dir_NIRS,'\resting'];
    
    % Load nirs (Resting)
    if ~exist([filename,'_bwards.nirs'],'file')
        snirf2nirs(filename);
    end
    nirs = load([filename,'_bwards.nirs'], '-mat');
    
    nirs.s = zeros(size(nirs.t));
    nirs.s(5000,1)=1; % event(dummy)
    
    % settings:
    flag_prune = false; % channel pruning turned on (add synthetic HRF only on half of good channels)
    utest = false; %unit test, saves hrf onto constant signal (ignores resting data);
    
    % add HRF(dummy)
    nirs_hrf_100 = addSimHRF(nirs, hrf_100.hrf, utest, flag_prune);
    
    dod_org = hmrIntensity2OD(nirs.d );
    dc_org_50Hz = hmrOD2Conc(dod_org,nirs_hrf_100.SD, [6,6]);
    dc_org_10Hz_30mm{k_each,1} = resample(dc_org_50Hz(:,1:2,nirs_hrf_100.lstLongAct),1,5);
end


%% preprocessings
rs_base_detrended_n = cell(size(dir_NIRS_each ,1),1);
for k_each = 1:size(dir_NIRS_each ,1)
    % HDMS
    % Here, please use POTATo
    % data_dc_hdms <- (permute(dc_org_10Hz_30mm{k_each,1},[1,3,2])
    
    % Detrending
    % Remove trend by Wavelet-MDL
    bias_HR = detrend_waveletMDL(data_dc_hdms(:,:,1),ones(size(data_dc_hdms,1),1));
    rs_base_detrended = data_dc_hdms(:,:,1)-bias_HR;
    
    % exclude bias
    rs_base_detrended_n{k_each,1} = rs_base_detrended-repmat(mean(rs_base_detrended,1),[size(rs_base_detrended,1) 1]);
end


%% combinations
% + Exclude data less than 300 s (5 min) (- > total 27 sets)
% + Will use only combination between channels from different set
idx_sim = zeros(N_sim,4);
rowSizes = cellfun(@(x) size(x, 1), rs_base_detrended_n);
colSizes = cellfun(@(x) size(x, 2), rs_base_detrended_n);
idx_okset = find(rowSizes>=L*Fs);
comb_sets = nchoosek(idx_okset,2);

k_sim = 0; % initialize
while k_sim<(N_sim+100)
    k_sim = k_sim+1;
    
    if rem(k_sim,size(comb_sets,1))==1 % randmize
        disp(k_sim)
        comb_sets = comb_sets(randperm(size(comb_sets,1)),:);
    end
    tmp_comb = comb_sets(rem(k_sim-1,size(comb_sets,1))+1,:);
    
    tmp_used= idx_sim(((idx_sim(:,1)==tmp_comb(1))+(idx_sim(:,3)==tmp_comb(2))==2),[2,4]); % used
    if size(tmp_used,1)==0
        tmp_used = [0,0];
    end
    tmp_nonused1 = setdiff(1:colSizes(tmp_comb(1)),tmp_used(1));
    tmp_nonused2 = setdiff(1:colSizes(tmp_comb(2)),tmp_used(2));
    
    if size(tmp_nonused,2) > 0 && size(tmp_nonused2,2) > 0
        idx_sim(k_sim,:) = [tmp_comb(1),pick1ele(tmp_nonused1),tmp_comb(2),pick1ele(tmp_nonused2)];
    end
end


%% Make social signal
Observ = zeros(L*Fs,N_sim,size(prob_gens,1));
for k_ptr = 1:size(prob_gens,1)
    for k_sim = 1:N_sim
        % Social signal
        for k_L = 1:(L*Fs-L_hrf+1)
            if sum(sum(Observ(unique(max((k_L-1),1):k_L),k_sim,k_ptr)))==0 % no HRF
                tmp_O = rand(1)<=prob_gens(k_ptr);
                if tmp_O==1
                    Observ(k_L:(k_L+L_hrf-1),k_sim,k_ptr) = hrf_s;
                end
            end
        end
    end
end


%% Simulation
N_band = 75;
FOI = 12:68;

mat2vec = @(x)(x(:));
pickcorr2 = @(x)(x(2,1));

FOI_GLMeach = 36:84; % (12 oct.) 0.5 - 0.03 Hz
[~,~,FreqHz]=wcoherence(Observ(:,1,k_ptr),Observ(:,2,k_ptr),Fs, 'amor',1,12,'AR1');

N_sim2 = 100;

%% WTC
WTC4GLM.WTCpw = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim2,4);
WTC4GLM.WTCpw_synth = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim2,4);
WTC4GLM.FreqHz = FreqHz(FOI_GLMeach);
for k_ptr = 2:5 % 4.7 hour@DELL
    tic; %
    for k_sim = 1:N_sim2
        if rem(k_sim,10)==1
            disp(['-->',num2str(k_sim),'/',num2str(N_sim2),' @ cond. ',num2str(k_ptr),'/',num2str(size(prob_gens,1))])
        end
        
        data_Oxy = repmat(Observ(:,k_sim,k_ptr),[1,2])+...
            [rs_base_detrended_n{idx_sim(k_sim,1),1}(1:L*Fs,idx_sim(k_sim,2)),...
            rs_base_detrended_n{idx_sim(k_sim,3),1}(1:L*Fs,idx_sim(k_sim,4))];
        
        % Prewhitening (100)
        data_Oxy_pw100 = zeros(size(data_Oxy));
        
        % Please modify the code to use with fixed AR order
        [~, data_Oxy_pw100(:,1)] = ar_fit( data_Oxy(:,1),Pmed);
        [~, data_Oxy_pw100(:,2)] = ar_fit( data_Oxy(:,2),Pmed);
        
        synth_Oxy_pw100 = zeros(size(data_Oxy,1),size(data_Oxy,2),N_sh);
        tmp_k_sets = randperm(size(idx_sim,1));
        tmp_k_sets(tmp_k_sets == k_sim)=[];
        tmp_k_sets = tmp_k_sets(1:N_sh);
        for k_sh = 1:N_sh
            synth_Oxy = repmat(Observ(:,k_sim,k_ptr),[1,2])+...
                [rs_base_detrended_n{idx_sim(tmp_k_sets(k_sh),1),1}(1:L*Fs,idx_sim(tmp_k_sets(k_sh),2)),...
                rs_base_detrended_n{idx_sim(tmp_k_sets(k_sh),3),1}(1:L*Fs,idx_sim(tmp_k_sets(k_sh),4))];
            
            [~, synth_Oxy_pw100(:,1,k_sh)] = ar_fit_hm( synth_Oxy(:,1),Pmed,1);
            [~, synth_Oxy_pw100(:,2,k_sh)] = ar_fit_hm( synth_Oxy(:,2),Pmed,1);
        end
        
        % cut
        data_Oxy_pw100 = data_Oxy_pw100(Pmed+1:end,:);
        synth_Oxy_pw100 = synth_Oxy_pw100(Pmed+1:end,:,:);
        
        % wtc
        [tmp_wtc,wcrosspec,f]=wcoherence_modified(data_Oxy_pw100(:,1),data_Oxy_pw100(:,2),Fs, 'amor',1,12,'AR1');
        WTC4GLM.WTCpw(:,:,k_sim,k_ptr-1) = tmp_wtc(FOI_GLMeach ,:);
        
        WTCpw_amor_synth_tmp = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sh);
        parfor k_sh = 1:N_sh
            [tmp_wtc,wcrosspec,f]=wcoherence_modified(synth_Oxy_pw100(:,1,k_sh),synth_Oxy_pw100(:,2,k_sh),Fs, 'amor',1,12,'AR1');
            WTCpw_amor_synth_tmp(:,:,k_sh) = tmp_wtc(FOI_GLMeach ,:);
        end
        WTC4GLM.WTCpw_synth(:,:,k_sim,k_ptr-1) = mean(WTCpw_amor_synth_tmp,3);
        
    end
    toc;
end

WTC4GLM_SN.WTCpw = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim2,4,10);
WTC4GLM_SN.FreqHz = FreqHz(FOI_GLMeach);
WTC4GLM_SN.SNratio = ((1:10)-1)/10; % 0: no HRF
for k_ptr = 2:5
    for k_SN = 1:10
        disp(['(SN ratio set)-->',num2str(k_SN),'/10 ',num2str(k_ptr-1),'/4'])
        
        tic; %
        for k_sim = 1:N_sim2            
            data_Oxy = repmat(((k_SN-1)/10)*Observ(:,k_sim,k_ptr),[1,2])+...
                [rs_base_detrended_n{idx_sim(k_sim,1),1}(1:L*Fs,idx_sim(k_sim,2)),...
                rs_base_detrended_n{idx_sim(k_sim,3),1}(1:L*Fs,idx_sim(k_sim,4))];
            
            % Prewhitening (100)
            data_Oxy_pw100 = zeros(size(data_Oxy));
            [~, data_Oxy_pw100(:,1)] = ar_fit_hm( data_Oxy(:,1),Pmed,1);
            [~, data_Oxy_pw100(:,2)] = ar_fit_hm( data_Oxy(:,2),Pmed,1);
            
            % cut
            data_Oxy_pw100 = data_Oxy_pw100(Pmed+1:end,:);
            synth_Oxy_pw100 = synth_Oxy_pw100(Pmed+1:end,:,:);
            
            % wtc
            [tmp_wtc,wcrosspec,f]=wcoherence_modified(data_Oxy_pw100(:,1),data_Oxy_pw100(:,2),Fs, 'amor',1,12,'AR1');
            WTC4GLM_SN.WTCpw(:,:,k_sim,k_ptr-1,k_SN) = tmp_wtc(FOI_GLMeach ,:);
        end
        toc;
    end
end

%% CWT
CWT4GLM.CWTpw = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim2,2,4);
CWT4GLM.CWTpw_synth = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim2,2,4);
CWT4GLM.FreqHz = FreqHz(FOI_GLMeach);
for k_ptr = 2:5 % 4.7 hour@DELL
    tic; %
    for k_sim = 1:N_sim2
        if rem(k_sim,10)==1
            disp(['-->',num2str(k_sim),'/',num2str(N_sim2),' @ cond. ',num2str(k_ptr),'/',num2str(size(prob_gens,1))])
        end
        
        data_Oxy = repmat(Observ(:,k_sim,k_ptr),[1,2])+...
            [rs_base_detrended_n{idx_sim(k_sim,1),1}(1:L*Fs,idx_sim(k_sim,2)),...
            rs_base_detrended_n{idx_sim(k_sim,3),1}(1:L*Fs,idx_sim(k_sim,4))];
        
        % Prewhitening (100)
        data_Oxy_pw100 = zeros(size(data_Oxy));
        [~, data_Oxy_pw100(:,1)] = ar_fit_hm( data_Oxy(:,1),Pmed,1);
        [~, data_Oxy_pw100(:,2)] = ar_fit_hm( data_Oxy(:,2),Pmed,1);
        
        % 140 s
        synth_Oxy_pw100 = zeros(size(data_Oxy,1),size(data_Oxy,2),N_sh);
        tmp_k_sets = randperm(size(idx_sim,1));
        tmp_k_sets(tmp_k_sets == k_sim)=[];
        tmp_k_sets = tmp_k_sets(1:N_sh);
        for k_sh = 1:N_sh
            synth_Oxy = repmat(Observ(:,k_sim,k_ptr),[1,2])+...
                [rs_base_detrended_n{idx_sim(tmp_k_sets(k_sh),1),1}(1:L*Fs,idx_sim(tmp_k_sets(k_sh),2)),...
                rs_base_detrended_n{idx_sim(tmp_k_sets(k_sh),3),1}(1:L*Fs,idx_sim(tmp_k_sets(k_sh),4))];
            
            [~, synth_Oxy_pw100(:,1,k_sh)] = ar_fit_hm( synth_Oxy(:,1),Pmed,1);
            [~, synth_Oxy_pw100(:,2,k_sh)] = ar_fit_hm( synth_Oxy(:,2),Pmed,1);
        end

        % cut
        data_Oxy_pw100 = data_Oxy_pw100(Pmed+1:end,:);
        synth_Oxy_pw100 = synth_Oxy_pw100(Pmed+1:end,:,:);
        
        % CWT
        for k_pat = 1:2
            [tmp_cwt,F2,coi,fb] = cwt(data_Oxy_pw100(:,k_pat),'amor',Fs,'VoicesPerOctave',12);
            tmp_cwt = smoothCFS(abs(tmp_cwt).^2,fb.Scales',12);
            area_available = repmat(log2(F2),[1 size(coi2,1)]) > repmat(log2(coi2'),[size(F2,1),1]);
            tmp_cwt(~area_available)=nan;
            CWT4GLM.CWTpw(:,:,k_sim,k_pat,k_ptr-1) = tmp_cwt(FOI_GLMeach ,:);
            
            
            CWTpw_amor_synth_tmp = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sh);
            parfor k_sh = 1:N_sh
                [tmp_cwt,F2,coi,fb] = cwt(synth_Oxy_pw100(:,k_pat,k_sh),'amor',Fs,'VoicesPerOctave',12);
                tmp_cwt = smoothCFS(abs(tmp_cwt).^2,fb.Scales',12);
                area_available = repmat(log2(F2),[1 size(coi2,1)]) > repmat(log2(coi2'),[size(F2,1),1]);
                tmp_cwt(~area_available)=nan;
                CWTpw_amor_synth_tmp(:,:,k_sh) = tmp_cwt(FOI_GLMeach ,:);
            end
            CWT4GLM.CWTpw_synth(:,:,k_sim,k_pat,k_ptr-1) = mean(CWTpw_amor_synth_tmp,3);
        end
    end
    toc;
end


CWT4GLM_SN.CWTpw = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim2,2,4,10);
CWT4GLM_SN.FreqHz = FreqHz(FOI_GLMeach);
CWT4GLM_SN.SNratio = ((1:10)-1)/10; % 0: no HRF
for k_ptr = 2:5
    for k_SN = 1:10
        disp(['(SN ratio set)-->',num2str(k_SN),'/10 ',num2str(k_ptr-1),'/4'])
        
        tic; %
        for k_sim = 1:N_sim2
            data_Oxy = repmat(((k_SN-1)/10)*Observ(:,k_sim,k_ptr),[1,2])+...
                [rs_base_detrended_n{idx_sim(k_sim,1),1}(1:L*Fs,idx_sim(k_sim,2)),...
                rs_base_detrended_n{idx_sim(k_sim,3),1}(1:L*Fs,idx_sim(k_sim,4))];
            
            % Prewhitening (100)
            data_Oxy_pw100 = zeros(size(data_Oxy));
            [~, data_Oxy_pw100(:,1)] = ar_fit_hm( data_Oxy(:,1),Pmed,1);
            [~, data_Oxy_pw100(:,2)] = ar_fit_hm( data_Oxy(:,2),Pmed,1);
            
            % cut
            data_Oxy_pw100 = data_Oxy_pw100(Pmed+1:end,:);
            synth_Oxy_pw100 = synth_Oxy_pw100(Pmed+1:end,:,:);
            
            % cwt
            for k_pat = 1:2
                [tmp_cwt,F2,coi,fb] = cwt(data_Oxy_pw100(:,k_pat),'amor',Fs,'VoicesPerOctave',12);
                tmp_cwt = smoothCFS(abs(tmp_cwt).^2,fb.Scales',12);
                area_available = repmat(log2(F2),[1 size(coi2,1)]) > repmat(log2(coi2'),[size(F2,1),1]);
                tmp_cwt(~area_available)=nan;
                CWT4GLM_SN.CWTpw(:,:,k_sim,k_pat,k_ptr-1,k_SN) = tmp_cwt(FOI_GLMeach ,:);
            end
        end
        toc;
    end
end


%%
N_sh = 100;
N_sim3 = 300;

WTC4GLM_1ev.WTCpw = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim3);
WTC4GLM_1ev.WTCpw_synth = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim3);
WTC4GLM_1ev.WTCpw_noHRF = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim3);
WTC4GLM_1ev.WTCpw_ev1_p = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sim3);
WTC4GLM_1ev.FreqHz = FreqHz(FOI_GLMeach);

Observ_ev1 = zeros(L*Fs,1);
k_L = L*Fs/2;
Observ_ev1(k_L:(k_L+L_hrf-1),1) = hrf_s;

Observ_ev1_pseudo = zeros(L*Fs,300);
r_seed = repmat([1:5,-1:-1:-5],[1,30]);
r_seed = r_seed(randperm(300));
for k_s = 1:N_sim3
    k_L = L*Fs/2 + 50*r_seed(k_s);
    Observ_ev1_pseudo(k_L:(k_L+L_hrf-1),k_s) = hrf_s;
end
tic; %
for k_sim = 1:N_sim3
    if rem(k_sim,10)==1
        disp(['-->',num2str(k_sim),'/',num2str(N_sim3),' @ 1ev'])
    end
    
    data_Oxy = repmat(Observ_ev1,[1,2])+...
        [rs_base_detrended_n{idx_sim(k_sim+200,1),1}(1:L*Fs,idx_sim(k_sim+200,2)),...
        rs_base_detrended_n{idx_sim(k_sim+200,3),1}(1:L*Fs,idx_sim(k_sim+200,4))];
    data_noHRF = [rs_base_detrended_n{idx_sim(k_sim+200,1),1}(1:L*Fs,idx_sim(k_sim+200,2)),...
        rs_base_detrended_n{idx_sim(k_sim+200,3),1}(1:L*Fs,idx_sim(k_sim+200,4))];
    data_Oxy_ev1_p = repmat(Observ_ev1_pseudo(:,k_sim),[1,2])+...
        [rs_base_detrended_n{idx_sim(k_sim+200,1),1}(1:L*Fs,idx_sim(k_sim+200,2)),...
        rs_base_detrended_n{idx_sim(k_sim+200,3),1}(1:L*Fs,idx_sim(k_sim+200,4))];
    
    data_Oxy_pw100 = zeros(size(data_Oxy));
    data_noHRF_pw100 = zeros(size(data_noHRF));
    data_Oxy_ev1_p_pw100 = zeros(size(data_Oxy_ev1_p));
    [~, data_Oxy_pw100(:,1)] = ar_fit_hm( data_Oxy(:,1),Pmed,1);
    [~, data_Oxy_pw100(:,2)] = ar_fit_hm( data_Oxy(:,2),Pmed,1);
    [~, data_noHRF_pw100(:,1)] = ar_fit_hm( data_noHRF(:,1),Pmed,1);
    [~, data_noHRF_pw100(:,2)] = ar_fit_hm( data_noHRF(:,2),Pmed,1);
    [~, data_Oxy_ev1_p_pw100(:,1)] = ar_fit_hm( data_Oxy_ev1_p(:,1),Pmed,1);
    [~, data_Oxy_ev1_p_pw100(:,2)] = ar_fit_hm( data_Oxy_ev1_p(:,2),Pmed,1);
    
    % 140 s
    synth_Oxy_pw100 = zeros(size(data_Oxy,1),size(data_Oxy,2),N_sh);
    tmp_k_sets = randperm(size(idx_sim,1));
    tmp_k_sets(tmp_k_sets == (k_sim+200))=[];
    tmp_k_sets = tmp_k_sets(1:N_sh);
    for k_sh = 1:N_sh
        synth_Oxy = repmat(Observ_ev1,[1,2])+...
            [rs_base_detrended_n{idx_sim(tmp_k_sets(k_sh),1),1}(1:L*Fs,idx_sim(tmp_k_sets(k_sh),2)),...
            rs_base_detrended_n{idx_sim(tmp_k_sets(k_sh),3),1}(1:L*Fs,idx_sim(tmp_k_sets(k_sh),4))];
        
        [~, synth_Oxy_pw100(:,1,k_sh)] = ar_fit_hm( synth_Oxy(:,1),Pmed,1);
        [~, synth_Oxy_pw100(:,2,k_sh)] = ar_fit_hm( synth_Oxy(:,2),Pmed,1);
    end
    
    
    % cut
    data_Oxy_pw100 = data_Oxy_pw100(Pmed+1:end,:);
    synth_Oxy_pw100 = synth_Oxy_pw100(Pmed+1:end,:,:);
    data_noHRF_pw100 = data_noHRF_pw100(Pmed+1:end,:);
    data_Oxy_ev1_p_pw100 = data_Oxy_ev1_p_pw100(Pmed+1:end,:);
    % wtc
    [tmp_wtc,~,~]=wcoherence_modified(data_Oxy_pw100(:,1),data_Oxy_pw100(:,2),Fs, 'amor',1,12,'AR1');
    WTC4GLM_1ev.WTCpw(:,:,k_sim) = tmp_wtc(FOI_GLMeach ,:);
    [tmp_wtc,~,~]=wcoherence_modified(data_noHRF_pw100(:,1),data_noHRF_pw100(:,2),Fs, 'amor',1,12,'AR1');
    WTC4GLM_1ev.WTCpw_noHRF(:,:,k_sim) = tmp_wtc(FOI_GLMeach ,:);
    [tmp_wtc,~,~]=wcoherence_modified(data_Oxy_ev1_p_pw100(:,1),data_Oxy_ev1_p_pw100(:,2),Fs, 'amor',1,12,'AR1');
    WTC4GLM_1ev.WTCpw_ev1_p(:,:,k_sim) = tmp_wtc(FOI_GLMeach ,:);
    
    WTCpw_amor_synth_tmp = nan(size(FOI_GLMeach,2),L*Fs-Pmed,N_sh);
    parfor k_sh = 1:N_sh
        [tmp_wtc,~,~]=wcoherence_modified(synth_Oxy_pw100(:,1,k_sh),synth_Oxy_pw100(:,2,k_sh),Fs, 'amor',1,12,'AR1');
        WTCpw_amor_synth_tmp(:,:,k_sh) = tmp_wtc(FOI_GLMeach ,:);
    end
    WTC4GLM_1ev.WTCpw_synth(:,:,k_sim) = mean(WTCpw_amor_synth_tmp,3);
end
toc;


%%
corr_x = nan(size(FOI_GLMeach,2),N_sim3);
corr_x_pseudo = nan(size(FOI_GLMeach,2),N_sim3);
corr_x_noHRF = nan(size(FOI_GLMeach,2),N_sim3);
for k_sim = 1:300
    for k=1:size(WTC4GLM_1ev.WTCpw_synth,1)
        corr_x(k,k_sim) = pickcorr2(corr([WTC4GLM_1ev.WTCpw_synth(k,:,k_sim);WTC4GLM_1ev.WTCpw(k,:,k_sim)]','rows','complete'));
        corr_x_pseudo(k,k_sim) = pickcorr2(corr([WTC4GLM_1ev.WTCpw_synth(k,:,k_sim);WTC4GLM_1ev.WTCpw_ev1_p(k,:,k_sim)]','rows','complete'));
        corr_x_noHRF(k,k_sim) = pickcorr2(corr([WTC4GLM_1ev.WTCpw_synth(k,:,k_sim);WTC4GLM_1ev.WTCpw_noHRF(k,:,k_sim)]','rows','complete'));
    end
end

figure(101) % Fig 6
options.color_area = [0.4,0.4,0.9];
options.color_line = [0.1, 0.1, 0.9];
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'c95';
shadedErrorBar_auto(atanh(corr_x)', options); hold on;
options.color_area = [0.9,0.4,0.4];
options.color_line = [0.9, 0.1, 0.1];
shadedErrorBar_auto(atanh(corr_x_pseudo)', options); hold on;
options.color_area = [0.9,0.9,0.4];
options.color_line = [0.9, 0.9, 0.1];
shadedErrorBar_auto(atanh(corr_x_noHRF)', options); hold off;
ylabel('Fisher Z (a.u.)');  xlabel('Frequency (Hz)');
xticks(5:5:45); xticklabels(round(FreqHz(FOI_GLMeach(1:5:46)),3));
legend({'','Congruent','','Incongruent','','Baseline'},'Location','Best')
set(gca,'FontSize',11)


%% GLM regression
% set range of WTC to be used ffor WTC-GLM
FOI4GLM = 20:34;
% FOI4GLM = 29:34;

% WTC-GLM (pre analysis)
beta_1ev = zeros(2,N_sim3,3);
for k_sim = 1:N_sim3
    y_target = WTC4GLM_1ev.WTCpw(FOI4GLM,:,k_sim)';
    y_target =y_target(:);
    y2_target = WTC4GLM_1ev.WTCpw_ev1_p(FOI4GLM,:,k_sim)';
    y2_target =y2_target(:);
    y3_target = WTC4GLM_1ev.WTCpw_noHRF(FOI4GLM,:,k_sim)';
    y3_target =y3_target(:);
    X_regressor = WTC4GLM_1ev.WTCpw_synth(FOI4GLM,:,k_sim)';
    X_regressor =X_regressor(:);
    idx_ok = (isnan(y_target))+(isnan(X_regressor))==0;
    [beta_1ev(:,k_sim,1),dev,stats] = glmfit(X_regressor(idx_ok),y_target(idx_ok),'gamma','constant','on','link','log');
    [beta_1ev(:,k_sim,2),dev,stats] = glmfit(X_regressor(idx_ok),y2_target(idx_ok),'gamma','constant','on','link','log');
    [beta_1ev(:,k_sim,3),dev,stats] = glmfit(X_regressor(idx_ok),y3_target(idx_ok),'gamma','constant','on','link','log');
end

figure(105) % Fig 5
FOI4GLM_plot_Set = [23,43];
for k_p = 1:2 %23; %43
    subplot(2,1,k_p);
    FOI4GLM_plot = FOI4GLM_plot_Set(k_p);
    y_target = WTC4GLM_1ev.WTCpw(FOI4GLM_plot,:,k_sim)';
    y2_target = WTC4GLM_1ev.WTCpw_ev1_p(FOI4GLM_plot,:,k_sim)';
    y3_target = WTC4GLM_1ev.WTCpw_noHRF(FOI4GLM_plot,:,k_sim)';
    X_regressor = WTC4GLM_1ev.WTCpw_synth(FOI4GLM_plot,:,k_sim)';
    % idx_ok = (isnan(y_target))+(isnan(X_regressor))==0;
    idx_ok = ones(size(y_target),'logical');
    plot(y_target(idx_ok),'-','color',0.3*ones(1,3),'LineWidth',1.5); hold on;
    plot(y3_target(idx_ok),'--','color',0.1*ones(1,3),'LineWidth',1);
    plot(X_regressor(idx_ok),'-','color',0.8*ones(1,3),'LineWidth',3); hold off;
    xticks(400:500:2900); xlim([1 2900]); xticklabels(50:50:300); xlabel('Time (s)')
    ylabel('WTC (a.u.)'); ylim([0 0.8]);
    legend({'Target','Baseline','Regressor'},'Location','Best')
    title([num2str(round(WTC4GLM_1ev.FreqHz(FOI4GLM_plot),3)),' Hz'])
    set(gca,'FontSize',11)
end

%% WTC-GLM (main part)
beta_cong = zeros(2,N_sim2,4);
beta_incong = zeros(2,(N_sim2-1),N_sim2,4);
beta_cong_SN = zeros(2,N_sim2,4,10);
beta_incong_SN = zeros(2,(N_sim2-1),N_sim2,4,10);
dev_cong = zeros(N_sim2,4);
dev_incong = zeros((N_sim2-1),N_sim2,4);
dev_cong_SN = zeros(N_sim2,4,10);
for k_ptri = 1:4
    for k_sim = 1:N_sim2
        if rem(k_sim,10)==1
            disp(['++',num2str(k_sim),'/',num2str(N_sim2),' @ ',num2str(k_ptri)]);
        end
        y_target = WTC4GLM.WTCpw(FOI4GLM,:,k_sim,k_ptri)';
        y_target =y_target(:);
        X_regressor = WTC4GLM.WTCpw_synth(FOI4GLM,:,k_sim,k_ptri)';
        X_regressor =X_regressor(:)-min(X_regressor(:));
        
        y_target_SN = permute(WTC4GLM_SN.WTCpw(FOI4GLM,:,k_sim,k_ptri,:),[2,1,5,3,4]);
        y_target_SN = reshape(y_target_SN,[size(y_target_SN,1)*size(y_target_SN,2),size(y_target_SN,3)]);
        
        idx_ok = (isnan(y_target))+(isnan(X_regressor))==0;
        
        % true HRF x true data (congruent)
        [beta_cong(:,k_sim,k_ptri),dev_cong(k_sim,k_ptri),stats] = glmfit(X_regressor(idx_ok),y_target(idx_ok),'gamma','constant','on','link','log');
        
        % pseudo HRF x true data (incongruent)
        tmp_idx = randperm(N_sim2);
        tmp_idx(tmp_idx==k_sim) = [];
        y_pseudo_SN = permute(WTC4GLM_SN.WTCpw(FOI4GLM,:,:,k_ptri,:),[2,1,5,3,4]);
        y_pseudo_SN = reshape(y_pseudo_SN,[size(y_pseudo_SN,1)*size(y_pseudo_SN,2),size(y_pseudo_SN,3),size(y_pseudo_SN,4)]);
        for k_p = 1:(N_sim2-1)
            y_pseudo = WTC4GLM.WTCpw(FOI4GLM,:,tmp_idx(k_p),k_ptri)';
            y_pseudo =y_pseudo(:);
            [beta_incong(:,k_p,k_sim,k_ptri),dev_incong(k_p,k_sim,k_ptri),stats] = glmfit(X_regressor(idx_ok),y_pseudo(idx_ok),'gamma','constant','on','link','log');
        end
        
        for k_SN = 1:10
            [beta_cong_SN(:,k_sim,k_ptri,k_SN),dev_cong_SN(k_sim,k_ptri,k_SN),stats] = ...
                glmfit(X_regressor(idx_ok),y_target_SN(idx_ok,k_SN),'gamma','constant','on','link','log');
            for k_p = 1:(N_sim2-1)
                y_pseudo = y_pseudo_SN(:,k_SN,tmp_idx(k_p));
                beta_incong_SN(:,k_p,k_sim,k_ptri,k_SN) = glmfit(X_regressor(idx_ok),y_pseudo(idx_ok),'gamma','constant','on','link','log');
            end
        end
    end
end

betaWTC.cong = beta_cong;
betaWTC.cong_SN = beta_cong_SN;
betaWTC.incong = beta_incong;
betaWTC.incong_SN = beta_incong_SN;


%% CWT-GLM
betaCWT_cong = zeros(2,N_sim2,2,4);
betaCWT_incong = zeros(2,(N_sim2-1),N_sim2,2,4);
betaCWT_cong_SN = zeros(2,N_sim2,2,4,10);
betaCWT_incong_SN = zeros(2,(N_sim2-1),N_sim2,2,4,10);
devCWT_cong = zeros(N_sim2,2,4);
devCWT_incong = zeros((N_sim2-1),N_sim2,2,4);
devCWT_cong_SN = zeros(N_sim2,2,4,10);
for k_ptri = 1:4
    for k_sim = 1:N_sim2
        if rem(k_sim,10)==1
            disp(['++',num2str(k_sim),'/',num2str(N_sim2),' @ ',num2str(k_ptri)]);
        end
        
        for k_pat = 1:2
            y_target = CWT4GLM.CWTpw(FOI4GLM,:,k_sim,k_pat,k_ptri)';
            y_target = y_target(:);
            X_regressor = CWT4GLM.CWTpw_synth(FOI4GLM,:,k_sim,k_pat,k_ptri)';
            X_regressor =X_regressor(:)-min(X_regressor(:));
            X_regressor = X_regressor;
            
            y_target_SN = permute(CWT4GLM_SN.CWTpw(FOI4GLM,:,k_sim,k_pat,k_ptri,:),[2,1,6,3,4,5]);
            y_target_SN = reshape(y_target_SN,[size(y_target_SN,1)*size(y_target_SN,2),size(y_target_SN,3)]);
            y_target_SN = y_target_SN;
            
            idx_ok = (isnan(y_target))+(isnan(X_regressor))==0;
            
            % true HRF x true data (congruent)
            [betaCWT_cong(:,k_sim,k_pat,k_ptri),devCWT_cong(k_sim,k_pat,k_ptri),stats] = glmfit(10^15*X_regressor(idx_ok),10^15*y_target(idx_ok),'gamma','constant','on','link','log');
            
            % pseudo HRF x true data (incongruent)
            tmp_idx = randperm(N_sim2);
            tmp_idx(tmp_idx==k_sim) = [];
            y_pseudo_SN = permute(CWT4GLM_SN.CWTpw(FOI4GLM,:,:,k_pat,k_ptri,:),[2,1,6,3,4,5]);
            y_pseudo_SN = reshape(y_pseudo_SN,[size(y_pseudo_SN,1)*size(y_pseudo_SN,2),size(y_pseudo_SN,3),size(y_pseudo_SN,4)]);
            
            for k_p = 1:(N_sim2-1)
                y_pseudo = CWT4GLM.CWTpw(FOI4GLM,:,tmp_idx(k_p),k_ptri)';
                y_pseudo = y_pseudo(:);
                [betaCWT_incong(:,k_p,k_sim,k_pat,k_ptri),devCWT_incong(k_p,k_sim,k_pat,k_ptri),stats] = glmfit(10^15*X_regressor(idx_ok),10^15*y_pseudo(idx_ok),'gamma','constant','on','link','log');
            end

            for k_SN = 1:10
                [betaCWT_cong_SN(:,k_sim,k_pat,k_ptri,k_SN),devCWT_cong_SN(k_sim,k_pat,k_ptri,k_SN),stats] = ...
                    glmfit(10^15*X_regressor(idx_ok),10^15*y_target_SN(idx_ok,k_SN),'gamma','constant','on','link','log');
                parfor k_p = 1:(N_sim2-1)
                    y_pseudo = y_pseudo_SN(:,k_SN,tmp_idx(k_p));
                    betaCWT_incong_SN(:,k_p,k_sim,k_pat,k_ptri,k_SN) = glmfit(10^15*X_regressor(idx_ok),10^15*y_pseudo(idx_ok),'gamma','constant','on','link','log');
                end
            end
        end
    end
end

betaCWT.cong = betaCWT_cong;
betaCWT.cong_SN = betaCWT_cong_SN;
betaCWT.incong = betaCWT_incong;
betaCWT.incong_SN = betaCWT_incong_SN;


%% Cohen_d %
Cohen_d = zeros(4,1);
for k_ptri = 1:4
    [Cohen_d(k_ptri),~] = func_hedges_g(beta_cong(2,:,k_ptri),extm(beta_cong_SN(2,:,k_ptri,1)));
    func_est_sample2cohenDTtest(beta_cong(2,:,k_ptri),extm(beta_cong_SN(2,:,k_ptri,1)),20,100);
end

Cohen_d_SN = zeros(10,4);
for k_ptri = 1:4
    for k_SN = 1:9
        Cohen_d_SN(k_SN,k_ptri) = func_hedges_g(extm(beta_cong_SN(2,:,k_ptri,k_SN+1)),extm(beta_cong_SN(2,:,k_ptri,1)));
    end
    Cohen_d_SN(10,k_ptri) = func_hedges_g(beta_cong(2,:,k_ptri),extm(beta_cong_SN(2,:,k_ptri,1)));
end

Cohen_d_SN_incong = zeros(10,4);
for k_ptri = 1:4
    for k_SN = 1:9
        Cohen_d_SN_incong(k_SN,k_ptri) = func_hedges_g(extm(beta_cong_SN(2,:,k_ptri,k_SN+1)),extm(beta_incong_SN(2,:,:,k_ptri,k_SN+1)));
    end
    Cohen_d_SN_incong(10,k_ptri) = func_hedges_g(extm(beta_cong(2,:,k_ptri)),extm(beta_incong(2,:,:,k_ptri)));
end

Cohen_d_CWT = zeros(4,1);
for k_ptri = 1:4
    [Cohen_d_CWT(k_ptri),~] = func_hedges_g(extm(betaCWT_cong(n_coeff,:,:,k_ptri)),extm(betaCWT_cong_SN(n_coeff,:,:,k_ptri,1)));
    func_est_sample2cohenDTtest(extm(betaCWT_cong(n_coeff,:,:,k_ptri)),extm(betaCWT_cong_SN(n_coeff,:,:,k_ptri,1)),20,100);
end

Cohen_d_SN_CWT = zeros(10,4);
for k_ptri = 1:4
    for k_SN = 1:9
        Cohen_d_SN_CWT(k_SN,k_ptri) = func_hedges_g(extm( cat(2,betaCWT_cong_SN(n_coeff,:,1,k_ptri,k_SN+1),betaCWT_cong_SN(n_coeff,:,2,k_ptri,k_SN+1)) ),...
            extm( cat(2,betaCWT_cong_SN(n_coeff,:,1,k_ptri,1),betaCWT_cong_SN(n_coeff,:,2,k_ptri,1))   ));        
    end
    Cohen_d_SN_CWT(10,k_ptri) = func_hedges_g(extm( cat(2,betaCWT_cong(n_coeff,:,1,k_ptri),betaCWT_cong(n_coeff,:,2,k_ptri)) ),...
        extm( cat(2,betaCWT_cong_SN(n_coeff,:,1,k_ptri,1),betaCWT_cong_SN(n_coeff,:,2,k_ptri,1)) ));
end


betaCWT_incong(abs(betaCWT_incong)>100) = nan;
% VS incongruent
Cohen_d_SN_CWT_incong = zeros(10,4);
for k_ptri = 1:4
    for k_SN = 1:9
        Cohen_d_SN_CWT_incong(k_SN,k_ptri) = func_hedges_g(extm( betaCWT_cong_SN(n_coeff,:,:,k_ptri,k_SN+1)),...
            extm( betaCWT_incong_SN(n_coeff,:,:,:,k_ptri,k_SN+1)));        
    end   
    Cohen_d_SN_CWT_incong(10,k_ptri) = func_hedges_g(extm( betaCWT_cong(n_coeff,:,:,k_ptri)),...
            extm( betaCWT_incong(n_coeff,:,:,:,k_ptri)));     
end



figure(206) % Fig S7
k_ptri = 1;
beta_minmax_WTC = [min(extm(beta_cong_SN(2,:,k_ptri,:))),max(extm(beta_cong_SN(2,:,k_ptri,:)));
    min(extm(beta_cong(2,:,k_ptri))),max(extm(beta_cong(2,:,k_ptri)))];
beta_minmax_WTC = max(abs([floor(min(beta_minmax_WTC(:,1))),ceil(max(beta_minmax_WTC(:,2)))]));

beta_minmax_CWT = [min(extm(betaCWT.cong_SN(n_coeff,:,:,k_ptri,:))),max(extm(betaCWT.cong_SN(n_coeff,:,:,k_ptri,:)));
    min(extm(betaCWT.cong(n_coeff,:,:,k_ptri))),max(extm(betaCWT.cong(n_coeff,:,:,k_ptri)))];
beta_minmax_CWT = max(abs([floor(min(beta_minmax_CWT(:,1))),ceil(max(beta_minmax_CWT(:,2)))]));

subplot(2,1,1)
histogram(beta_cong(2,:,k_ptri),[-beta_minmax_WTC:1:beta_minmax_WTC],'FaceColor',[0,0,1],'FaceAlpha',0.5,'Normalization','probability'); hold on;
histogram(beta_cong_SN(2,:,k_ptri,1),[-beta_minmax_WTC:1:beta_minmax_WTC],'FaceColor',[1,0,0],'FaceAlpha',0.4,'Normalization','probability');
hold off;
ylabel('Probability'); xlabel('\beta'); ylim([0 0.21]);
title(['prob = ',num2str(prob_gens(k_ptri+1)),' | Cohen d = ',num2str(round(Cohen_d(k_ptri),2))]); set(gca,'FontSize',11);
n_coeff = 2;
[Cohen_d_CWT_tmp,~] = func_hedges_g(extm(betaCWT.cong(n_coeff,:,:,k_ptri)),extm(betaCWT.cong_SN(n_coeff,:,:,k_ptri,1)));

subplot(2,1,2)
histogram(extm(betaCWT.cong(n_coeff,:,:,k_ptri)),[-beta_minmax_CWT:0.25:beta_minmax_CWT],'FaceColor',[0,0,1],'FaceAlpha',0.5,'Normalization','probability'); hold on;
histogram(extm(betaCWT.cong_SN(n_coeff,:,:,k_ptri,1)),[-beta_minmax_CWT:0.25:beta_minmax_CWT],'FaceColor',[1,0,0],'FaceAlpha',0.4,'Normalization','probability'); 
hold off;
legend({'Congruent','Rest'},'Location','best');
ylabel('Probability'); xlabel('\beta'); ylim([0 0.42]);
title(['prob = ',num2str(prob_gens(k_ptri+1)),' | Cohen d = ',num2str(round(Cohen_d_CWT_tmp,2))]); set(gca,'FontSize',11);


% to draw Fig 7 and Fig S8, combine the results of different settings of FOI4GLM