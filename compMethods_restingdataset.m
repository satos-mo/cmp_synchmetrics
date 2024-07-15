close all

%% ADDPATH
% Please add path to ...
%  Homer2 (hmrOD2Conc.m, hmrBandpassFilt.m) from https://homer-fnirs.org/
%  NIRS-SPM (detrend_wavelet_MDL.m) from https://www.nitrc.org/projects/nirs_spm/
%  huppertt-nirs-toolbox (ar_fit.m) from https://github.com/huppertt/nirs-toolbox
%   to fix the AR order, modify the function


%% Load File
%%% SET Directories of datasets
dir_name = '****\s045\NIRS\'; % change

type = {'Rest','AV','Face','Lang'};

%% Set paramters
Fs = 4; % 4 Hz

ttype = 1; % Rest
fname = [dir_name,type{1,ttype},'\allruns.snirf'];
data_raw=loadh5(fname);

%% Convert data (raw to od) using function in Homer2 (hmrOD2Conc) or other function
data_od = OD2Conc( data_raw ); % This OD2Conc is dummy

%% pre-processing
k_subj=1;%:size(data_od,1)

% HDMS
%  Please apply hdms using hdms_sep.p & POTETo
%   https://github.com/hkwgc/open-potato
%   https://unit.aist.go.jp/hiiri/nrehrg/download/dl001.html
% tmp_data <- data_od{k_subj}

% Remove trend by Wavelet-MDL
bias_Oxy = detrend_wavelet_MDL(tmp_data (:,:,1),ones(size(tmp_data,1),1));
bias_Deoxy = detrend_wavelet_MDL(tmp_data (:,:,2),ones(size(tmp_data,1),1));
tmp_data(:,:,1) = tmp_data (:,:,1)-bias_Oxy;
tmp_data(:,:,2) = tmp_data (:,:,2)-bias_Deoxy;
tmp_data(:,:,3) = tmp_data(:,:,1)+tmp_data(:,:,2);

% Pre-whitening
Pmed = 10*Fs; % 10 s
tmp_data_Oxy_pw = nan(size(tmp_data,1),size(tmp_data,2));
for k_ch = 1:size(tmp_data,2)
    [~, tmp_data_Oxy_pw(:,k_ch)] = ar_fit(tmp_data(:,k_ch),Pmed);
end
tmp_data_Oxy_pw = tmp_data_Oxy_pw(Pmed+1:end,:);

%%
Comb_ch = nchoosek( 1:size(tmp_data,2) , 2 );




%% Band-pass (resting)
lpf = 0.08;
hpf = 0.009;

y2 = hmrBandpassFilt( tmp_data(:,:,1), Fs, hpf, lpf );
y2_pw = hmrBandpassFilt(  tmp_data_Oxy_pw , Fs, hpf, lpf );

% cut
tmp_data_Oxy = tmp_data(Pmed+1:end,:,1);
y2 = y2(Pmed+1:end,:);

%% Correlation
tmp_corr = corrcoef(y2);
FC = nan(1,size(Comb_ch,1));
for k_ch_comb = 1:size(Comb_ch,1)
    k_ch_A = Comb_ch(k_ch_comb,1);
    k_ch_B = Comb_ch(k_ch_comb,2);
    FC(1,k_ch_comb) = tmp_corr(k_ch_A,k_ch_B);
end
tmp_corr = corrcoef(y2_pw);
FC_pw = nan(1,size(Comb_ch,1));
for k_ch_comb = 1:size(Comb_ch,1)
    k_ch_A = Comb_ch(k_ch_comb,1);
    k_ch_B = Comb_ch(k_ch_comb,2);
    FC_pw(1,k_ch_comb) = tmp_corr(k_ch_A,k_ch_B);
end

%% Dynamic Correlation
w_DFC = [20,80,120,240];
for k_w = 1:size(w_DFC,2)
    dur = w_DFC(k_w);
    
    DFC.(['w',num2str(w_DFC(k_w))]) = nan(size(y2,1),size(Comb_ch,1));
    
    for k=1:(size(y2,1)-dur+1)
        if rem(k,500) == 1
            disp(['DFC w',num2str(w_DFC(k_w)),': ',num2str(round(k/((size(y2,1)-dur+1))*100)),'%'])
        end
        
        tmp_corr = corr(y2(k:k+dur-1,:));
        for k_ch_comb = 1:size(Comb_ch,1)
            k_ch_A = Comb_ch(k_ch_comb,1);
            k_ch_B = Comb_ch(k_ch_comb,2);
            DFC.(['w',num2str(w_DFC(k_w))])(k+floor(dur/2),k_ch_comb) = tmp_corr(k_ch_A,k_ch_B);
        end
    end
end


%% Dynamic Correlation w/ prewhitening

% Dynamic FC
w_DFC_pw = [20,80,120,240];
for k_w = 1:size(w_DFC_pw,2)
    dur = w_DFC_pw(k_w);
    
    DFC_pw.(['w',num2str(w_DFC_pw(k_w))]) = nan(size(y2,1),size(Comb_ch,1));
    
    for k=1:(size(y2_pw,1)-dur+1)
        if rem(k,500) == 1
            disp(['DFC w',num2str(w_DFC_pw(k_w)),': ',num2str(round(k/((size(y2,1)-dur+1))*100)),'%'])
        end
        
        tmp_corr = corr(y2_pw(k:k+dur-1,:));
        for k_ch_comb = 1:size(Comb_ch,1)
            k_ch_A = Comb_ch(k_ch_comb,1);
            k_ch_B = Comb_ch(k_ch_comb,2);
            DFC_pw.(['w',num2str(w_DFC_pw(k_w))])(k+floor(dur/2),k_ch_comb) = tmp_corr(k_ch_A,k_ch_B);
        end
    end
end



%% Wavelet coherence
N_band = 75;
FOI = 12:68;


WTC_amor = nan(N_band,size(tmp_data_Oxy_pw,1),size(Comb_ch,1));
WTC_amor_meanabs = nan(size(FOI,2),size(Comb_ch,1));
WTCpw_amor = nan(N_band,size(tmp_data_Oxy_pw,1),size(Comb_ch,1));
WTCpw_amor_meanabs = nan(size(FOI,2),size(Comb_ch,1));
WTCpw_cgau1 = nan(N_band,size(tmp_data_Oxy_pw,1),size(Comb_ch,1));
WTCpw_cgau1_meanabs = nan(size(FOI,2),size(Comb_ch,1));
WTCpw_cgau2 = nan(N_band,size(tmp_data_Oxy_pw,1),size(Comb_ch,1));
WTCpw_cgau2_meanabs = nan(size(FOI,2),size(Comb_ch,1));

WTCpw_amor_ci95_upper = nan(N_band,size(tmp_data_Oxy_pw,1),size(Comb_ch,1));
WTCpw_amor_ci95_lower = nan(N_band,size(tmp_data_Oxy_pw,1),size(Comb_ch,1));
WTCpw_amor_meanabs_ci95_upper = nan(size(FOI,2),size(Comb_ch,1));
WTCpw_amor_meanabs_randmed = nan(size(FOI,2),size(Comb_ch,1));
tic;
for k_ch_comb = 1:size(Comb_ch,1)
    if rem(k_ch_comb,100) == 1
        disp(['*: ',num2str(round(k_ch_comb/size(Comb_ch,1)*100)),'%'])
    end
    k_ch_A = Comb_ch(k_ch_comb,1);
    k_ch_B = Comb_ch(k_ch_comb,2);
    
    [tmp_wtc,wcrosspec]=wcoherence_modified(tmp_data_Oxy(:,k_ch_A,1),tmp_data_Oxy(:,k_ch_B,1),Fs, 'amor',1,12);
    WTC_amor(:,:,k_ch_comb) = tmp_wtc(1:N_band,:);
    WTC_amor_meanabs(:,k_ch_comb) = abs(mean(wcrosspec(FOI,sum(isnan(wcrosspec(FOI,:)),1)==0),2));
    [tmp_wtc,wcrosspec,~,wtc_rp,wcrosspec_rp]=wcoherence_modified(tmp_data_Oxy_pw(:,k_ch_A,1),tmp_data_Oxy_pw(:,k_ch_B,1),Fs, 'amor',1000,12);
    WTCpw_amor(:,:,k_ch_comb) = tmp_wtc(1:N_band,:);
    WTCpw_amor_meanabs(:,k_ch_comb) = abs(mean(wcrosspec(FOI,sum(isnan(wcrosspec(FOI,:)),1)==0),2));
    WTCpw_amor_ci95_upper(:,:,k_ch_comb) = prctile(wtc_rp(1:N_band,:,:), 97.5, 3);
    WTCpw_amor_ci95_lower(:,:,k_ch_comb) = prctile(wtc_rp(1:N_band,:,:), 2.5, 3);
    WTCpw_amor_meanabs_ci95_upper(:,k_ch_comb) = prctile(abs(mean(wcrosspec_rp(FOI,sum(isnan(wcrosspec(FOI,:)),1)==0,:),2)), 97.5, 3);
    WTCpw_amor_meanabs_randmed(:,k_ch_comb) = median(abs(mean(wcrosspec_rp(FOI,sum(isnan(wcrosspec(FOI,:)),1)==0,:),2)), 3);
    
    [tmp_wtc,wcrosspec]=wcoherence_modified(tmp_data_Oxy_pw(:,k_ch_A,1),tmp_data_Oxy_pw(:,k_ch_B,1),Fs, 'cgau1',1,12);
    WTCpw_cgau1(:,:,k_ch_comb) = tmp_wtc(1:N_band,:);
    WTCpw_cgau1_meanabs(:,k_ch_comb) = abs(mean(wcrosspec(FOI,sum(isnan(wcrosspec(FOI,:)),1)==0),2));
    [tmp_wtc,wcrosspec]=wcoherence_modified(tmp_data_Oxy_pw(:,k_ch_A,1),tmp_data_Oxy_pw(:,k_ch_B,1),Fs, 'cgau2',1,12);
    WTCpw_cgau2(:,:,k_ch_comb) = tmp_wtc(1:N_band,:);
    WTCpw_cgau2_meanabs(:,k_ch_comb) = abs(mean(wcrosspec(FOI,sum(isnan(wcrosspec(FOI,:)),1)==0),2));%
end
[~,~,f1]=wcoherence_modified(y2(:,1,1),y2(:,2,1),Fs, 'amor',1,12);

WTC.amor = WTC_amor;
WTC.amor_meanabs = WTC_amor_meanabs;

WTC.amor_pw = WTCpw_amor;
WTC.amor_pw_meanabs = WTCpw_amor_meanabs;
WTC.amor_pw_ci95_upper= WTCpw_amor_ci95_upper;
WTC.amor_pw_ci95_lower = WTCpw_amor_ci95_lower;
WTC.amor_pw_meanabs_ci95_upper = WTCpw_amor_meanabs_ci95_upper;
WTC.amor_pw_meanabs_randmed = WTCpw_amor_meanabs_randmed;

WTC.cgau1_pw = WTCpw_cgau1;
WTC.cgau1_pw_meanabs = WTCpw_cgau1_meanabs;
WTC.cgau2_pw = WTCpw_cgau2;
WTC.cgau2_pw_meanabs = WTCpw_cgau2_meanabs;
WTC.Freq = f1(1:N_band);
WTC.FOI_meanabs = f1(FOI);
toc;


%% phase mutual information (w/o prewhitening)

pMIseq.MI = nan(1,size(Comb_ch,1)  );
pMIseq.w120 = nan(size(y2,1),size(Comb_ch,1)  );
pMIseq.w240 = nan(size(y2,1),size(Comb_ch,1)  );
for k_ch_comb = 1:size(Comb_ch,1)
    if rem(k_ch_comb,100) == 1
        disp(['*: ',num2str(round(k_ch_comb/size(Comb_ch,1)*100)),'%'])
    end
    k_ch_A = Comb_ch(k_ch_comb,1);
    k_ch_B = Comb_ch(k_ch_comb,2);
    
    
    hilbert_A = func_hilbert_nonnan(y2(:,k_ch_A,1));
    hilbert_B = func_hilbert_nonnan(y2(:,k_ch_B,1));
    phase_data1 = angle(hilbert_A);
    phase_data1 = phase_data1+pi;
    phase_data2 = angle(hilbert_B);
    phase_data2 = phase_data2+pi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % binsize as in Scott et al. %
    binsize1 = 3.49* std(phase_data1,[],1,'omitnan')* sum(~isnan(phase_data1))^(-1/3);
    binsize2 = 3.49* std(phase_data2,[],1,'omitnan')* sum(~isnan(phase_data2))^(-1/3);
    Nbins1 = length(0:binsize1:2*pi);
    Nbins2 = length(0:binsize2:2*pi);
    binsize1 = 2*pi/Nbins1; % update
    binsize2 = 2*pi/Nbins2; % update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_th1 = 0:2*pi/Nbins1:2*pi;
    state_data1 = sum(repmat(phase_data1,[1,Nbins1+1]) >=  repmat(tmp_th1,[size(phase_data1,1),1]),2);
    tmp_th2 = 0:2*pi/Nbins2:2*pi;
    state_data2 = sum(repmat(phase_data2,[1,Nbins2+1]) >=  repmat(tmp_th2,[size(phase_data2,1),1]),2);
    
    pMIseq.MI(1,k_ch_comb) = func_seqNMI_cluster(state_data1,state_data2,size(state_data1,1),size(state_data1,1));
    MIseq_30 = func_seqNMI_cluster(state_data1,state_data2,Fs*30,1,100);
    pMIseq.w120(ceil(Fs*30/2):ceil(Fs*30/2)+size(MIseq_30,1)-1,k_ch_comb) = MIseq_30;
    MIseq_60 = func_seqNMI_cluster(state_data1,state_data2,Fs*60,1,100);
    pMIseq.w240(ceil(Fs*60/2):ceil(Fs*60/2)+size(MIseq_60,1)-1,k_ch_comb) = MIseq_60;
end




%% phase mutual information (w/ prewhitening)

pMIseq_pw.MI = nan(1,size(Comb_ch,1)  );
pMIseq_pw.w120 = nan(size(y2_pw,1),size(Comb_ch,1)  );
pMIseq_pw.w240 = nan(size(y2_pw,1),size(Comb_ch,1)  );

for k_ch_comb = 1:size(Comb_ch,1)
    if rem(k_ch_comb,100) == 1
        disp(['*: ',num2str(round(k_ch_comb/size(Comb_ch,1)*100)),'%'])
    end
    k_ch_A = Comb_ch(k_ch_comb,1);
    k_ch_B = Comb_ch(k_ch_comb,2);
    
    
    hilbert_A = func_hilbert_nonnan(y2_pw(:,k_ch_A,1));
    hilbert_B = func_hilbert_nonnan(y2_pw(:,k_ch_B,1));
    phase_data1 = angle(hilbert_A);
    phase_data1 = phase_data1+pi;
    phase_data2 = angle(hilbert_B);
    phase_data2 = phase_data2+pi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % binsize as in Scott et al. %
    % Scott, D. W. (1992). Multivariate density estimation: theory, practice, and visualization. John Wiley & Sons.
    % https://doi.org/10.1002/9780470316849
    binsize1 = 3.49* std(phase_data1,[],1,'omitnan')* sum(~isnan(phase_data1))^(-1/3);
    binsize2 = 3.49* std(phase_data2,[],1,'omitnan')* sum(~isnan(phase_data2))^(-1/3);
    Nbins1 = length(0:binsize1:2*pi);
    Nbins2 = length(0:binsize2:2*pi);
    binsize1 = 2*pi/Nbins1; % update
    binsize2 = 2*pi/Nbins2; % update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_th1 = 0:2*pi/Nbins1:2*pi;
    state_data1 = sum(repmat(phase_data1,[1,Nbins1+1]) >=  repmat(tmp_th1,[size(phase_data1,1),1]),2);
    tmp_th2 = 0:2*pi/Nbins2:2*pi;
    state_data2 = sum(repmat(phase_data2,[1,Nbins2+1]) >=  repmat(tmp_th2,[size(phase_data2,1),1]),2);
    
    pMIseq_pw.MI(1,k_ch_comb) = func_seqNMI_cluster(state_data1,state_data2,size(state_data1,1),size(state_data1,1));
    
    MIseq_30_pw = func_seqNMI_cluster(state_data1,state_data2,Fs*30,1,100);
    pMIseq_pw.w120(ceil(Fs*30/2):ceil(Fs*30/2)+size(MIseq_30_pw,1)-1,k_ch_comb) = MIseq_30_pw;
    MIseq_60_pw = func_seqNMI_cluster(state_data1,state_data2,Fs*60,1,100);
    pMIseq_pw.w240(ceil(Fs*60/2):ceil(Fs*60/2)+size(MIseq_60_pw,1)-1,k_ch_comb) = MIseq_60_pw;
end



FOI_pMIwt = 51:4:75;
% phase Mutual Information based on WT

pMIwtseq.MI = nan(size(tmp_data,1),size(Comb_ch,1),size(FOI_pMIwt,2));  % fixed

[~,F1,~] = cwt(tmp_data(:,1,1),'amor',Fs,'VoicesPerOctave',12);
c_delay = round(Fs./(F1(FOI_pMIwt)'));
pMIwtseq.c_delay = c_delay;
pMIwtseq.FOI = F1(FOI_pMIwt);

tic; % 18.7 hours
for k_ch_comb = 1:size(Comb_ch,1)
    if rem(k_ch_comb,100) == 1
        disp(['MI @ WT *: ',num2str(round(k_ch_comb/size(Comb_ch,1)*100)),'%'])
    end
    k_ch_A = Comb_ch(k_ch_comb,1);
    k_ch_B = Comb_ch(k_ch_comb,2);
    
    [DA,F1,coi1] = cwt(tmp_data(:,k_ch_A,1),'amor',Fs,'VoicesPerOctave',12);
    DB = cwt(tmp_data(:,k_ch_B,1),'amor',Fs,'VoicesPerOctave',12);
    %area_available = repmat(log2(F1),[1 size(coi1,2)]) > repmat(log2(coi1),[size(F1,1),1]);
    area_available = repmat(log2(F1),[1 size(coi1,1)]) > repmat(log2(coi1)',[size(F1,1),1]);
    DA(~area_available)=nan;
    DB(~area_available)=nan;
    
    phase_data1 = angle(DA(FOI_pMIwt,:)');
    phase_data1 = phase_data1+pi;
    phase_data2 = angle(DB(FOI_pMIwt,:)');
    phase_data2 = phase_data2+pi;
    
    N_ignore_COI = sum(~area_available(FOI_pMIwt,:),2)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % binsize as in Scott et al. %
    binsize1 = 3.49.* std(phase_data1,[],1,'omitnan').* sum(~isnan(phase_data1),1).^(-1/3);
    binsize2 = 3.49.* std(phase_data2,[],1,'omitnan').* sum(~isnan(phase_data2),1).^(-1/3);
    Nbins = zeros(2,size(FOI_pMIwt,2));
    for k_N = 1:size(FOI_pMIwt,2)
        Nbins(1,k_N) = length(0:binsize1(k_N):2*pi);
        Nbins(2,k_N) = length(0:binsize2(k_N):2*pi);
    end
    binsize = 2.*pi./Nbins; % update
    
    
    for k_N = 1:size(FOI_pMIwt,2)
        tmp_th1 = 0:2*pi/Nbins(1,k_N):2*pi;
        state_data1 = sum(repmat(phase_data1(:,k_N),[1,Nbins(1,k_N)+1]) >=  repmat(tmp_th1,[size(phase_data1(:,k_N),1),1]),2);
        tmp_th2 = 0:2*pi/Nbins(2,k_N):2*pi;
        state_data2 = sum(repmat(phase_data2(:,k_N),[1,Nbins(2,k_N)+1]) >=  repmat(tmp_th2,[size(phase_data2(:,k_N),1),1]),2);
        state_data1(state_data1==0)=nan;
        state_data2(state_data2==0)=nan;
        %            [MIseq_tmp,MIseq_rand_tmp,MIseq_max] = func_seqNMI_cluster(state_data1,state_data2,c_delay(k_N),1,100);
        %             pMIwtseq.NMI(ceil(c_delay(k_N)/2):ceil(c_delay(k_N)/2)+size(MIseq_tmp,1)-1,k_ch_comb,k_N) = (MIseq_tmp-mean(MIseq_rand_tmp,2))./MIseq_max;
        MIseq_tmp = func_seqNMI_cluster(state_data1,state_data2,c_delay(k_N),1,100);
        MIseq_tmp(1:N_ignore_COI(k_N)) = nan;
        MIseq_tmp(end-N_ignore_COI(k_N)+1:end) = nan;
        pMIwtseq.MI(ceil(c_delay(k_N)/2):ceil(c_delay(k_N)/2)+size(MIseq_tmp,1)-1,k_ch_comb,k_N) = MIseq_tmp;
    end
end


%%

figure(12)
freqsets = [17,41,49,57];
for k_f = 1:size(freqsets,2)
    subplot(size(freqsets,2),4,1 + 4*(k_f-1));
    plot(WTC.amor_pw_meanabs(freqsets(k_f),:),WTC.amor_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel(['amor at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    ylabel(['amor w/o pw at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']); set(gca,'FontSize',11);
    pbaspect([1 1 1]);
    
    subplot(size(freqsets,2),4,2 + 4*(k_f-1));
    plot(WTC.amor_pw_meanabs(freqsets(k_f),:),WTC.cgau1_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel(['amor at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    ylabel(['cgau1 at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']); set(gca,'FontSize',11);
    pbaspect([1 1 1]);
    subplot(size(freqsets,2),4,3 + 4*(k_f-1));
    plot(WTC.amor_pw_meanabs(freqsets(k_f),:),WTC.cgau2_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel(['amor at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    ylabel(['cgau2 at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']); set(gca,'FontSize',11);
    pbaspect([1 1 1]);
    subplot(size(freqsets,2),4,4*k_f);
    plot(WTC.cgau1_pw_meanabs(freqsets(k_f),:),WTC.cgau2_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel(['cgau1 at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    ylabel(['cgau2 at ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']); set(gca,'FontSize',11);
    pbaspect([1 1 1]);
end


k_sim = 1;
tmp_dynm = [DFC.w120(:,k_sim).^2,DFC.w240(:,k_sim).^2,...
    pMIseq.w120(:,k_sim),pMIseq.w240(:,k_sim),...
    squeeze(pMIwtseq.MI(Pmed+1:end,k_sim,:)),...
    permute(WTC.amor_pw(:,:,k_sim),[2,1,3])];

corr_dynm = zeros(size(tmp_dynm,2),size(tmp_dynm,2),5565);
tic; %
for k_sim = 1:5565
    if rem(k_sim,100) == 1
        disp(['corr *: ',num2str(round(k_sim/5565 *100)),'%'])
    end
    tmp_dynm = [DFC.w120(:,k_sim).^2,DFC.w240(:,k_sim).^2,... % 1:2
        pMIseq.w120(:,k_sim),pMIseq.w240(:,k_sim),... % 3:4
        squeeze(pMIwtseq.MI(Pmed+1:end,k_sim,:)),... % 5:11
        permute(WTC.amor(:,:,k_sim),[2,1,3])]; % 13:86
    corr_dynm(:,:,k_sim) = corr(tmp_dynm,'Type','Spearman','Rows','pairwise');
end


corr_dynm_pw = zeros(size(tmp_dynm,2),size(tmp_dynm,2),5565);
tic; %
for k_sim = 1:5565
    if rem(k_sim,100) == 1
        disp(['corr *: ',num2str(round(k_sim/5565 *100)),'%'])
    end
    tmp_dynm = [DFC_pw.w120(:,k_sim).^2,DFC.w240(:,k_sim).^2,... % 1:2
        pMIseq_pw.w120(:,k_sim),pMIseq.w240(:,k_sim),... % 3:4
        squeeze(pMIwtseq.MI(Pmed+1:end,k_sim,:)),... % 5:11
        permute(WTC.amor_pw(:,:,k_sim),[2,1,3])]; % 13:86
    corr_dynm_pw(:,:,k_sim) = corr(tmp_dynm,'Type','Spearman','Rows','pairwise');
end
toc;




figure(18) % Fig 1
subplot(1,3,1);
plot(mean(WTC.amor_pw_meanabs(52:end,:),1),FC_pw.^2,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
xlim([0 1]); xticks(0:0.25:1); xlabel('tsWTC (a.u.)');
ylim([0 1]); yticks(0:0.25:1); ylabel('sqCC \itr\rm^2 (a.u.)');
set(gca,'FontSize',11)
pbaspect([1 1 1])
subplot(1,3,2);
plot(mean(WTC.amor_pw_meanabs(52:end,:),1),pMIseq_pw.MI,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
xlim([0 1]); xticks(0:0.25:1); xlabel('tsWTC (a.u.)');
ylim([0 lim_MI]); yticks(0:0.25:lim_MI);ylabel('pMI (bit)');
set(gca,'FontSize',11)
pbaspect([1 1 1])
subplot(1,3,3);
plot(FC_pw.^2,pMIseq_pw.MI,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
xlim([0 1]); xticks(0:0.25:1); xlabel('sqCC \itr\rm^2 (a.u.)');
ylim([0 lim_MI]); yticks(0:0.25:lim_MI);ylabel('pMI (bit)');
set(gca,'FontSize',11)
pbaspect([1 1 1])
svfig(18,800,260,'scatterplots',['comp_tsWTC03to08xFCxpMI_PW_',svfilehead],['SIM_compMethods_03_ATR4Hz',type{1,ttype}],1)
disp(['Spearman CC: [tsWTC x sqCC] ',num2str(corr(mean(WTC.amor_pw_meanabs(52:end,:),1)',FC_pw'.^2,'type','Spearman'))])
disp(['Spearman CC: [tsWTC x pMI] ',num2str(corr(mean(WTC.amor_pw_meanabs(52:end,:),1)',pMIseq_pw.MI','type','Spearman'))])
disp(['Spearman CC: [sqCC x pMI] ',num2str(corr(FC_pw'.^2,pMIseq_pw.MI','type','Spearman'))])

figure(19) % Fig S1
subplot(2,1,1);
plot(FC_pw,mean(WTC.amor_pw_meanabs(52:end,:),1),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
xlim([-1 1]); xticks(-1:0.25:1); xlabel('CC \itr \rm(a.u.)');
ylim([0 1]); yticks(0:0.25:1); ylabel('tsWTC (a.u.)');
set(gca,'FontSize',11)
pbaspect([2 1 1])
subplot(2,1,2);
plot(FC_pw,pMIseq_pw.MI,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
xlim([-1 1]); xticks(-1:0.25:1); xlabel('CC \itr \rm(a.u.)');
ylim([0 lim_MI]); yticks(0:0.25:lim_MI);ylabel('pMI (bit)');
set(gca,'FontSize',11)
pbaspect([2 1 1])


figure(26) % Fig S2
subplot(2,4,1);
plot(mean(WTC.amor_meanabs(41:end,:),1),mean(WTC.amor_pw_meanabs(52:end,:),1),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
hold on; plot([0 1],[0 1],'k--');
xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
xlabel('w/o PW'); ylabel('w/ PW');
title('tsWTC (0.08 - 0.03 Hz)')
% h1 = lsline; h1.Color = 'r'; h1.LineStyle = '-';
hold off;
set(gca,'FontSize',11)
pbaspect([1 1 1])
subplot(2,4,2);
plot(FC,FC_pw,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
hold on; plot([-1 1],[-1 1],'k--');
xlim([-1 1]); xticks(-1:0.5:1); xlabel('w/o PW');
ylim([-1 1]); yticks(-1:0.5:1); ylabel('w/ PW');
title('CC')
% h1 = lsline; h1.Color = 'r'; h1.LineStyle = '-';
hold off;
set(gca,'FontSize',11)
pbaspect([1 1 1])
subplot(2,4,3);
plot(FC.^2,FC_pw.^2,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
hold on; plot([-1 1],[-1 1],'k--');
xlim([0 1]); xticks(0:0.25:1); xlabel('w/o PW');
ylim([0 1]); yticks(0:0.25:1); ylabel('w/ PW');
title('sqCC')
% h1 = lsline; h1.Color = 'r'; h1.LineStyle = '-';
hold off;
set(gca,'FontSize',11)
pbaspect([1 1 1])
subplot(2,4,4);
lim_MI = ceil(max([pMIseq.MI,pMIseq_pw.MI])*10)/10;
plot(pMIseq.MI,pMIseq_pw.MI,'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
hold on; plot([0 lim_MI],[0 lim_MI],'k--');
xlim([0 lim_MI]); xticks(0:0.25:lim_MI); xlabel('w/o PW');
ylim([0 lim_MI]); yticks(0:0.25:lim_MI); ylabel('w/ PW');
title('pMI')
hold off;
set(gca,'FontSize',11)
pbaspect([1 1 1])
freqsets = [17,41,49,57];
for k_f = 1:size(freqsets,2)
    subplot(2,4,k_f+4);
    plot(WTC.amor_meanabs(freqsets(k_f),:),WTC.amor_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel('w/o PW'); ylabel('w/ PW');
    title(['tWTC ',num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz'])
    set(gca,'FontSize',11);
    pbaspect([1 1 1]);
end


figure(24) % Fig 2
plot(corr(FC_pw'.^2,(WTC.amor_pw_meanabs'),'Type','Spearman'),'+','MarkerSize',5,'MarkerEdgeColor',0.3*ones(1,3)); hold on;
plot(corr(pMIseq_pw.MI'.^2,(WTC.amor_pw_meanabs'),'Type','Spearman'),'^','MarkerSize',5,'MarkerEdgeColor',0.3*ones(1,3)); hold off;
xlim([1 inf]); xticks(5:5:55); xticklabels(round(WTC.FOI_meanabs(5:5:55),3)); xlabel('Frequency band (Hz)');
ylabel('Correlation (\rho)'); ylim([-0.1 1]);
legend({'sqCC','pMI'},'Location','best'); set(gca,'FontSize',11);


figure(25) % Fig S3
plot(corr(FC'.^2,(WTC.amor_meanabs'),'Type','Spearman'),'+','MarkerSize',5,'MarkerEdgeColor',0.3*ones(1,3)); hold on;
plot(corr(pMIseq.MI'.^2,(WTC.amor_meanabs'),'Type','Spearman'),'^','MarkerSize',5,'MarkerEdgeColor',0.3*ones(1,3)); hold off;
xlim([1 inf]); xticks(5:5:55); xticklabels(round(WTC.FOI_meanabs(5:5:55),3)); xlabel('Frequency band (Hz)');
ylabel('Correlation (\rho)'); ylim([-0.1 1]);
legend({'sqCC','pMI'},'Location','best'); set(gca,'FontSize',11);



figure(15) % Fig 3
options.color_area = [0.5,0.5,0.9];
options.color_line = [0.2,0.2,0.9];
options.alpha      = 0.2;
options.line_width = 2;
options.error      = 'std';
shadedErrorBar_auto(atanh(permute(corr_dynm_pw(1,13:86,:),[3,2,1])),options); hold on;
options.color_area = [0.9,0.5,0.5];
options.color_line = [0.9,0.2,0.2];
shadedErrorBar_auto(atanh(permute(corr_dynm_pw(2,13:86,:),[3,2,1])),options);  hold on;
options.color_area = [0.5,0.9,0.9];
options.color_line = [0.2,0.7,0.9];
shadedErrorBar_auto(atanh(permute(corr_dynm_pw(3,13:86,:),[3,2,1])),options);  hold on;
options.color_area = [0.9,0.9,0.5];
options.color_line = [0.9,0.7,0.2];
shadedErrorBar_auto(atanh(permute(corr_dynm_pw(4,13:86,:),[3,2,1])),options);   hold off;
xticks(38:5:75); xticklabels(round(WTC.Freq(38:5:75),3));
xlim([35 70]); xlabel('Frequency band (Hz)'); ylabel('Fisher Z (a.u.)');
legend({'','sqswCC 30','','sqswCC 60','','swpMI 30','','swpMI 60'},'Location','best')
set(gca,'FontSize',11)


figure(17) % Fig 4
tmp_corr_dynm = [atanh(squeeze(corr_dynm_pw(1,3,:))),atanh(squeeze(corr_dynm_pw(2,4,:)))];
pd1 = fitdist(tmp_corr_dynm(:,1),'Normal');
ci95_1 = pd1.paramci;
pd2 = fitdist(tmp_corr_dynm(:,2),'Normal');
ci95_2 = pd2.paramci;
barwitherr(std(tmp_corr_dynm,[],1),mean(tmp_corr_dynm,1),'facecolor',0.8*ones(1,3));
% barwitherr([diff(ci95_1(:,1))/2,diff(ci95_2(:,1))/2],mean(tmp_corr_dynm));
ylabel('Fisher Z (a.u.)'); ylim([-0.05 0.8]); xticklabels({'30 s','60 s'}); xlim([0.2 2.8])
set(gca,'FontSize',11);

corr_mWTCsqCCpMI_pw = corr([mean(WTC.amor_pw_meanabs(39:end,:),1);FC_pw.^2;pMIseq_pw.MI]','Type','Spearman');
disp({'mmWTC x sqCC','mmWTC x pMI','sqCC x pMI';corr_mWTCsqCCpMI_pw(2,1),corr_mWTCsqCCpMI_pw(3,1),corr_mWTCsqCCpMI_pw(3,2)})



figure(21) % Fig S4
freqsets = [17,41,49,57];
for k_f = 1:size(freqsets,2)
    subplot(size(freqsets,2),3,1 + 3*(k_f-1));
    plot(WTC.amor_pw_meanabs(freqsets(k_f),:),WTC.cgau1_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel('amor');
    ylabel('cgau1');
    title([num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    set(gca,'FontSize',11);
    pbaspect([1 1 1]);
    subplot(size(freqsets,2),3,2 + 3*(k_f-1));
    plot(WTC.amor_pw_meanabs(freqsets(k_f),:),WTC.cgau2_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel('amor');
    ylabel('cgau2');
    title([num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    set(gca,'FontSize',11);
    pbaspect([1 1 1]);
    subplot(size(freqsets,2),3,3*k_f);
    plot(WTC.cgau1_pw_meanabs(freqsets(k_f),:),WTC.cgau2_pw_meanabs(freqsets(k_f),:),'.','MarkerSize',2,'MarkerEdgeColor',0.3*ones(1,3));
    hold on; plot([0 1],[0 1],'k--'); hold off;
    xlim([0 1]); ylim([0 1]); xticks(0:0.25:1); yticks(0:0.25:1);
    xlabel('cgau1');
    ylabel('cgau2');
    title([num2str(round(WTC.FOI_meanabs(freqsets(k_f)),3)),' Hz']);
    set(gca,'FontSize',11);
    pbaspect([1 1 1]);
end



figure(515) % Fig S5
options.color_area = [0.5,0.5,0.9];
options.color_line = [0.2,0.2,0.9];
options.alpha      = 0.2;
options.line_width = 2;
options.error      = 'std';
shadedErrorBar_auto(atanh(permute(corr_dynm(1,13:86,:),[3,2,1])),options); hold on;
options.color_area = [0.9,0.5,0.5];
options.color_line = [0.9,0.2,0.2];
shadedErrorBar_auto(atanh(permute(corr_dynm(2,13:86,:),[3,2,1])),options);  hold on;
options.color_area = [0.5,0.9,0.9];
options.color_line = [0.2,0.7,0.9];
shadedErrorBar_auto(atanh(permute(corr_dynm(3,13:86,:),[3,2,1])),options);  hold on;
options.color_area = [0.9,0.9,0.5];
options.color_line = [0.9,0.7,0.2];
shadedErrorBar_auto(atanh(permute(corr_dynm(4,13:86,:),[3,2,1])),options);   hold off;
xticks(38:5:75); xticklabels(round(WTC.Freq(38:5:75),3));
xlim([35 70]); xlabel('Frequency band (Hz)'); ylabel('Fisher Z (a.u.)');
legend({'','sqswCC 30','','sqswCC 60','','swpMI 30','','swpMI 60'},'Location','best')
set(gca,'FontSize',11)


figure(521) % Fig S6
tmp_corr_dynm = [atanh(squeeze(corr_dynm(1,3,:))),atanh(squeeze(corr_dynm(2,4,:)))];
pd1 = fitdist(tmp_corr_dynm(:,1),'Normal');
ci95_1 = pd1.paramci;
pd2 = fitdist(tmp_corr_dynm(:,2),'Normal');
ci95_2 = pd2.paramci;
barwitherr(std(tmp_corr_dynm,[],1),mean(tmp_corr_dynm,1),'facecolor',0.8*ones(1,3));
% barwitherr([diff(ci95_1(:,1))/2,diff(ci95_2(:,1))/2],mean(tmp_corr_dynm));
ylabel('Fisher Z (a.u.)'); ylim([-0.05 0.8]); xticklabels({'30 s','60 s'}); xlim([0.2 2.8])
set(gca,'FontSize',11);

