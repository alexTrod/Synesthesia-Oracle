clear all; close all;clc

% define pathways
masterdir = '/research/cisc2/projects/ward_connectomesyn';

% load data
load('MSMall_data_all.mat')

datasetlist = {'MSMall_data_hcpmain', 'MSMall_data_hcpDA', 'MSMall_data_mcpf5076'};
Pnetmats_icov_prep_all = [];

for idataset = 1:numel(datasetlist)

    MSMall_data = eval(datasetlist{idataset}); 

    %% preprocess, calculate full and partial FC matrices, plot

%     MSMall_data=MSMall_data(:,1:360,:); %take just cortical
    % data wrangle for FSLnets
    % concatenate all subjects' timepoints per node

    % Bandpass filter settings
    TR = 0.8;                                       % Temporal resolution
    fnq = 1/(2*TR);                                 % Nyquist frequency
    flp = 1/(60*1);                               % lowpass frequency of filter (Hz)
    fhi = 0.15;                                   % highpass
    Wn = [flp/fnq fhi/fnq];                         % butterworth bandpass non-dimensional frequency
    k = 2;                                          % 2nd order butterworth filter
    [bfilt,afilt] = butter(k,Wn);


    MSMall_data_concat_raw = [];
    MSMall_data_concat_prep = [];

    for ii = 1:size(MSMall_data,1)
        tmp_raw = squeeze(MSMall_data(ii,:,:))';
        tmp = detrend(tmp_raw);
        tmp = filtfilt(bfilt,afilt,tmp); % Bandpass filtering
        MSMall_data_concat_raw = [MSMall_data_concat_raw; tmp_raw];
        MSMall_data_concat_prep = [MSMall_data_concat_prep; tmp];

    end

    %make struct for fslnets
    ts_raw = [];
    ts_raw.Nsubjects = size(MSMall_data,1);
    ts_raw.Nnodes = size(MSMall_data,2);
    ts_raw.Ntimepoints = size(MSMall_data,3)*size(MSMall_data,1);
    ts_raw.NtimepointsPerSubject = ts_raw.Ntimepoints/ts_raw.Nsubjects;
    ts_raw.DD = 1:ts_raw.Nnodes;
    ts_raw.tr = 0.8;
    ts_raw.ts = MSMall_data_concat_raw;  

    ts_prep = [];
    ts_prep.Nsubjects = size(MSMall_data,1);
    ts_prep.Nnodes = size(MSMall_data,2);
    ts_prep.Ntimepoints = size(MSMall_data,3)*size(MSMall_data,1);
    ts_prep.NtimepointsPerSubject = ts_prep.Ntimepoints/ts_prep.Nsubjects;
    ts_prep.DD = 1:ts_prep.Nnodes;
    ts_prep.tr = 0.8;
    ts_prep.ts = MSMall_data_concat_prep;  

    ts_spectra_raw = nets_spectra(ts_raw);
    ts_spectra_prep = nets_spectra(ts_prep);

    Fnetmats_raw = nets_netmats(ts_raw,1,'corr'); % full correlation with r-to-z
    Fnetmats_prep = nets_netmats(ts_prep,1,'corr'); % full correlation with r-to-z

    Pnetmats_icov_raw = nets_netmats(ts_raw,1,'icov'); % partial correlation with r-to-z, inverse covariance
    Pnetmats_icov_prep = nets_netmats(ts_prep,1,'icov'); % partial correlation with r-to-z, inverse covariance

    % group-level analysis to look at the mean network matrix across all subjects:
    % simple average of netmats accross all subjects (Mnet)
    % simple one-group t-test (against zero) across subjects as Z values (Znet)
    [Znet_F_raw,Mnet_F_raw] = nets_groupmean(Fnetmats_raw,1);
    [Znet_F_prep,Mnet_F_prep] = nets_groupmean(Fnetmats_prep,1);

    [Znet_P_icov_raw,Mnet_P_icov_raw] = nets_groupmean(Pnetmats_icov_raw,1);
    [Znet_P_icov_prep,Mnet_P_icov_prep] = nets_groupmean(Pnetmats_icov_prep,1);
    
    
    Pnetmats_icov_prep_all = [Pnetmats_icov_prep_all; Pnetmats_icov_prep];    
    
    
end

save('parcor_all.mat', 'Pnetmats_icov_prep_all');



% % plot correlation matrix itself
% figure('position',[100 100 1100 400]);
% subplot(1,2,1); 
% imagesc(Mnet_F_raw);  colormap('jet');  colorbar; axis square;
% title('Mean Pearson correlation - "raw" data');
% subplot(1,2,2); 
% imagesc(Mnet_F_prep);  colormap('jet');  colorbar; axis square;
% title('Mean Pearson correlation - preprocessed data');
% 
% figure('position',[100 100 1100 400]);
% subplot(1,2,1); 
% imagesc(Mnet_P_icov_raw);  colormap('jet');  colorbar; axis square;
% title('Mean Pearson correlation - "raw" data');
% subplot(1,2,2); 
% imagesc(Mnet_P_icov_prep);  colormap('jet');  colorbar; axis square;
% title('Mean Pearson correlation - preprocessed data');


