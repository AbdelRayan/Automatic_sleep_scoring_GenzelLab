% acceleration data 
% the aim of this script is to do experimentation on the accelerometer data
% and how to analyze it the best way 

clear all; close all; clc 
% 1. uploading the accelerometer files 
[file,path]=uigetfile('*.cont*','Select the accelerometer file(s) (select multiple if possible)', 'MultiSelect', 'on');
patIDacc = {};
for jj=1:length(file)
    tmp = file(jj);
    patIDacc{jj} = strcat([path,tmp{1}]);
end

% 2. reading the data from the openephys 
ACC_RES = [];
disp('Loading and downsampling of the accelerometer data')
for jj=1:length(patIDacc)
    [Dataacc, TimeVect, ~] = load_open_ephys_data_faster(patIDacc{jj});
%     Dataacc=filtfilt(b,a,Dataacc);
%     Dataacc=downsample(Dataacc,acq_fhz/500);
    ACC_RES = [ACC_RES  Dataacc]; 
end
disp('Finished reading the accelerometer data')
%% working with downsampled data to reduce noise 
ACC_RES1 = decimate(ACC_RES(:,1),200,'FIR');
ACC_RES2 = decimate(ACC_RES(:,2),200,'FIR');
ACC_RES3 = decimate(ACC_RES(:,3),200,'FIR');
%% important parameters for further computation 
Ts      = mean(diff(TimeVect)); % this the period 
Fs      = 1/Ts;                 % the sampling frequency
%% extracting the x, y and z coordinates 
g_x     = ACC_RES1;%(:,1);
g_y     = ACC_RES2;%(:,2);
g_z     = ACC_RES3;%(:,3);
%% calculate the g-level norm
g = sqrt(g_x.^2 + g_y.^2 + g_z.^2);
%% performin a Welch power spectrum analysis of the g-norm signal 
[accPow,f]=pwelch(g,[],[],[1:0.1:50],100); % because I did not downsample the signal so I am limiting the frequency range 
% from the power spectrum, I notice a fundmental frequency up 3 Hz
%% filtering the g-norm signal I am following the method from Carr et al in Nature 
% Zero-Delay Filter Options for G-level filtering
HPF = 0.01;     % Half Power Frequency (HPF) [Hz]
% Other Filter options
FilterOptions = {'lowpassiir','FilterOrder',12,'DesignMethod','butter'};
%% designing a filter 
% Design the filter
d1 = designfilt(FilterOptions{:},'HalfPowerFrequency',HPF);
% Do the filtering
g_x_filt = filtfilt(d1,g_x);
g_y_filt = filtfilt(d1,g_y);
g_z_filt = filtfilt(d1,g_z);
%
%% Calculate the filtered G-level (norm of the filtered g vector)
g_filt = sqrt(g_x_filt.^2 + g_y_filt.^2 + g_z_filt.^2);
%% I will be trying to use HMM to classify sleep and awake 

dataMat = [NewAcceleration envelopes'];
% The amplitude envelope is computed using the hilbert transform
envelopes = abs(hilbert(dataMat)');

% Here we generate an amplitude threshold using the Shin 2018 method
threshold = 2*median(envelopes,2);

%% HMM
 options = struct();
    options.initrep = 2; % Set to be large as this is a short simulation
    options.K = 2;
    options.standardise = 0;
    options.verbose = 1;
    options.Fs = 100;
    options.useMEX = 1;
    options.zeromean = 0;
    options.dropstates = 1;
    options.order = 0;
    options.DirichletDiag = 1000; % set very large as we don't have much data

    % HMM inference - we only store the Gamma time-course of posterior probabilities
    T = length(NewAcceleration);
    [hmm, Gamma_emb{1},~,vpath] = hmmmar(NewAcceleration,T,options);
    disp('Finished')
%% I would like here to remove the drifts from the data for better analysis
NewAcceleration = locdetrend(g_filt,20000,[.1 .05]); 
envelopes = abs(hilbert(NewAcceleration)');
%%

plot(NewAcceleration)
% histogram(envelopes,0:0.01:2)
% threshold = 2*median(zscore(NewAcceleration));
% plot(NewAcceleration(1:5000000))
hold on
yline(threshold(1),'r--','LineWidth',2)
yline(-threshold(1),'r--','LineWidth',2)
%%
area(1:1:T,Gamma_emb{1})
hold on
plot(1:1:T,NewAcceleration*10,'k','LineWidth',2)


xlim([1 90000])
%%
matDecode = Gamma_emb{1};
area(matDecode)