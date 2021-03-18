% The aim of this script is to compute the theta delta ration from the
% ephys data and correlate it with the accelerometer data 

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
% the data is downsampled to 100 Hz --> could be reduced more 
NDown = 200; % how many times you want to reduce the sampling of the data 
FsDown = 20000/NDown;
g_x = decimate(ACC_RES(:,1),NDown,'FIR');
g_y = decimate(ACC_RES(:,2),NDown,'FIR');
g_z = decimate(ACC_RES(:,3),NDown,'FIR');
%% important parameters for further computation 
Ts      = mean(diff(TimeVect)); % this the period 
Fs      = 1/Ts;                 % the sampling frequency
%% calculate the g-level norm
g = sqrt(g_x.^2 + g_y.^2 + g_z.^2);
TimeVectDown = linspace(0,max(TimeVect),numel(g));
%% the next step is to remove the drifts from the data to make the analysis easier 
g_detrend = locdetrend(g,FsDown,[.1 .01]); 
%%
PowAccel = powerEMG( g_detrend', 100,  4, 1 );
%% creating a histogram for the detrended data 
[h,n] = hist(g_detrend,-0.03:0.001:0.03);
figure
bar(n,h/sum(h))
xlabel('Accelerometer')
ylabel('Probability')
box off 
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
% export_fig('HistogramAccelerometerData','-pdf','-r300','-q70','-transparent')
%% uploading the ephys data
%In each case, there is at least one file to load
disp("Loading and downsampling of the data for the channels of interest")
[DataHPC, ~, ~] = load_open_ephys_data_faster('100_CH2.continuous');
[DataPFC, ~, ~] = load_open_ephys_data_faster('100_CH33.continuous');
%% downsampling the ephys data 
DataHPCDown = decimate(DataHPC,16,'FIR');
DataPFCDown = decimate(DataPFC,16,'FIR');
% locally detrending the data 
DataPFCDownD = locdetrend(DataPFCDown,1250,[.1 .05]); 
% creating a downsampled time vector 
TimeVectDownLFP = linspace(0,max(TimeVect),numel(DataPFCDown));
%% plotting the ephys data 
figure 
plot(TimeVectDownLFP,DataPFCDown)
%% computing the delta/theta ration in windows of 4 seconds 
% divide the data into epochs 
fs = 1250;
epochLenWant = 4;
epochLendp = epochLenWant*fs;
EpochNum = floor(length(DataPFCDown)/epochLendp);
dataEpochs = zeros(epochLendp ,EpochNum); % preloc

startI = 1;
EndI = epochLendp;

for iEpoch = 1:EpochNum
    dataEpochs(:,iEpoch) = DataPFCDown(startI:EndI);
    startI = startI +epochLendp;
    EndI = EndI + epochLendp;
end
%% verification step --> plotting one epoch 
figure
plot(dataEpochs(:,200))
%% computing the spectra and the ration 
winsize                  = 1250;                      % Equivalent to 0.25 s
overlap                 = ceil(winsize/2);                      % Equivalent to 50% of the window
nfft                    = 2^nextpow2(winsize);
deltaIdx = 3:9;
thetaIdx = 10:19;
% thetaIdx = 121:201;
for iEpochS = 1:size(dataEpochs,2)
    [~,F,~,P] = spectrogram(dataEpochs(:,iEpochS),hanning(winsize),overlap,nfft,1250);
    meanPowr(:,iEpochS) = mean(P,2);
    meandelta = mean(meanPowr(deltaIdx,iEpochS));
    meantheta = mean(meanPowr(thetaIdx,iEpochS));
    deltaTheta(iEpochS,1) = meandelta / meantheta;
end
%% the average accelerometer data in 4 secods 
fsA = 100;
epochLenWant = 4;
epochLendpA = epochLenWant*fsA;
EpochNumA = floor(length(g_detrend)/epochLendpA);
AcceleroEpochs = zeros(epochLendpA ,EpochNumA); % preloc

startIA = 1;
EndIA = epochLendpA;

for iEpoch = 1:EpochNumA
    AcceleroEpochs(:,iEpoch) = g_detrend(startIA:EndIA);
    startIA = startIA +epochLendpA;
    EndIA = EndIA + epochLendpA;
end
% computing the mean per Epoch
meanAcceleroEpoch = mean(AcceleroEpochs,1);
%% interpolation 
spk_x       =interp1(TimeVectDown,vpath,...
            MatCell.SpikeMatrix.TimeStampSpike);  
%% making a scatter 
figure 
scatter(log(deltaTheta),log(PowAccel),'.')
xlabel('delta/theta')
ylabel('Accelerometer')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
% export_fig('ScatterACCeleroDeltaTheta','-pdf','-r300','-q70','-transparent')
% ylim([-0.03 0.03])
%% 
histogram(log(deltaTheta),50)

%% 
histogram(log(meanAcceleroEpoch+0.1),50)