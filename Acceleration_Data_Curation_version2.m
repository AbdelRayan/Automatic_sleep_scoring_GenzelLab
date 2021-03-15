% This is the second version of acceleration data curation --> more
% organized features and work flow 

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
%% Now we need to downsample the data to reduce the noise and flickering of the baseline 
% the data is downsampled to 100 Hz --> could be reduced more 
NDown = 200; % how many times you want to reduce the sampling of the data 
FsDown = 20000/NDown;
g_x = decimate(ACC_RES(:,1),NDown,'FIR');
g_y = decimate(ACC_RES(:,2),NDown,'FIR');
g_z = decimate(ACC_RES(:,3),NDown,'FIR');
%% calculate the g-level norm
g = sqrt(g_x.^2 + g_y.^2 + g_z.^2);
% create time vector 
TimeVectDown = linspace(0,max(TimeVect),numel(g));
%% the next step is to remove the drifts from the data to make the analysis easier 
g_detrend = locdetrend(g,FsDown,[.1 .01]); 
envelopes = abs(hilbert(g_detrend)');
%% creating a data matrix for HMM 
dataMat = [g_detrend envelopes'];
%% generating a threshold for the data 
% Here we generate an amplitude threshold using the Shin 2018 method
threshold = 2*mean(envelopes,2);
%%
figure 
plot(TimeVectDown,g,'LineWidth',2)
hold on 
plot(TimeVectDown,g_detrend,'r','LineWidth',2)
% plot(TimeVectDown,envelopes,'k','LineWidth',2)
yline(threshold,'b--','LineWidth',2)
yline(-threshold,'b--','LineWidth',2)
% legend('Raw Signal','Detrended Signal','Envelope')
% legend('Raw Signal','Detrended Signal','Location','Best')
xlabel('Time [s]')
ylabel('g-level value')
box off 
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
xlim([480 520])
axis off
export_fig('DetrendedAccelerometerDataCloseUp','-jpg','-r300','-q70','-transparent')
%% running HMM model --> determining the parameters for the model 

options = struct();
options.initrep = 2; % Set to be large as this is a short simulation
options.K = 2;
options.standardise = 0;
options.verbose = 1;
options.Fs = 100;
options.useMEX = 1;
options.zeromean = 0; % changed it from zero 
options.dropstates = 1;
options.order = 0;
% options.timelag = 100; % added a time lag in the form of how many samples 
options.DirichletDiag = 1000; % set very large as we don't have much data
%%
% HMM inference - we only store the Gamma time-course of posterior probabilities
T = length(g_detrend);
[hmm, Gamma_emb{1},~,vpath] = hmmmar(g_detrend,T,options);
disp('Finished')
%% plotting the results of the classification 
figure
norm_g = (g_detrend - min(g_detrend)) / ( max(g_detrend) - min(g_detrend) );
hold on
area(TimeVectDown,Gamma_emb{1})
plot(TimeVectDown,norm_g,'k','LineWidth',2)
yline(threshold,'b--','LineWidth',2)
yline(-threshold,'b--','LineWidth',2)
% plot(TimeVectDown,vpath,'b.')
% xlim([200 1000])
%% comparing with real data 
state1 = TimeVectDown(vpath ==1);
state2 = TimeVectDown(vpath == 2);
NewStates = zeros(1,numel(states));
NewStates(states ==1 ) = 2;
NewStates(states > 1 ) = 1;
%% 
vpathDown = downsample(vpath,100);
%% 
timeVectStates = 1:1:numel(states);
plot(timeVectStates,NewStates,'r.')
% hold on 
% plot(TimeVectDown,vpath,'k.')
ylim([0 3])
xlim([200 1000])
%% making plots for the meeting on Wednesday 
g_x = ACC_RES(:,1);
g_y = ACC_RES(:,2);
g_z = ACC_RES(:,3);
g = sqrt(g_x.^2 + g_y.^2 + g_z.^2);
%% plotting the raw data 
figure;
plot(TimeVect,g,'LineWidth',2)
xlabel('Time [s]')
ylabel('g-level [a.u.]')
box off 
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
export_fig('RawDataNotDownSampled','-jpg','-r300','-q70','-transparent')
