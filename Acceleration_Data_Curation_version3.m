%% The aim of this script is to try the approach reported in Santos Lima et al
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
%% calculate the g-level using the new method
g_x_mean = mean(g_x);
g_y_mean = mean(g_y);
g_z_mean = mean(g_z);

g_sum = (g_x - g_x_mean) + (g_y - g_y_mean) + (g_z - g_z_mean);
%% locally detrending the data and creating a time vector 
g_detrend = locdetrend(g_sum,FsDown,[.1 .01]); 
TimeVectDown = linspace(0,max(TimeVect),numel(g_sum));
threshold = 2*mean(g_sum);
%% making a plot to compare and confirm the results 
figure 
plot(TimeVectDown,g_sum,'LineWidth',2)
hold on 
plot(TimeVectDown,g_detrend,'r','LineWidth',2)
% plot(TimeVectDown,envelopes,'k','LineWidth',2)
% yline(threshold,'b--','LineWidth',2)
% yline(-threshold,'b--','LineWidth',2)
% legend('Raw Signal','Detrended Signal','Envelope')
% legend('Raw Signal','Detrended Signal','Location','Best')
xlabel('Time [s]')
ylabel('g-level value')
box off 
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
% xlim([480 520])
% axis off
% export_fig('DetrendedAccelerometerDataCloseUp','-jpg','-r300','-q70','-transparent')
%% running HMM model --> determining the parameters for the model 

options = struct();
options.initrep = 5; % Set to be large as this is a short simulation
options.K = 2;
options.standardise = 0;
options.verbose = 1;
options.Fs = 100;
options.useMEX = 1;
options.zeromean = 0; % changed it from zero 
options.dropstates = 1;
options.order = 0;
options.timelag = 10; % added a time lag in the form of how many samples 
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
a = area(TimeVectDown,Gamma_emb{1})
plot(TimeVectDown,norm_g,'k','LineWidth',2)
xlabel('Time [s]')
ylabel('Normalized g-level')
% axis off 
box off 
a(1).FaceColor = 'r';
a(1).FaceAlpha = 0.4;
a(2).FaceColor = 'b';
a(2).FaceAlpha = 0.4;
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
set(gcf,'Color','w')
% xlim([950 1050])
% export_fig('AccelerometerSleepWakeHMMCloseUp','-pdf','-r300','-q70','-transparent')
% print(gcf, '-dpdf', 'AccelerometerSleepWakeHMMCloseUp.pdf');
% yline(threshold,'b--','LineWidth',2)
% yline(-threshold,'b--','LineWidth',2)
% plot(TimeVectDown,vpath,'b.')
% xlim([200 1000])