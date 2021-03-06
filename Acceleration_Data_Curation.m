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
%% important parameters for further computation 
Ts      = mean(diff(TimeVect)); % this the period 
Fs      = 1/Ts;                 % the sampling frequency
%% extracting the x, y and z coordinates 
g_x     = ACC_RES(:,1);
g_y     = ACC_RES(:,2);
g_z     = ACC_RES(:,3);
%% calculate the g-level norm
g = sqrt(g_x.^2 + g_y.^2 + g_z.^2);
%% performin a Welch power spectrum analysis of the g-norm signal 
[accPow,f]=pwelch(g,[],[],[1:0.1:10],Fs); % because I did not downsample the signal so I am limiting the frequency range 
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
