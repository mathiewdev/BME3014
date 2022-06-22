% Mathiew Daniels-Diehl

% importing the data file
raw = importdata('hw1_data.dat');

% lowpass filter design
lpfilt = designfilt('lowpassfir','PassbandFrequency',250,...
 'StopbandFrequency', 325,'StopbandAttenuation',65,'SampleRate',1000);

% assigning x values (1 second for every 1000 data entries since sample rate is 1000 Hz) and plotting the raw data
time1 = linspace(1,2, numel(y));

plot(time1, raw)
title('Raw Data Over Time')
xlabel('Time (seconds)')
ylabel('mV')

% applying low-pass filter to the raw data
dataFiltered = filter(lpfilt, y);

% how many samples the filter delays the signal
gd = grpdelay(lpfilt);
% gd = 15 = 0.015s
delay1 = 0.015;
delay2 = 15;

%delaying the data and time signature, and plotting the new delayed data
timeDelay1 = delayseq(time1, delay1);
signalDelay1 = delayseq(dataFiltered, delay2);
plot (timeDelay1, signalDelay1)
title('Filtered Data Over Time')
xlabel('Time (seconds)')
ylabel('mV')

%Designing bandpass FIR filter
bpFilt = designfilt('bandpassfir', 'PassbandFrequency1', 125,...
    'StopbandFrequency1', 50, 'PassbandFrequency2', 250, 'StopbandFrequency2', 300,...
    'SampleRate', 1000);

% filtering raw data with new bp filter
dataFiltered2 = filter(bpFilt, y);

% Delay of bp filter
gd2 = grpdelay(bpFilt);
delay3 = 0.023;
delay4 = 23;

% Accounting for delay
timeDelay2 = delayseq(time1, delay3);
signalDelay2 = delayseq(dataFiltered2, delay4);
plot (x3, y3)
title('Band-Pass Filtered Data Over Time')
xlabel('Time (seconds)')
ylabel('mV')


% Mean-squared error for noise and actual signal
actual = importdata('signal.mat');
err = raw - actual;
mse = (sum(err, 'all'))/numel(err);

% MSE for low pass and actual
err2 = signalDelay1 - actual;
mse2 = (sum(err2, 'all'))/numel(err2);

% MSE bandpass filter
err3 = signalDelay2 - actual;
mse3 = (sum(err3, 'all'))/numel(err3);





