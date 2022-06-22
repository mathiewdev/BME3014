
%% GROUP INFORMATION
%--> In these comments input group number and member names
% Group number: 5
% Group members: Mathiew Daniels-Diehl, Naisargi Mehta, Jack Rothenberg

%% DEFINE INPUTS
% (1) Set the filename here to match your baseline file
% (2) Set the 'Fs' variable to match your sampling frequency in Hz
fnameRest = "BME3014_Lab1Data_group5_Exercising.txt";
Fs = 200;

%% Load data into matlab
% (3) Use the 'importdata' or 'dlmread' matlab function to read the text from your data
%   file. For help type 'doc importdata' in the command window
% (4) Extract the raw data from the 3rd column of imported data and put it
%   in the 'rawdata' variable
% (5) Extract the time from the 1st column of imported data and put it
%   in the 'rawdata' variable


rawdata = importdata(fnameRest).data(:,3);
time =  importdata(fnameRest).data(:,1);

%datafile = importdata(fnameRest);
%rawdata = datafile.data(:,3);
%time = datafile.data(:,1);

%% Initialize plots
close all
%Plot raw data
   samples = 1:3000;% Number of samples to graph

    % Generate plots
    n = 6;
    figure
    subplot(n,1,1)
    plot(time(samples),rawdata(samples),'b-')
    title('ECG signal')
    xlabel('Time (s)')
    
%% BANDPASS FILTER
% APPLY LOWPASS FILTER - Get Coefficients from Pan and Tompkins
% (6) Using the Pan and Tomkins paper, extract the B coefficients from
%   the transfer function (Equation 1). Define them using the same method
%   shown for the A coefficients.
B = [1 zeros(1,5) -2 zeros(1,5) 1];
A = [1 -2 1];

lpdata = filter(B,A,rawdata); % Apply the lowpass filter to the data

% Account for the data delay caused by the filter (check Pan and Tompkins)
delay = 5;
lpdata = lpdata(delay:end);

% APPLY HIGHPASS FILTER -- Use the example above to apply Pan and Tompkins'
%   high pass filter to the data: FIND THE CORRECT EQUATION IN THE ERRATA. 
%   Save the filtered data as 'hpdata'
% (7) Define the high-pass filter coefficients for A and B by extracting
%   them from eqn 4 from pan and tompkins
% (8) Use the matlab 'filter' function as above to apply the high pass
%   filter to the low pass filtered data
% (9) Read through the Pan and Tompkins paper and determine the delay caused
%   by the filter (Explicitely stated in the paper)--account for this delay
%   as was done for the low pass filter

B = [-1/32 zeros(1,15) 1 -1 zeros(1,14) 1/32];

A = [1 -1];

hpdata = filter(B,A,lpdata);%apply the high pass filter to the data

delay = 16;
hpdata = hpdata(delay:end);% apply the filter delay

% Plot bandpass filtered data
    subplot(n,1,2)
    plot(time(samples), hpdata(samples),'b-')
    title('Band pass filtered signal')
    xlabel('Time (s)')


%% DERIVATIVE FILTER
% Apply the derivative filter defined by Pan and Tompkins
B = (1/8)*[-1 -2 0 2 1];
A = 1;

derdata = filter(B,A,hpdata);

% Account for the data delay caused by the filter (check Pan and Tompkins)
delay = 2;
derdata = derdata(delay:end);

% Square the derivative
% (10) Create a new variable (dersquared) that contains the derivative
%    squared
dersquared = derdata.^2;

% Plot data after derivative filter and squaring
subplot(n,1,3)
    plot(time(samples),derdata(samples),'b-')
    title('Output of derivative filter')
    xlabel('Time (s)')
    
    subplot(n,1,4)
    plot(time(samples),dersquared(samples),'b-')
    title('Output of derivative filter squared')
    xlabel('Time (s)')

%% Moving window integration
% A moving window integration is the equivalent of a moving average filter. 
% (11)Create filter coefficients for a moving average filter (30
%   sample rectangular window).
% (12) Use the filter function to apply the moving window integration filter 


B = 1/30*ones(1,30);
A = 1;

intdata = filter(B,A,dersquared);

intdata = intdata/max(intdata); % Normalize intdata

delay = 15;
intdata = intdata(delay:end);

% Plot data after integration filter
    subplot(n,1,5)
    plot(time(samples),intdata(samples),'b-')
    title('Output of integration  filter')
    xlabel('Time (s)')


    % Plot data after integration filter
   % P = time(samples);
    %Q = [0 ; diff(intdata(samples))];
    %subplot(n,1,6)
    %plot(P,Q,'b-')
    %title('Output of integration  filter')
    %xlabel('Time (s)')



%% Threshold data
% (13) Run the code up to this point on your baseline dataset. The above line
%   of code will produce a graph of the first few beats after
%   implementation of the moving window integration.
% (14) Using this graph, choose an appropriate threshold level to isolate the
%   QRS waves
level = 0.05;

threshdata = false(size(intdata)); % Create a thresholded array
threshdata(intdata > level) = true; % Create a thresholded array

%Get points where data changes
% (15) Run the code up to this point. The above plot command will produce a
%   plot of the first few beats of thresholded data
% (16) Find a way to identify the start and end of each rwave 
%   (HINT: type 'doc diff' into the matlab command window)
% (17) Create two variables (qrsstart and qrsend) which define the first
%   index and last index of each qrs wave

%Get points where data changes

diffdata = diff(threshdata);
qrsstart = find(diffdata == 1);
%len = length(qrsstart);
%qrsstart2 = len - 1;
qrsend = find(diffdata == -1);

%% Code to find mxima of each QRS wave
rpeak = zeros(size(qrsstart)-1);
amp = zeros(size(qrsstart)-1);
% find peaks and troughs of qrs waves
for i = 1:length(qrsstart)- 1
  [qrsmax,curind] = max(hpdata(qrsstart(i):qrsend(i)));
  rpeak(i) = curind+qrsstart(i)-1;
  [qrsmin,~] = min(hpdata(qrsstart(i):qrsend(i)));
  amp(i) = qrsmax-qrsmin;
end

% Determine rr-intervals
rrint = diff(rpeak)/Fs; % convert intervals from frames to seconds

% (18) Convert the rrint to instantaneous heart rate in beats per minute
HR = 60 ./ rrint;

%% Plot final results
     
    subplot(n,1,6)
    plot(time(samples),hpdata(samples),'b-')
    title('Selected ECG sections from threshold')
    xlabel('Time (s)')
    
    regiondata = hpdata(samples).*threshdata(samples);
    regiondata(regiondata==0)=NaN;
    hold on
    plot(time(samples),regiondata,'r-')
    figure
    
    plot(time(samples),hpdata(samples),'b-')
    title('Band pass filtered signal with fiducial marks')
    xlabel('Time (s)')
    ylabel('Filtered ECG signal (V)')
    hold on
    plot(time(rpeak(rpeak<max(samples))),hpdata(rpeak(rpeak<max(samples))),'ro')

%% IF YOU HAVE TIME AT THE END OF DAY 1
% -Determine a method to track heart rate changes over time
% -Refer to the papers for lab 1B (3/23) on myWPI and begin attempting to
%   extract respiratory rate data from the ECG signal
    
%% 1B Part 1 

fnameAct = "BME3014_Lab1Data_group5_Exercising.txt";

timeRest = importdata(fnameRest).data(:,1);
dataRest = importdata(fnameRest).data(:,4);
timeAct = importdata(fnameRest).data(:,1);
dataAct = importdata(fnameAct).data(:,4);

n = 7;
figure
subplot(n,1,1);
plot(timeRest, dataRest);
title('Rawdata Plot')
xlabel('Time (sec)')
ylabel('V')

lpFilt = designfilt('lowpassfir','PassbandFrequency',0.2,...
 'StopbandFrequency', 0.55,'StopbandAttenuation',60,'SampleRate',200);

dataRestlp = filter(lpFilt, dataRest);

lpdelay = grpdelay(lpFilt,1); 
dataRestlp = dataRestlp(lpdelay:end);
lptime = linspace(1,180,size(dataRestlp,1));

subplot(n,1,2);
plot(lptime, dataRestlp);
title('Lowpass Filtered Plot')
xlabel('Time (sec)')
ylabel('V')

subplot(n,1,3);
pwelch(dataRestlp,[],[],[],200)
xlim([0 5]);

%% Part 2 RSA Method

xq = (0: 1/200 : 180); 
x = time(rpeak(1:end-1));
y = rrint;
subplot(n,1,4);
hold on
plot(x,y,'b');
yinterp = interp1(x,y,xq,'pchip');
scatter(xq,yinterp,'r.');
title('RSA Method')
xlabel('Time (sec)')
ylabel('V')
hold off

subplot(n,1,5);
pwelch(yinterp,[],[],[],200)
xlim([0,15]);

%% Part 3 EDR Method 

xq = (0: 1/200 : 180); 
x = time(rpeak(1:end));
y = amp;
subplot(n,1,6);
hold on

plot(x,y,'b');
yinterp = interp1(x,y,xq,'pchip');
scatter(xq,yinterp,'r.');
title('EDR Method')
xlabel('Time (sec)')
ylabel('V')
hold off 

subplot(n,1,7)
pwelch(yinterp,[],[],[],200)
xlim([0,15]);