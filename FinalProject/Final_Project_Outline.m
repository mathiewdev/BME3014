clear all
close all
%% Import Data
% Please remember that this code must prompt the user for the file name.Use
% funvtion inputdlg(). Having done so, use the filename input to import
% your data set. Plot the data set to assure that it is importing the
% correct data. You are allowed to ask the user for if the data set came
% from a healthy or infracted heart.
% 
%This code imports the original files and removes the header
prompt=inputdlg('What is the filename?');
%user selects file name
fname=char(prompt);
%Eliminates the 23 lines of the header
rawdata=dlmread(fname,',',23,0);
%List creates and assigned the number 1 = Healthy, 2 = Infarcted
list={'Sham','Infarcted'};
isHealthy=listdlg('PromptString', {'What type of heart is the data from?'},...
    'SelectionMode','single','ListString',list);
time = rawdata(:,1); %Isolating time from the 1st column of the CSV
heartwaveform = rawdata(:,2); %Isolating the Heart Wave Form Data from the 2nd column of the CSV
%% Set sampling frequency
Fs = 250; % Hz
%% Design and Apply Low-Pass Filter to Raw Data Set
% Here, use designfilt() to design a lowpass filter. You may use an if
% statement to switch between different lowpass filters for infracter and
% healthy hearts.Adjust the passband frequency in order that you know that
% the signal is properly filtered. Also, plot the unfiltered and filtered
% signal to display the improvement made by filtering.

%First, plot the raw data for visual inspection and later comparison with
%lowpass filter
figure
plot(time,heartwaveform)
xlabel('Time (Seconds)')
ylabel('Pressure (mmHg)')
title('Raw Unfiltered Heart Condition Data')
%% Stop
if isHealthy == 1  %filter for healthy hearts with higher passband since there is less noise in the data set.
    LP = designfilt('lowpassfir','PassbandFrequency',12,...
    'StopbandFrequency',60,'StopbandAttenuation',70,'SampleRate',Fs);
    filtdata = filter(LP,heartwaveform); %Applying the filter to the data set
elseif isHealthy == 2  %filter for infarcted hearts with lower passband since there is more noise from this data set.
      LP = designfilt('lowpassfir','PassbandFrequency',8,...
    'StopbandFrequency',40,'StopbandAttenuation',60,'SampleRate',Fs);
    filtdata = filter(LP,heartwaveform); %Applying the filter to the data set
else
    disp('Invalid Heart State input. Please try again.') %If the Data is not Healthy or Infarcted Invalid State will be output
end
%% Plotting the Filtered Data

timedelay = grpdelay(LP); % find delay associated with low pass filter
disp(timedelay(1)); %Displaying the delay to the user
delay = timedelay(1); %Assigning the delay value to what the delay is for each filter
filtdata = filtdata(delay:end); % account for this delay in dataset
delaytime = time(1:length(filtdata));% Time of dataset accounting for time delay form filter
%% Plot of LowPass Data & Raw Data

figure
plot(time,heartwaveform, 'b-')  %Plot of original data in a blue line
hold on  %Overlay two plots on top of each other
plot(delaytime,filtdata, 'r-') %Plot of lowpass data in a red line
xlabel('Time (Seconds)','FontSize',16)
ylabel('Pressure (mmHg)','FontSize',16)
title('Low Pass vs. Raw Data Filtered Heart Condition','FontSize',18)
legend('Raw Data','Filtered Data','FontSize',12)
hold off %Deletes the overlay
%% Identify Local Maxima and Count for Heartbeats
% You may use previous code you have created in earlier labs to perfrom
% local maxima detection on the filtered signal. You should set your
% thresholds here. However, please remember to CHANGE YOUR VARIABLE
% NAMES. Tip - preallocate vectors of zeros to save time and processing power.
% It takes longer for the CPU to append to a vector than to change a vector
% value.

level = mean(filtdata); %Defining the level to be the mean of filtdata after the lowpass filtering
threshdata = false(size(filtdata)); %creating 0 and 1 as the code identifies points above the level (1) below the level (0)
thresdata(filtdata > level) = true; %Assigns a 1 if the filtdata is above the level for that point.
threshdataTrue = threshdata(threshdata == true);

peakdata = islocalmax(filtdata); %Finding local Maxima of the peakdata for heart count
maxlocal = find(peakdata); %Finding the max locations or indexes of the peakdata
disp(maxlocal) %Displaying these max locations to the user

for i = 1: length(maxlocal) %Creating a for loop to identify where the systolicpeaks index is less than the level 
    if maxlocal(i) < level
        maxlocal(i) = 0;
    end
end
 peakdata(peakdata==0) = []; %Allocating for the zeros
 maxlocal(maxlocal ==0) = []; 
 
relativepeaks = filtdata(maxlocal);

for i = 1: length(relativepeaks)
    if relativepeaks(i) < 50
      relativepeaks (i) = 0;
      maxlocal(i) = 0;
    end
end
 relativepeaks(relativepeaks==0) = []; %Allocating for the zeros
 maxlocal(maxlocal ==0) = []; 

maximumvalues = [];
for i = 1:length(maxlocal)
    maximumvalues(i) = filtdata(maxlocal(i));
end

%Plotting the max local data to visually verify how well the command is
%working. 

figure
plot(delaytime,filtdata, 'b-')
hold on
plot(delaytime(maxlocal),relativepeaks, 'or', 'MarkerSize',12)
xlabel('Time (Seconds)', 'FontSize',16)
ylabel('Pressure (mmHg)', 'FontSize',16)
title('Local Maxima - Filtered Heart Condition', 'FontSize',18)
%% Find peaks (Systolic)
% Use the findpeaks() function to find the peaks of the cleaned signal.
% Plot the peaks over the cleaned signal to prove that your threshold is
% correct.

%Thresholding the peaks of the filtered data
avgdata = mean(filtdata); %Using the avg of the data as the threshold/level value
[peaks,loc] = findpeaks(filtdata); %Utilizing the findpeaks function to find the peaks of each waveform period
systolicpeaks = (peaks); %Reassigning the peaks for clarity that these represent the systolicpeaks 
systolicloc = (loc); %Reassigning the variable for loc for clarity to represent this value as systolicloc
for i = 1: length(systolicpeaks) %Creating a for loop to identify where the systolicpeaks index is less than the level 
    if systolicpeaks(i) < avgdata
        systolicpeaks(i) = 0;
        systolicloc(i) = 0;
    end
end
 systolicpeaks(systolicpeaks==0) = []; %Allocating for the zeros
 systolicloc(systolicloc ==0) = []; %Allocating for the zeros

maxlocations = systolicloc; %Defining the max locations as the systolic locations since this is before the data is inverted.
disp(maxlocations); %Displaying these values to the user

%% Plotting the Systolic Pressure
% figure
% plot(delaytime(systolicloc),systolicpeaks, 'o', delaytime,filtdata);
% xlabel('Time(s)') 
% ylabel('Pressure (mmHg)')
% title('Systolic Peaks of Heart Pressure Waveform')

figure
plot(delaytime(systolicloc),systolicpeaks, 'o', delaytime,filtdata, 'MarkerSize',12); 
xlabel('Time(s)', 'FontSize',16)
ylabel('Pressure (mmHg)', 'FontSize',16)
title('Systolic Peaks of Heart Pressure Waveform', 'FontSize',18)
%% Find Minima (Diastolic) (inverted data set)
% Do the same as with the systolic, however invert the signal in order to
% find the diastolic minima occurance which now looks like a peak and thus you are able to use findpeaks(). Plot the occurances of the minima on
% the original filtered signal to prove that your threshold is correct.

%if and elseif statements created to find the correct diastolic peak
%locations. For most of the Sham data the best distance was +25 but for the
%Infarcted it was +30. There were two outliers where Sham 3 needed +30 and
%Infarct 1 needed +55. Additionally, since this is diastolic the inverse of
%the filtdata was used or -filtdata. 
if strcmp( fname, 'Sham 3.csv')
    [peaks1,loc1] = findpeaks(-filtdata,'MinPeakDistance',+30);
elseif isHealthy == 1
    [peaks1,loc1] = findpeaks(-filtdata,'MinPeakDistance',+25); %Sham 3 Data likes +30
elseif strcmp( fname, 'Infarct 1.csv')
     [peaks1,loc1] = findpeaks(-filtdata,'MinPeakDistance',+55);
elseif isHealthy == 2
    [peaks1,loc1] = findpeaks(-filtdata,'MinPeakDistance',+30);
 %55 for Infarct 1 and 30 for the rest of Infarct data
end

diastolicloc = []; %Creating the arrays for diastolic loc and disastolic peaks.
diastolicpeaks = [];
level = mean(-filtdata); %Creating another level but now with the inverse of the filtered data. 
k = 1; %assigning another indexing variable of k = 1
for i = 1:length(loc1) %Using a for loop for the locations of the diastolic peaks
    if peaks1(i)<level %Specifying if peaks index are less than the level to continue searching
        continue
    else
        diastolicloc(k) = loc1(i); %Otherwise, if peaks1 > level than diastolicloc is equal to the index of loc1
        diastolicpeaks(k) = peaks1(i); %Otherwise, if peaks1 > level than diastolicpeaks is equal to the index of peaks
        k = k+1; %keep cycling through the loop
    end
end

peakdata = islocalmax(-filtdata); %Finding localmax of the inverse filtdata
maxlocal = find(peakdata);
disp(maxlocal)

 minvalues = []; %Finding the min values of where the "max locations" are which is the min locations. 
for i = 1:length(maxlocal)
    minvalues(i) = -filtdata(maxlocal(i)); 
end
 minlocations = diastolicloc;
disp(minlocations); %Setting min locations as the diastoliclocs
%% Plotting Distolic Pressure Waveform

figure
plot(delaytime(diastolicloc),diastolicpeaks, 'o', delaytime,-filtdata, 'MarkerSize',12);  
xlabel('Time(s)', 'FontSize',16) 
ylabel('Pressure (mmHg)', 'FontSize',16)
title('Minima (Diastolic) of Heart Pressure Waveform', 'FontSize',18)
%% Maximum Developed Pressure
% Maximum developed pressure is the mean of the difference between the
% systolic and diastolic pressures. However, please remember that you may
% have more diastolic points than systolic points depending on when the
% recording starts during the heart beat! Use an if statement to adjust
% which systolic pressure to use (first recorded value or second)!

disp(maxlocations);
disp(minlocations);
%maxDP = average(systolic - diastolic)

if maxlocations > 10
  maxDP = mean((filtdata(maxlocations(1))-filtdata(minlocations))); %Identifiying if the maxDP calculation should occur at the 1st or 2nd based on how many maxlocations there are.
else 
   maxDP = mean((filtdata(maxlocations(2))-filtdata((minlocations)))); 
end

%% Maximum rate of pressure increase 
% Take the derivative of the filtered signal and find the peaks using the
% findpeaks() function once more. Please plot the differentiated signal and
% the peaks in order to prove that your are finding the peaks.

derivolt=diff(filtdata); %Taking the derivative of the filtdata 
level = 5; %Setting a level to compare the derivative of the filtdata to.
[peaks2,loc2] = findpeaks(derivolt); %Finding rhe peaks of derivative of the filtdata
Pmaxpeaks = (peaks2); %Assigning the variables with more clear names for future use in other code. 
Pmaxloc = (loc2); %Assigning the variables with more clear names for future use in other code. 
for i = 1:length(peaks2) %If peaks are less than the level then the index variables are equal to zero. 
    if Pmaxpeaks(i) < level
        Pmaxpeaks(i) = 0;
        Pmaxloc(i) = 0;
    end
end
Pmaxpeaks(Pmaxpeaks==0) = [];
Pmaxloc(Pmaxloc ==0) = []; 

derivolt=diff(filtdata);
i = 1; %Defining the index i = 1
for value = 1:length(Pmaxpeaks) %For loop to find the index values of the Pmaxpeaks from the derivolt value
    PmaxIndex(i) = find(derivolt == Pmaxpeaks(value));
    i = i+1;
end
%% Minimum rate of pressure increase
% Do the same as above, however you would apply the findpeaks() function to
% the inverted derivative vector to find the minima. Plot the minimum rates
% of pressure increase on the derivative graph to show that your threshold
% was adequate.

derivolt= diff(-filtdata);  %Taking the derivative of the -filtdata 
% avgdata = mean(-derivolt);
level = 5; %Setting a level to compare the derivative of the filtdata.
[peaks3,loc3] = findpeaks(derivolt); %Finding rhe peaks of derivative of the filtdata
Pminpeaks = (peaks3); %Assigning the variables with more clear names for future use in other code.
Pminloc = (loc3); %Assigning the variables with more clear names for future use in other code.
for i = 1: length(-Pminpeaks) %If peaks are less than the level then the index variables are equal to zero. 
    if Pminpeaks(i) < level
        Pminpeaks(i) = 0;
       Pminloc(i) = 0;
    end
end
Pminpeaks(Pminpeaks==0) = [];
Pminloc(Pminloc ==0) = []; 
 derivolt=diff(-filtdata);
i = 1;
for value = 1:length(Pminpeaks)  %For loop to find the index values of the Pminpeaks from the derivolt value
    PminIndex(i) = find(derivolt == Pminpeaks(value));
    i = i+1;
end
%% Validation of minima dp/dt and minima
% Plot the original filtered signal, but now with where the max and minimum
% change in pressures noted. Best way to do so is to take the occurances of
% the minima and maxima (which should be samples) and plot it against the
% original signal values at those occurances(aka samples).
% 

figure
plot(delaytime,filtdata)
hold on
plot(delaytime(Pmaxloc),filtdata(PmaxIndex), 'o', 'MarkerSize',12)
plot(delaytime(Pminloc),filtdata(PminIndex), 'bo', 'MarkerSize',12)
xlabel('Time(s)','FontSize',16) 
ylabel('Pressure (mmHg)','FontSize',16)
title('dp/dt of Heart Pressure Waveform','FontSize',18)
%% Diastolic Time Constant
% % Find the diastolic time constant over a time range as noted in lecture.
% % Please see the pressureerror and the pressureeqn Matlab functions and
% % scripts provided by the Professor. Remeber to acount for if the first
% % diastolic value occurs after the first minimum dp/dt value (use an if
% % loop). Plot to show how well the curve fits the original signal (or if
% % it works at all!) This is the hardest part of the final project, so don't
% % get discouraged if you have issues in this section.

overalltime = []; %Creating array for the overalltime values
overallmag = []; %Creating array for the overallmag values
tao_estimate = []; %Creating array for the tao estimate 
i = 1; %Setting index to 1 
if diastolicloc(i) < Pminloc(i) %Making sure these arrays are the same size before finding the region
diastolicloc = diastolicloc(2:end);
end
minima = diastolicpeaks; %Setting the minima to the diastolic peaks
voltage = filtdata; %Voltage is defined as the filtdata
for i = 1:length(minima)-1 %Using this for loop to ensure the array's are the same sizes and are compatible. 
    if Pminloc(i) < diastolicloc(i)
        region = (Pminloc(i): diastolicloc(i)); %specifying the region for the tao constant to run through which is the lower portion of the graph
    elseif Pminloc(i) > diastolicloc(i)
        region = (diastolicloc(i): Pminloc(i));
    elseif Pminloc(i) == 0 || diastolicloc(i) == 0
      Pminloc(i) = [ ] ;
      diastolicloc(i) = [ ];
    end
 timex = time(region);
overalltime = [overalltime timex' NaN]; %NaN to filter out NaN and only plot the time of the diastolic relaxation pressure. 
% Define starting point
% [Po,P1,tau]
P0 = [1 1 1];
% Lower bounds
lb = [0 0 .00001];
% Upper bounds
ub = [Inf Inf Inf];

anonfunc = @(P) pressureerror(P,timex,voltage(region));
% fitted_pressure = pressureeqn(Pest,timex);
Pest = fmincon(anonfunc,P0,[],[],[],[],lb,ub);

tao_estimate = [tao_estimate Pest(3)];

fin = pressureeqn(Pest,timex);
overallmag = [overallmag fin' NaN]; %NaN to filter out NaN and only plot the magnitude of the diastolic relaxation pressure. 
end

figure
plot(delaytime,filtdata,'b-') %Plotting the filtered data
hold on %Overlay on filtered data
plot(overalltime,overallmag,'r-')  %Plotting the diastolic pressure decrease overlayed onto the filtered data
xlabel('Time(s)','FontSize',16) 
ylabel('Pressure (mmHg)','FontSize',16)
title('Diastolic Time Constants of Heart Pressure Waveform','FontSize',14)
%% Final Display of all Parameters to perform t and p tests on 
% %Finally display your average diastolic and systolic pressures, your
% %maximum deveoped pressure, your tau, and your maximum and minimum dp/dt
% %values for the user to see on the command window. And thats it :D
% 
% 
disp(['Heartrate: ', num2str(((length(systolicpeaks))/delaytime(length(delaytime)))*60), ' Beats/min ']) 

disp(['Average Diastolic Pressure: ', num2str(-mean(minvalues)), ' mmHg      ', 'Average Systolic Pressure: ', num2str(mean(systolicpeaks)), ' mmHg '])
% 
disp(['Maximum Developed Pressure: ', num2str(mean(maxDP)), ' mmHg '])
% 
disp(['Tao: ', num2str(mean(tao_estimate)), ' Units N/A'])
% 
disp(['Maximum dp/dt Value: ', num2str(max(Pmaxpeaks)), ' mmHg/s ','        Minimum dp/dt Value: ', num2str(min(-Pminpeaks)), ' mmHg/s '])