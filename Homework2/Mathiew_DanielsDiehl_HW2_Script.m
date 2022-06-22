f = csvread('vpressure_45Hz.csv');
f2 = csvread('vpressure_200Hz.csv');


intf = interp1(1:numel(f), f, 1:0.225:numel(f));

%RMSE Calculation and resize f2 data
fnorm = f2(1:end-4);
rmse = sqrt(mean(intf-fnorm).^2);

% Filter Design
diff1 = designfilt('differentiatorfir', 'FilterOrder', 5);
diff2 = designfilt('differentiatorfir', 'FilterOrder', 13);
diff3 = designfilt('differentiatorfir', 'FilterOrder', 19);

% filt filt application 45 Hz
filt1 = (filtfilt(diff1, f))*45;
filt2 = (filtfilt(diff2, f))*45;
filt3 = (filtfilt(diff3, f))*45;

% Filtfilt application 200Hz
filt4 = (filtfilt(diff1, f2))*200;
filt5 = (filtfilt(diff2, f2))*200;
filt6 = (filtfilt(diff3, f2))*200;

% plotting
x1 = linspace(0,5.91, 266);
y1 = filt1;
y2 = filt2;
y3 = filt3;

figure(1)
plot(x1, y1, 'r')
hold on
plot(x1, y2, 'g')
hold on
plot(x1, y3, 'b')
title('dP/dt 45Hz over Time')
xlabel('Time (s)')
ylabel('dP/dt')
hold off

figure(2)
x2 = linspace(0, 5.91, 1182);
y4 = filt4;
y5 = filt5;
y6 = filt6;
plot(x2, y4, 'r')
hold on
plot(x2, y5, 'g')
hold on
plot(x2, y6, 'b')
title('dP/dt 200Hz over Time')
xlabel('Time (s)')
ylabel('dP/dt')
hold off










