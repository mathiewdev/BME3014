%% statistical analysis
filename = input("Please input the final file name: ", 's');
B = readmatrix(filename)
heartbeat = B(:,2);
distolicP = B(:,3);
systolicP = B(:,4);
maxDevP = B(:,5);
tao = B(:,6);
maxDPDT = B(:,7);
minDPDT = B(:,8);

%% calculation for sham data 
    SavgHR = mean(heartbeat(1:6, :));
    SavgDP = mean(distolicP(1:6, :));
    SavgSP = mean(systolicP(1:6, :));
    SavgMDP = mean (maxDevP(1:6, :));
    SavgT = mean (tao(1:6, :));
    SavgMaDPDT = mean (maxDPDT(1:6, :));
    SavgMiDPDT = mean (minDPDT(1:6, :));
    
    SdevHR = std(heartbeat(1:6, :));
    SdevDP = std(distolicP(1:6, :));
    SdevSP = std(systolicP(1:6, :));
    SdevMDP = std (maxDevP(1:6, :));
    SdevT = std (tao(1:6, :));
    SdevMaDPDT = std (maxDPDT(1:6, :));
    SdevMiDPDT = std (minDPDT(1:6, :));
 
 disp('AVERAGE OF SHAM DATA W/ STANDARD DEVIATION')
 disp(['Heart rate ', num2str(SavgHR), ' +/- ',num2str(SdevHR)])
 disp(['Diast P ', num2str(SavgDP), ' +/- ',num2str(SdevDP)])
 disp(['Syst P ', num2str(SavgSP), ' +/- ',num2str(SdevSP)])
 disp(['Max dev P ', num2str(SavgMDP), ' +/- ',num2str(SdevMDP)])
 disp(['tau ', num2str(SavgT), ' +/- ',num2str(SdevT)])
 disp(['max dpdt ', num2str(SavgMaDPDT), ' +/- ',num2str(SdevMaDPDT)])
 disp(['min dpdt ', num2str(SavgMiDPDT), ' +/- ',num2str(SdevMiDPDT)])
 
 %% Infarcted data
  
    IavgHR = mean(heartbeat(7:12, :));
    IavgDP = mean(distolicP(7:12, :));
    IavgSP = mean(systolicP(7:12, :));
    IavgMDP = mean (maxDevP(7:12, :));
    IavgT = mean (tao(7:12, :));
    IavgMaDPDT = mean (maxDPDT(7:12, :));
    IavgMiDPDT = mean (minDPDT(7:12, :));
    
    IdevHR = std(heartbeat(7:12, :));
    IdevDP = std(distolicP(7:12, :));
    IdevSP = std(systolicP(7:12, :));
    IdevMDP = std (maxDevP(7:12, :));
    IdevT = std (tao(7:12, :));
    IdevMaDPDT = std (maxDPDT(7:12, :));
    IdevMiDPDT = std (minDPDT(7:12, :));
 
 disp('AVERAGE OF INFARCTED DATA W/ STANDARD DEVIATION')
 disp(['Heart rate ', num2str(IavgHR), ' +/- ',num2str(IdevHR)])
 disp(['Diast P ', num2str(IavgDP), ' +/- ',num2str(IdevDP)])
 disp(['Syst P ', num2str(IavgSP), ' +/- ',num2str(IdevSP)])
 disp(['Max dev P ', num2str(IavgMDP), ' +/- ',num2str(IdevMDP)])
 disp(['tau ', num2str(IavgT), ' +/- ',num2str(IdevT)])
 disp(['max dpdt ', num2str(IavgMaDPDT), ' +/- ',num2str(IdevMaDPDT)])
 disp(['min dpdt ', num2str(IavgMiDPDT), ' +/- ',num2str(IdevMiDPDT)])  
 
 
 %% T-tests
 disp('T-Test')
  disp('Heart rate')
 [h1, p1] = ttest2(heartbeat(1:6),heartbeat(7:12))
 disp('Diastp');
  [h2, p2] = ttest2(distolicP(1:6) , distolicP(7:12),'Alpha',0.05)
 disp('Systp ');
 [h3, p3] = ttest2(systolicP(1:6) , systolicP(7:12),'Alpha',0.05)
 disp('maxDevP ');
  [h4, p4] = ttest2(maxDevP(1:6) , maxDevP(7:12),'Alpha',0.05)
 disp('tao ')
   [h5, p5] = ttest2(tao(1:6) , tao(7:12),'Alpha',0.05)
 disp('maxDPDT ');
  [h6, p6] = ttest2(maxDPDT(1:6) , maxDPDT(7:12),'Alpha',0.05)
 disp('minDPDT ');
   [h7, p7] = ttest2(minDPDT(1:6) , minDPDT(7:12),'Alpha',0.05)
