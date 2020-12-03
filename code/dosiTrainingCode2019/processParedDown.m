%This is a script to validate the pared down version of dosigui processing.
%The motivation was that dosigui performs some operations to the raw data
%prior to fitting which we thought might be important for getting accurate
%OPs. We wanted to perform those same operations to transmission
%measurements, so I replicated the dosigui processing in a pared down
%format. This might also be a helpful training tool and is a more
%easy-to-edit version of dosigui that provides the same output. Note:
%dosigui has a ton of other options, this script only replicates options
%selected for the "Streaming Chromophores" paper (Yanyu Zhao, Mattew B.
%Applegate, Raeef Istfan, Ashvin Pande, and Darren Roblyer, 
%"Quantitative real-time pulse oximetry with ultrafast frequency-domain 
%diffuse optics and deep neural network processing," Biomed. Opt. Express 
%9, 5997-6008 (2018)

%Location of ASC files containing raw amplitudes and phases
dataDir = '.\testData\181030\';
%ACRIN9 is the calibration phantom
%Get file names of all the measurements of ACRIN9
fnames = dir(fullfile(dataDir,'ACRIN9*miniLBS.asc'));
n=1.4; %Index of refraction
%Preallocate a cell array containing the locations of all of the
%calibration measurements
f = cell(1,length(fnames));
%Fill the array
for i = 1:length(fnames)
    f{i} = fullfile(dataDir,fnames(i).name);
end
%Average together all of the calibration measurements
avgCal=averageASCData(f,6,[.03,.3]);
%Calculate the system response based on the calibration data
calDat=getCalibrationData(avgCal,'ACRIN9.txt',n);

%Get file names of all the samples to measure
fnames2 = dir(fullfile(dataDir,'*-miniLBS.asc'));
f2 = cell(1,length(fnames2));
for i = 1:length(fnames2)
    f2{i} = fullfile(dataDir,fnames2(i).name);
end

tic
%Process each file
for i = 1:length(f2)
   
    thisC = f2(i); %Name of current file. If you want to average raw data make this a cell array with all the datasets to be averaged
    rawDat = averageASCData(thisC,6,[0.03,.3]); %Read raw data
    
    %Replace any zeros in the amplitude with the previous value. Will fail
    %if the first amplitude is zero.
    if ~isempty(find(rawDat.AC==0,1))
        warning('Amplitude value of zero replaced')
        idX = find(rawDat.AC == 0);
        rawDat.AC(idX) = rawDat.AC(idX-1);
    end
    %Get the name of the sample
    splitsies = strsplit(f2{i},'\');
    file_name = splitsies{end};
    splitsies2 = strsplit(file_name,'.');
    name = splitsies2{1};
    fprintf('Working on scan %d of %d: %s\n',i,length(f2),name)
    
    %Calibrate the raw data
    cal=calibratePD_FDPM(calDat,rawDat);
    %Transform from amplitude and phase to real and imaginary parts
    realAtEnd = cal.AC .* cos(cal.phase);
    imagAtEnd = cal.AC .* sin(cal.phase);
    %Concatenate the data
    DATA = [realAtEnd;imagAtEnd];
    %Calculate the weighting function
    ERR_REAL1=sqrt((cal.damp.*cos(cal.phase)).^2 + (cal.dphi.*cal.AC.*sin(cal.phase)).^2);
    ERR_IMAG1=sqrt((cal.damp.*sin(cal.phase)).^2 + (cal.dphi.*cal.AC.*cos(cal.phase)).^2);
    WT_REIM1 = abs([1./ERR_REAL1;1./ERR_IMAG1]);
    
    %Replace any NaNs in the weighting function with 1s
    if sum(isnan(WT_REIM1(:))) ~= 0
        warning('NaN in the weighting function replaced with 1')
        idxs = find(isnan(WT_REIM1(:))==1);
        WT_REIM1(idxs) = 1;
    end
    %Apply the weight
    YDATA = DATA.*WT_REIM1;
    %Fit for mua and musp
    fitDat(i) = fitMu(YDATA,WT_REIM1,cal.freq,n,cal.dist,name);
  
end
toc
%%%%%%%%%%%%This part compares the Training code with DOSIGUI
%Read in ASC file
guiFname = '.\PROCESSED\testData_181030\pttestData_181030_allFreqs__dBMU.asc';
x = read_dbMU(guiFname);

%Organize data more sensibly
guiDat = struct('name','','wv',[],'mua',[],'musp',[]);
iter = 1;
iterFast = 1;
for i = 1:length(fitDat)
    guiDat(iter).name = x(iterFast).name;
    for j = 1:6
        guiDat(iter).wv(j) = x(iterFast).wv;
        guiDat(iter).mua(j) = x(iterFast).mua;
        guiDat(iter).musp(j) = x(iterFast).mus;
        iterFast = iterFast+1;
    end
    iter = iter+1;
end

%Make sure the pared down and GUI data are in the same order
pd = struct2table(fitDat);
gd = struct2table(guiDat);
sortedPD = sortrows(pd,'name');
sortedGD = sortrows(gd,'name');
sortedPD = table2struct(sortedPD);
sortedGD = table2struct(sortedGD);

%Extract the mua and musp from each method
clear muaVec musVec musVecGui muaVecGui;
muaVecGUI = [];
musVecGUI = [];
muaVec = [];
musVec = [];
iter = 1;
for j= 1:length(guiDat)
    muaVec = [muaVec,sortedPD(iter).mua];
    musVec = [musVec,sortedPD(iter).mus];
    
    musVecGUI = [musVecGUI,sortedGD(iter).musp];
    muaVecGUI = [muaVecGUI,sortedGD(iter).mua];
    iter=iter+1;
end
%Fitting returns a "goodness of fit" statistic for discriminating good and
%bad fits
chiSq = [fitDat(:).chi];
goodFits = find(chiSq<50);

%%Plot the two processing results
figure
plot(muaVecGUI)
hold on
plot(muaVec)
xlabel('Sample number')
ylabel('\mu_a (mm^{-1})')

figure
plot(musVecGUI)
hold on
plot(musVec)
xlabel('Sample number')
ylabel('\mu_s'' (mm^{-1})')

figure
plot(muaVecGUI,muaVec,'o')
hold on
plot([0,.25],[0,.25],'k--')
xlabel('\mu_a DOSIGUI (mm^{-1})')
ylabel('\mu_a Training Code (mm^{-1})')

figure
plot(musVecGUI,musVec,'o')
hold on
plot([0,3],[0,3],'k--')
xlabel('\mu_s'' DOSIGUI (mm^{-1})')
ylabel('\mu_s'' Training (mm^{-1})')

%Only the "good" fits

figure
plot(muaVecGUI(goodFits),muaVec(goodFits),'o')
hold on
plot([0,.2],[0,.2],'k--')
xlabel('\mu_a DOSIGUI (mm^{-1})')
ylabel('\mu_a Training Code (mm^{-1})')

figure
plot(musVecGUI(goodFits),musVec(goodFits),'o')
hold on
plot([0,1],[0,1],'k--')
xlabel('\mu_s'' DOSIGUI (mm^{-1})')
ylabel('\mu_s'' Training (mm^{-1})')