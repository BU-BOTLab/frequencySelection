%compareInversionTiming.m
%
%PURPOSE: Measure the inversion time for the three different inversion
%models tested. That is, the optimized iterative model, the defualt
%iterative model, and the LUT model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Setup some constants needed for the simulations

%Set random seed for reproducibility
rng(14850)
%Add forward and inverse model to path
addpath(genpath('dosiTrainingCode2019'))
addpath(genpath('utilities'))
load('../generatedData/OPSet10k.mat')
%Use the all digital noise model
dataDir='../experimentalData/allDigitalSystem/Gen2_30mm/8192/200614';
sampName = 'b4h2*.asc';
%Acquisition parameters for noise model samples
numDiodes = 6;     %Number of lasers used in the experimental data
numMeasFreqs = 70; %number of frequencies used in experimental data
nind = 1.4;        %Index of refraction for forward model

%Parameters to use for simulations
sepRange= [30]; %Source-detector separation
wt =0;          %Weight parameter for forward model
reImFlag = 1;   %Do all tests with real/imaginary parts
%Only test a subset of the total number of OPs for speed
numTrials = 300;

%Allocate memory for OPs
recMuaLUT = zeros(1,numTrials);
recMusLUT = zeros(1,numTrials);
recMuaHDLUT = zeros(1,numTrials);
recMusHDLUT = zeros(1,numTrials);
recMuaDefault = zeros(1,numTrials);
recMusDefault = zeros(1,numTrials);
recMuaOpt = zeros(1,numTrials);
recMusOpt = zeros(1,numTrials);

%Frequencies to test inversion speed with
freqList = {50,50:99,50:149,50:199,50:249,50:299,50:349,50:399,50:449, 50:499};
numFreqs = [1,50,100,150,200,250,300,350,400,450];
%Allocate memory for recovered OPs
saveMuaTrue = zeros(numTrials,length(freqList),length(sepRange));
saveMusTrue = zeros(numTrials,length(freqList),length(sepRange));
saveMuaRecLUT = zeros(numTrials,length(freqList),length(sepRange));
saveMusRecLUT = zeros(numTrials,length(freqList),length(sepRange));
saveMuaRecDefault = zeros(numTrials,length(freqList),length(sepRange));
saveMusRecDefault = zeros(numTrials,length(freqList),length(sepRange));
saveMuaRecOpt = zeros(numTrials,length(freqList),length(sepRange));
saveMusRecOpt = zeros(numTrials,length(freqList),length(sepRange));

%Allocate memory for timing variables
timingLUT = zeros(1,length(freqList));
timingDefault = zeros(1,length(freqList));
timingOpt = zeros(1,length(freqList));

%%
%Test the LUT method
%Range of OPs for the LUT
muaAll = [0.001,0.05];
musAll = [.1,3];
%Grid of OPs to calculate LUT
allMua = linspace(muaAll(1),muaAll(2),500);
allMus = linspace(musAll(1),musAll(2),500);

sepIdx = 1; %SD Separation index
waveIdx = 6; %Wavelength index
filterLen = 5; %Moving average filter to smooth the noise model
%Iterate over each frequency set
for p = 1:length(freqList)
    fa= freqList{p}; %This frequency set
    %Generate or load the necessary LUT
    LUT = generateLUT(allMua,allMus,sepRange(sepIdx),fa);
    %Load the noise model
    [~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,fa,filterLen);
    tic %start timing
    for i=1:numTrials
        %Load the current OP pair
        randMua = trueMua(i);
        randMus = trueMus(i);
        %Generate the raw data
        rawDat = p1seminf_mba([randMua,randMus],fa',nind,sepRange(sepIdx),wt,reImFlag);
        %extract real and imaginary parts
        rp = rawDat(1:length(fa));
        ip = rawDat(length(fa)+1:end);
        %Convert to amp/phase
        amp = sqrt(rp.^2 + ip.^2);
        phase = unwrap(atan2(ip,rp));
        %Add noise to each frequency
        noiseAmp = zeros(size(amp));
        noisePhase = zeros(size(phase));
        for thisFreq = 1:length(amp)
            noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
            noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
        end
       
        %Convert back to real/imaginary
        noiseReal = noiseAmp .* cos(noisePhase);
        noiseImag = noiseAmp .* sin(noisePhase);
        %Run inverse model
        [recMuaLUT(i),recMusLUT(i)] = getLUTOPs(LUT,[noiseReal';noiseImag'],0);
    end
    %Save the amount of time it took to run all of the inversions
    timingLUT(p) = toc;
    %Save the recovered data (not used for this experiment)
    saveMuaTrue(:,p,sepIdx) = trueMua(1:numTrials);
    saveMusTrue(:,p,sepIdx) = trueMus(1:numTrials);
    saveMuaRecLUT(:,p,sepIdx) = recMuaLUT;
    saveMusRecLUT(:,p,sepIdx) = recMusLUT;
end

%%
%Test the optimized iterative method
for p = 1:length(freqList)
    fprintf('Working freq set %d optimized settings\n',p)
    fa= freqList{p}; %This frequency set
    %Get noise model standard deviations
    [~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,fa,5);
    %Start timing
    tic
    for i=1:numTrials
        %Get this OP pair
        randMua = trueMua(i);
        randMus = trueMus(i);
        %Calculate raw data
        rawDat = p1seminf_mba([randMua,randMus],fa',nind,sepRange(sepIdx),wt,reImFlag);
        %extract real and imaginary parts
        rp = rawDat(1:length(fa));
        ip = rawDat(length(fa)+1:end);
        %Calculate amp/phase
        amp = sqrt(rp.^2 + ip.^2);
        phase = unwrap(atan2(ip,rp));
        %Add noise to each frequency
        noiseAmp = zeros(size(amp));
        noisePhase = zeros(size(phase));
        for thisFreq = 1:length(amp)
            noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
            noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
        end
        %Convert back to real/imaginary
        noiseReal = noiseAmp .* cos(noisePhase);
        noiseImag = noiseAmp .* sin(noisePhase);
        %Run inverse model
        recOPOpt = fitMu_newOpts([noiseReal;noiseImag],wt, fa', nind, sepRange(sepIdx),'test');
        %Save recovered OP
        recMuaOpt(i) = recOPOpt.mua;
        recMusOpt(i) = recOPOpt.mus;
    end
    %Get timing data
    timingOpt(p) = toc;
    saveMuaRecOpt(:,p,sepIdx) = recMuaOpt;
    saveMusRecOpt(:,p,sepIdx) = recMusOpt;
end

%%
%Test the default method
for p = 1:length(freqList)
    fprintf('Working freq set %d default settings\n',p)
    fa= freqList{p}; %This OP set
    %Get noise model for these frequencies
    [~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,fa,5);
    tic %Start timing
    for i=1:numTrials
        %Get this OP pair
        randMua = trueMua(i);
        randMus = trueMus(i);
        %Calculate forward data
        rawDat = p1seminf_mba([randMua,randMus],fa',nind,sepRange(sepIdx),wt,reImFlag);
        %extract real and imaginary parts
        rp = rawDat(1:length(fa));
        ip = rawDat(length(fa)+1:end);
        %Calculate amplitude/phase
        amp = sqrt(rp.^2 + ip.^2);
        phase = unwrap(atan2(ip,rp));
        %Add noise to each frequency
        noiseAmp = zeros(size(amp));
        noisePhase = zeros(size(phase));
        for thisFreq = 1:length(amp)
            noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
            noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
        end
       
        %Convert back to real/imaginary
        noiseReal = noiseAmp .* cos(noisePhase);
        noiseImag = noiseAmp .* sin(noisePhase);
        %Run inverse model
        recOPDefault = fitMu([noiseReal;noiseImag],wt, fa', nind, sepRange(sepIdx),'test');
        recMuaDefault(i) = recOPDefault.mua;
        recMuaDefault(i) = recOPDefault.mus;
    end
    %Stop timing
    timingDefault(p) = toc;
    saveMuaRecDefault(:,p,sepIdx) = recMuaDefault;
    saveMusRecDefault(:,p,sepIdx) = recMusDefault;
end
%%
%Test the HD-LUT method
%Range of OPs for the LUT
muaAll = [0.001,0.05];
musAll = [.1,3];
%Grid of OPs to calculate LUT
allMua = linspace(muaAll(1),muaAll(2),2000);
allMus = linspace(musAll(1),musAll(2),2000);
%Frequencies to test inversion speed with
freqList = {50,50:54,50:59,50:69};
numFreqsHD = [1,5,10,20];
sepIdx = 1; %SD Separation index
waveIdx = 6; %Wavelength index
filterLen = 5; %Moving average filter to smooth the noise model
saveMuaRecHDLUT = zeros(numTrials,length(freqList),length(sepRange));
saveMusRecHDLUT = zeros(numTrials,length(freqList),length(sepRange));
timingHDLUT = zeros(1,length(freqList));
%Iterate over each frequency set
for p = 1:length(freqList)
    fa= freqList{p}; %This frequency set
    %Generate or load the necessary LUT
    LUT = generateLUT(allMua,allMus,sepRange(sepIdx),fa);
    %Load the noise model
    [~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,fa,filterLen);
    tic %start timing
    for i=1:numTrials
        %Load the current OP pair
        randMua = trueMua(i);
        randMus = trueMus(i);
        %Generate the raw data
        rawDat = p1seminf_mba([randMua,randMus],fa',nind,sepRange(sepIdx),wt,reImFlag);
        %extract real and imaginary parts
        rp = rawDat(1:length(fa));
        ip = rawDat(length(fa)+1:end);
        %Convert to amp/phase
        amp = sqrt(rp.^2 + ip.^2);
        phase = unwrap(atan2(ip,rp));
        %Add noise to each frequency
        noiseAmp = zeros(size(amp));
        noisePhase = zeros(size(phase));
        for thisFreq = 1:length(amp)
            noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
            noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
        end
       
        %Convert back to real/imaginary
        noiseReal = noiseAmp .* cos(noisePhase);
        noiseImag = noiseAmp .* sin(noisePhase);
        %Run inverse model
        [recMuaHDLUT(i),recMusHDLUT(i)] = getLUTOPs(LUT,[noiseReal';noiseImag'],0);
    end
    %Save the amount of time it took to run all of the inversions
    timingHDLUT(p) = toc;
    %Save the recovered data (not used for this experiment)
    saveMuaTrue(:,p,sepIdx) = trueMua(1:numTrials);
    saveMusTrue(:,p,sepIdx) = trueMus(1:numTrials);
    saveMuaRecHDLUT(:,p,sepIdx) = recMuaLUT;
    saveMusRecHDLUT(:,p,sepIdx) = recMusLUT;
end
%%
save('../generatedData/inversionTimingData.mat','numTrials','timingHDLUT','timingLUT','timingOpt','timingDefault','numFreqs','numFreqsHD')

%%
load('../generatedData/inversionTimingData.mat')
invPerSec = numTrials./[timingLUT; timingOpt; timingDefault];
invPerSecHD = numTrials./timingHDLUT;

figure
l1=semilogy(numFreqs, invPerSec(1,:)','-o');
hold on
l2=semilogy(numFreqs, invPerSec(2,:)','-o');
l3=semilogy(numFreqs, invPerSec(3,:)','-o');
l4 = semilogy(numFreqsHD, invPerSecHD,'-o');
c1=get(l1,'color');
c2=get(l2,'color');
c3=get(l3,'color');
set(l1,'markerfacecolor',c1)
set(l2,'markerfacecolor',c2)
set(l3,'markerfacecolor',c3)
xlabel('Number of frequencies')
ylabel('Inversions/second')
title('Inversion speed comparison')
legend('LUT','Optimized','Default','HDLUT')
set(gca,'fontsize',18)
print('../plots/inversionTimingPlot.png','-dpng')