%distilled bar plots for Darren
rng(14850)
%Add forward and inverse model to path
addpath(genpath('dosiTrainingCode2019'))
load('../generatedData/OPSet10k.mat')
%Data directories containing raw data for building noise model
dataDir10='../experimentalData/allDigitalSystem/Gen2_10mm/8192/200614';
dataDir20='../experimentalData/allDigitalSystem/Gen2_20mm/8192/200614';
dataDir30='../experimentalData/allDigitalSystem/Gen2_30mm/8192/200614';
%Samples to use for noise models
sampName =  'bpav4*.asc';

%Acquisition parameters for noise model samples
numDiodes = 6;
numMeasFreqs = 70;
nind = 1.4;

%Parameters to use for simulations
sepRange = [30];
wt =0;
reImFlag = 1;
numTrials = 3000;

%Range of OPs
muaAll = [0.001,0.05];
musAll = [.1,3];

%Allocate memory for OPs
% trueMua = zeros(1,numTrials);
% trueMus = zeros(1,numTrials);

numReps = 50;

%Frequencies to explore
% freqList = {[50],[60],[70],[80],[90],[100],[110],[120],[130],[140],[150],[160],[170],[180],[190],...
%     [200],[210],[220],[230],[240],[250],[260],[270],[280],[290],[300],[310],[320],[330],[340],...
%     [350],[360],[370],[380],[390],[400],[410],[420],[430],[440],[450],[460],[470],[480],[490]};

freqList = {[50],[70],[90],[110],[130],[150],[170],[190],...
    [210],[230],[250],[270],[290],[310],[330],...
    [350],[370],[390],[410],[430],[450],[470],[490]};
%Number of frequencies in each test set
numFreqs = ones(1,length(freqList));
%Maximum number of frequencies
maxFreqs = max(numFreqs);
timVec = zeros(1,length(numFreqs));
it = 1;
wavelength = 850;
waveIdx = 6;
sep = 1;
recMua = zeros(1,numTrials);
recMus = zeros(1,numTrials);
%Allocate memory for MAT file variables
saveMuaTrue = zeros(numTrials,length(freqList),length(sepRange));
saveMusTrue = zeros(numTrials,length(freqList),length(sepRange));
saveMuaRec = zeros(numTrials,length(freqList),length(sepRange));
saveMusRec = zeros(numTrials,length(freqList),length(sepRange));
%Iterate over all separation ranges
dataDir = dataDir20;

for f=1:length(freqList)
    fprintf("\tworking on num freqs %d of %d\n", f,length(freqList));              
    fa = freqList{f};
    %Calculate noise model
    [~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,fa,5);
    %Pick a bunch of OP pairs to test
    tic
    parfor i=1:numTrials
   
        randMua = trueMua(i);
        randMus = trueMus(i);
        %Generate raw data
        rawDat = p1seminf_mba([randMua,randMus],fa',nind,sepRange(sep),wt,reImFlag);
        %extract real and imaginary parts
        rp = rawDat(1:length(fa));
        ip = rawDat(length(fa)+1:end);
        %Convert to amp/phase
        amp = sqrt(rp.^2 + ip.^2);
        phase = unwrap(atan2(ip,rp));

    %Allocate memory for 'equal data scheme'
        runSumAmp = zeros(size(rp));
        runSumPhase = zeros(size(rp));

        %To ensure that all conditions have 'equal data'
        %simulate multiple measurements for conditions with
        %small numbers of frequencies and average the
        %amplitude and phase
        for rep=1:numReps
           %Allocate memory for noisy measurement
           noiseAmp = zeros(size(amp));
           noisePhase = zeros(size(phase));
           %Uncomment for correlated amplitude
%                            if wavelength == 850
%                                thisRandom = randn(1);
%                            end
           %Add noise to each frequency
           for thisFreq = 1:length(amp)
               %Uncomment for correlated amplitude
%                                if wavelength == 850
%                                    noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*thisRandom*ampfit(thisFreq,waveIdx);
%                                else
%                                noiseAmp(thisFreq) = amp(thisFreq) + amp(thisFreq).*randn(1)*noiseLevel;%ampfit(thisFreq,waveIdx);
%                                end
               %Add noise to amp/phase
               noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
               noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
           end

           runSumAmp = runSumAmp + noiseAmp;
           runSumPhase = runSumPhase + noisePhase;
        end
        %Average amp/phase
        avgAmp = runSumAmp/numReps;
        avgPhase = runSumPhase/numReps;
        %Convert back to real/imaginary
        noiseReal = avgAmp .* cos(avgPhase);
        noiseImag = avgAmp .* sin(avgPhase);
        recOP = fitMu_newOpts([noiseReal;noiseImag],wt, fa', nind, sepRange(1),'test');
        %Save recovered OP pair
        recMua(i) = recOP.mua;
        recMus(i) = recOP.mus;

    end
    timVec(it) = toc;
    it = it+1;
    %Save values in array
    saveMuaTrue(:,f,sep) = trueMua(1:numTrials);
    saveMuaRec(:,f,sep) = recMua;
    saveMusTrue(:,f,sep) = trueMus(1:numTrials);
    saveMusRec(:,f,sep) = recMus;

end 


save(sprintf('../generatedData/digital_%dnm_optimalFreq_30mm.mat',wavelength),'saveMuaTrue','saveMusTrue', 'saveMuaRec','saveMusRec','freqList');



%Iterate over all separation ranges
rng(14850);
%Data directories containing raw data for building noise model
dataDir = '../experimentalData/networkAnalyzerSystem/mDOSI_Probe2020_Drift/200304';
%Samples to use for noise models
sampName =  'INO9*miniLBS.asc';
%Acquisition parameters for noise model samples
numDiodes = 6;
numMeasFreqs = 401;
nind = 1.4;

%Parameters to use for simulations
sep = 28;
wt =0;
reImFlag = 1;
numTrials = 3000;

%Allocate memory for MAT file variables
saveMuaTrue = zeros(numTrials,length(freqList),1);
saveMusTrue = zeros(numTrials,length(freqList),1);
saveMuaRec = zeros(numTrials,length(freqList),1);
saveMusRec = zeros(numTrials,length(freqList),1);
recMua = zeros(1,numTrials);
recMus = zeros(1,numTrials);

for f=1:length(freqList)
%        LUT = generateLUT(allMua,allMus,sep,fa);
    fprintf("\tworking on num freqs %d of %d\n", f,length(freqList));              
    fa = freqList{f};
    %Calculate noise model
    [~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,fa,15);
    %Pick a bunch of OP pairs to test
    parfor i=1:numTrials
        %fprintf('%d\n',i);
        %Save the true OP pair
        randMua=trueMua(i);
        randMus=trueMus(i);
        %Generate raw data
          rawDat = p1seminf_mba([randMua,randMus],fa',nind,sep,wt,reImFlag);
        %extract real and imaginary parts
        rp = rawDat(1:length(fa));
        ip = rawDat(length(fa)+1:end);
        %Convert to amp/phase
        amp = sqrt(rp.^2 + ip.^2);
        phase = unwrap(atan2(ip,rp));

        %Allocate memory for 'equal data scheme'
        runSumAmp = zeros(size(rp));
        runSumPhase = zeros(size(rp));

        %To ensure that all conditions have 'equal data'
        %simulate multiple measurements for conditions with
        %small numbers of frequencies and average the
        %amplitude and phase
        for rep=1:numReps
           %Allocate memory for noisy measurement
           noiseAmp = zeros(size(amp));
           noisePhase = zeros(size(phase));
           %Uncomment for correlated amplitude
%                            if wavelength == 850
%                                thisRandom = randn(1);
%                            end
           %Add noise to each frequency
           for thisFreq = 1:length(rp)
               %Uncomment for correlated amplitude
%                                if wavelength == 850
%                                    noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*thisRandom*ampfit(thisFreq,waveIdx);
%                                else
%                                noiseAmp(thisFreq) = amp(thisFreq) + amp(thisFreq).*randn(1)*noiseLevel;%ampfit(thisFreq,waveIdx);
%                                end
               %Add noise to amp/phase
               noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
               noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
           end

           runSumAmp = runSumAmp + noiseAmp;
           runSumPhase = runSumPhase + noisePhase;
        end
        %Average amp/phase
        avgAmp = runSumAmp/numReps;
        avgPhase = runSumPhase/numReps;
        %Convert back to real/imaginary
        noiseReal = avgAmp .* cos(avgPhase);
        noiseImag = avgAmp .* sin(avgPhase);
        recOP = fitMu_newOpts([noiseReal;noiseImag],wt, fa', nind, sep,'test');
        %Save recovered OP pair
        recMua(i) = recOP.mua;
        recMus(i) = recOP.mus;

    end
    %Save values in array
    saveMuaTrue(:,f,1) = trueMua(1:numTrials);
    saveMuaRec(:,f,1) = recMua;
    saveMusTrue(:,f,1) = trueMus(1:numTrials);
    saveMusRec(:,f,1) = recMus;

end 

%Save MAT FILES

save(sprintf('../generatedData/NA_%dnm_optimalFreq.mat',wavelength),'saveMuaTrue','saveMusTrue', 'saveMuaRec','saveMusRec','freqList');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%850 low mua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%850 high mua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../generatedData/digital_850nm_optimalFreq_30mm.mat')
wavelength = 850;
eucErrordDOS = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_850nm_optimalFreq.mat')
eucErrormDOS = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

distMua = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2);
distMus = sqrt(((saveMusTrue-saveMusRec)./saveMusTrue).^2);
avgErrDOS = mean(eucErrordDOS,1);
avgErrMDOS = mean(eucErrormDOS,1);
medErrDOS = median(eucErrordDOS,1);
medErrMDOS = median(eucErrormDOS,1);
stdErr = std(eucErrordDOS,[],1);
stdErrM = std(eucErrormDOS,[],1);
freqs = 50:20:490;

figure
e=plot(freqs,log10(squeeze(medErrDOS(1,1:length(freqs),1))),'o');
c=get(e,'color');
set(e,'markerfacecolor',c)
hold on
f=plot(freqs,log10(squeeze(medErrMDOS(1,1:length(freqs),1))),'o');
c2 = get(f,'color');
set(f,'markerfacecolor',c2);
xlabel('Modulation Frequency (MHz)')
ylabel('log_{10}(Median Error)')
title('Comparison of single frequencies')
xlim([40,425])
legend('all-digital','network analyzer')
set(gca,'fontsize',18)
print('../plots/optimalFrequencies.png','-dpng')


f=figure;
e=errorbar(freqs,log10(avgErrDOS),log10(stdErr),'o-','linewidth',1.5);
c=get(e,'color');
set(e,'markerfacecolor',c)
hold on
f=errorbar(freqs+5,log10(avgErrMDOS),log10(stdErrM),'o-','linewidth',1.5);
c2 = get(f,'color');
set(f,'markerfacecolor',c2);
xlabel('Modulation Frequency (MHz)')
ylabel('log_{10}(Median Error)')
title('Comparison of single frequencies')
xlim([40,410])
legend('all-digital','network analyzer')
set(gca,'fontsize',18)
print('../plots/optimalFrequencies_errBar.png','-dpng')

f=figure;
e=boxplot(log10(eucErrordDOS),'symbol','');
hold on
f=boxplot(log10(eucErrormDOS),'symbol','');
xlabel('Modulation Frequency (MHz)')
ylabel('log_{10}(Error)')
title('Comparison of single frequencies')
set(gca,'fontsize',18)
print('../plots/optimalFrequencies_boxPlot.png','-dpng')

