%Digital_noisemodel_equalData_LUT
%
%PURPOSE: Test the LUT inversion algorithm with noise 
%added according to the network analyzer noise model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Set random seed for reproducibility
close all
clear all

rng(14850);
%Add forward and inverse model to path
addpath(genpath('dosiTrainingCode2019'))
addpath(genpath('utilities'))
load('../generatedData/OPSet10k.mat')
%Data directories containing raw data for building noise model
dataDir10 = '../experimentalData/allDigitalSystem/Gen2_10mm/8192/200614';
dataDir20 = '../experimentalData/allDigitalSystem/Gen2_20mm/8192/200614';
dataDir30 = '../experimentalData/allDigitalSystem/Gen2_30mm/8192/200614';
%Samples to use for noise models
sampNames =  {'bpav4*.asc','b4h2*.asc'};
sampMuas = [0.003,0.02];
%Acquisition parameters for noise model samples
numDiodes = 6;
numMeasFreqs = 70;
nind = 1.4;

%Parameters to use for simulations
sep = [30];
wt =0;
reImFlag = 1;
numTrials = 10000;
smoothing = 5;
%Range of OPs
muaAll = [0.001,0.05];
musAll = [.1,3];
%Grid of OPs to calculate LUT
allMua = linspace(muaAll(1),muaAll(2),500);
allMus = linspace(musAll(1),musAll(2),500);

%Allocate memory for OPs
recMuaLUT = zeros(1,numTrials);
recMusLUT = zeros(1,numTrials);

%Frequencies to explore
freqList = {[70],[110],[500],[50,500],[50:7:253], [50:499]};
%Number of frequencies in each test set
numFreqs = [1,1,1,2,30,450];
%Maximum number of frequencies
maxFreqs = max(numFreqs);

%Iterate over all separation ranges
wavelength = 850;
waveIdx = 6;
%Allocate memory for MAT file variables
saveMuaTrue = zeros(numTrials,length(freqList),length(sep));
saveMusTrue = zeros(numTrials,length(freqList),length(sep));
saveMuaRec = zeros(numTrials,length(freqList),length(sep));
saveMusRec = zeros(numTrials,length(freqList),length(sep));

%Iterate through each SD separation
for sepIdx = 1:length(sep)
    if sep(sepIdx) == 10
        dataDir = dataDir10;
    elseif sep(sepIdx) == 20
        dataDir = dataDir20;
    elseif sep(sepIdx) == 30
        dataDir = dataDir30;
    else
        error('Unknown SD separation')
    end
    %Iterate through each modulation frequency set
    for f=1:length(freqList)
        fprintf("\tworking on num freqs %d of %d\n", f,length(freqList));              
        fa = freqList{f};
        %Generate the LUT
        LUT = generateLUT(allMua,allMus,sep(sepIdx),fa);
        %Calculate noise model
        [asd, psd, ampeqn,phaseeqn,f2]=getAmpPhaseSD_v2(dataDir,sampNames,sampMuas,numMeasFreqs,numDiodes,fa, smoothing);
        %Pick a bunch of OP pairs to test
        if length(fa) < 400
            modIdx = 1000;
        else
            modIdx = 100;
        end
        parfor i=1:numTrials
            if mod(i,modIdx) == 0
                fprintf("On trial %d\n",i)
            end
            %This OP pair
            randMua=trueMua(i);
            randMus=trueMus(i);
            %Generate raw data
            rawDat = p1seminf_mba([randMua,randMus],fa',nind,sep(sepIdx),wt,reImFlag);
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
            for rep=1:floor(maxFreqs/numFreqs(f))
                %Allocate memory for noisy measurement
                noiseAmp = zeros(size(amp));
                noisePhase = zeros(size(phase));
                 %Add noise to each frequency
                for thisFreq = 1:numFreqs(f)
                    noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*(ampeqn(thisFreq,waveIdx,1)*randMua+ampeqn(thisFreq,waveIdx,2));%ampfit(thisFreq,waveIdx);
                    noisePhase(thisFreq) = phase(thisFreq) + randn(1)*(phaseeqn(thisFreq,waveIdx,1)*randMua+phaseeqn(thisFreq,waveIdx,2));%phasefit(thisFreq,waveIdx);        
                end
                runSumAmp = runSumAmp + noiseAmp;
                runSumPhase = runSumPhase + noisePhase;
            end
            %Average amp/phase
            avgAmp = runSumAmp/floor(maxFreqs/numFreqs(f));
            avgPhase = runSumPhase/floor(maxFreqs/numFreqs(f));
            %Convert back to real/imaginary
            noiseReal = avgAmp .* cos(avgPhase);
            noiseImag = avgAmp .* sin(avgPhase);
            %Run inverse model
            [recMuaLUT(i),recMusLUT(i)] = getLUTOPs(LUT,[noiseReal,noiseImag]',0);
        end
        %Save values in array
        saveMuaTrue(:,f,sepIdx) = trueMua(1:numTrials);
        saveMuaRec(:,f,sepIdx) = recMuaLUT;
        saveMusTrue(:,f,sepIdx) = trueMus(1:numTrials);
        saveMusRec(:,f,sepIdx) = recMusLUT;

    end 
end

%%
    %Save MAT FILES
    save(sprintf('../generatedData/Digital_noise_%dnm_equalData_LUT_v2.mat',wavelength),'saveMuaTrue','saveMusTrue', 'saveMuaRec','saveMusRec','freqList');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%850nm digital
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../generatedData/Digital_noise_850nm_equalData_LUT_v2.mat')
wavelength = 850;
eucError = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

strFreqList = {'70','110','500','50, 500', '50:7:253','50:499'};
sepRange = [30];
%%%Euclidean error
f=figure;
set(f,'Position',[230,260,1470,900])
sgtitle('Digital noise model error single freq LUT')
it = 1;
for thisFreqIdx = 1:3
    for thisSepIdx = 1:length(sepRange)
        subplot(3,3,it)
        scatter3(saveMuaTrue(:,thisFreqIdx,thisSepIdx),saveMusTrue(:,thisFreqIdx,thisSepIdx),eucError(:,thisFreqIdx,thisSepIdx),4,log10(eucError(:,thisFreqIdx,thisSepIdx)),'filled')
        xlabel('\mu_a (1/mm)')
        ylabel('\mu_s'' (1/mm)')
        c=colorbar;
        caxis([-15,2])
        ylabel(c,'log_{10}(Error)')
        title(sprintf('%d mm, %s MHz',sepRange(thisSepIdx), strFreqList{thisFreqIdx}))
        view(0,90)
        it = it + 1;
    end
end
print('../plots/ErrorLandscapeSingleFreq_digital_850nm_LUT_v3.png','-dpng')


f=figure;
set(f,'Position',[230,260,1470,900])
sgtitle('Digital noise model error multi-freq LUT')
it = 1;
for thisFreqIdx = 4:6
    for thisSepIdx = 1:length(sepRange)
        subplot(3,3,it)
        scatter3(saveMuaTrue(:,thisFreqIdx,thisSepIdx),saveMusTrue(:,thisFreqIdx,thisSepIdx),eucError(:,thisFreqIdx,thisSepIdx),4,log10(eucError(:,thisFreqIdx,thisSepIdx)),'filled')
        xlabel('\mu_a (1/mm)')
        ylabel('\mu_s'' (1/mm)')
        c=colorbar;
        caxis([-15,2])
        ylabel(c,'log_{10}(Error)')
        title(sprintf('%d mm, %s MHz',sepRange(thisSepIdx), strFreqList{thisFreqIdx}))
        view(0,90)
        it = it + 1;
    end
end
print('../plots/ErrorLandscapeMultiFreq_digital_850nm_LUT_v3.png','-dpng')