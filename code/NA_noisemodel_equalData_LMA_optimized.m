%NA_noisemodel_equalData_LMA_optimized
%
%PURPOSE: Test the OPTIMIZED inversion algorithm with noise 
%added according to the network analyzer noise model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Set random seed for reproducibility
rng(14850);
%Add forward and inverse model to path
addpath(genpath('dosiTrainingCode2019'))
addpath(genpath('utilities'))
load('../generatedData/OPSet10k.mat')
%Data directories containing raw data for building noise model
dataDir ='../experimentalData/networkAnalyzerSystem/SystemNoiseOP_Review/201027';
%Samples to use for noise models
%Samples to use for noise models
sampNames =  {'bpav4*.asc','b4h2*.asc'};
sampMuas = [0.003,0.02];
%Acquisition parameters for noise model samples
numDiodes = 5;
numMeasFreqs = 401;
nind = 1.4;

%Parameters to use for simulations
sep = 28;
wt =0;
reImFlag = 1;
numTrials = length(trueMua);
smoothing = 15;
%Range of OPs
muaAll = [0.001,0.05];
musAll = [.1,3];

%Allocate memory for OPs
recMua = zeros(1,numTrials);
recMus = zeros(1,numTrials);

%Frequencies to explore
freqList = {[70],[110],[500],[50,500],[50:7:253], [50:499]};
%Number of frequencies in each test set
numFreqs = [1,1,1,2,30,450];
%Maximum number of frequencies
maxFreqs = max(numFreqs);

%Iterate over all separation ranges
wavelength = 850;
waveIdx = 5;
%Allocate memory for MAT file variables
saveMuaTrue = zeros(numTrials,length(freqList),1);
saveMusTrue = zeros(numTrials,length(freqList),1);
saveMuaRec = zeros(numTrials,length(freqList),1);
saveMusRec = zeros(numTrials,length(freqList),1);

for f=1:length(freqList)
	fprintf("\tworking on num freqs %d of %d\n", f,length(freqList));              
	fa = freqList{f};
	%Calculate noise model
	[asd, psd, ampeqn,phaseeqn,f2]=getAmpPhaseSD_v2(dataDir,sampNames,sampMuas,numMeasFreqs,numDiodes,fa, smoothing);
	%Pick a bunch of OP pairs to test
	parfor i=1:numTrials

		%This OP pair
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
		recOP = fitMu_newOpts([noiseReal;noiseImag],wt, fa', nind, sep,'test');
		%Save recovered OP pair
		recMua(i) = recOP.mua;
		recMus(i) = recOP.mus;
	end
	%Save values in array
	saveMuaTrue(:,f,1) = trueMua;
	saveMuaRec(:,f,1) = recMua;
	saveMusTrue(:,f,1) = trueMus;
	saveMusRec(:,f,1) = recMus;

end 
%%
%Save MAT FILES
save(sprintf('../generatedData/NA_noise_%dnm_equalData_optimized.mat',wavelength),'saveMuaTrue','saveMusTrue', 'saveMuaRec','saveMusRec','freqList');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%850
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../generatedData/NA_noise_850nm_equalData_optimized.mat')
%wavelength = 850;
eucError = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

strFreqList = {'70','110','500','50, 500', '50:7:253','50:499'};
sep = 28;
%%%Euclidean error
f=figure;
set(f,'Position',[230,260,1470,600])
sgtitle('mDOSI noise model equal data optimized 850 nm')
it = 1;
for thisFreqIdx = 1:length(strFreqList)
    subplot(2,3,it)
    scatter3(saveMuaTrue(:,thisFreqIdx,1),saveMusTrue(:,thisFreqIdx,1),eucError(:,thisFreqIdx,1),4,log10(eucError(:,thisFreqIdx,1)),'filled')
    xlabel('\mu_a (1/mm)')
    ylabel('\mu_s'' (1/mm)')
    c=colorbar;
    caxis([-15,2])
    ylabel(c,'log_{10}(Error)')
    title(sprintf('%d mm, %s MHz',sep, strFreqList{thisFreqIdx}))
    view(0,90)
    it = it + 1;
  
end
print('../plots/ErrorLandscape_NA_850nm_equalData_optimized_v2.png','-dpng')

