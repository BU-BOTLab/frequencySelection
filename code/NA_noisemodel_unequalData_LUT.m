%NA_noisemodel_unequalData_LUT
%
%PURPOSE: Test the default inversion algorithm with noise 
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
waveIdx = 5;
%Allocate memory for MAT file variables
saveMuaTrue = zeros(numTrials,length(freqList),1);
saveMusTrue = zeros(numTrials,length(freqList),1);
saveMuaRec = zeros(numTrials,length(freqList),1);
saveMusRec = zeros(numTrials,length(freqList),1);

for f=1:length(freqList)
	fprintf("\tworking on num freqs %d of %d\n", f,length(freqList));              
	fa = freqList{f};
    %Generate the LUT
	LUT = generateLUT(allMua,allMus,sep,fa);
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
	
        noiseAmp = zeros(size(amp));
        noisePhase = zeros(size(phase));
        for thisFreq = 1:length(amp)
             noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*(ampeqn(thisFreq,waveIdx,1)*randMua+ampeqn(thisFreq,waveIdx,2));%ampfit(thisFreq,waveIdx);
             noisePhase(thisFreq) = phase(thisFreq) + randn(1)*(phaseeqn(thisFreq,waveIdx,1)*randMua+phaseeqn(thisFreq,waveIdx,2));%phasefit(thisFreq,waveIdx);          
        end

		%Convert back to real/imaginary
		noiseReal = noiseAmp .* cos(noisePhase);
		noiseImag = noiseAmp .* sin(noisePhase);
		%Run inverse model
		[recMuaLUT(i),recMusLUT(i)] = getLUTOPs(LUT,[noiseReal,noiseImag]',0);
	end
	%Save values in array
	saveMuaTrue(:,f,1) = trueMua(1:numTrials);
	saveMuaRec(:,f,1) = recMuaLUT;
	saveMusTrue(:,f,1) = trueMus(1:numTrials);
	saveMusRec(:,f,1) = recMusLUT;

end 
%%
%Save MAT FILES
save(sprintf('../generatedData/NA_noise_%dnm_unequalData_LUT.mat',wavelength),'saveMuaTrue','saveMusTrue', 'saveMuaRec','saveMusRec','freqList');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%850
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../generatedData/NA_noise_850nm_unequalData_LUT.mat')
%wavelength = 850;
eucError = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

strFreqList = {'70','110','500','50, 500', '50:7:253','50:499'};
sep = 28;
%%%Euclidean error
f=figure;
set(f,'Position',[230,260,1470,600])
sgtitle('mDOSI noise model unequal data LUT 850 nm')
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
print('../plots/ErrorLandscape_NA_850nm_unequalData_LUT_v2.png','-dpng')

