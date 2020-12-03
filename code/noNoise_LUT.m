%noNoise_LMA_default.m
%
%PURPOSE: Test the accuracy of the default inversion algorithm on data without any
%added noise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rng(14850) %Set random seed for reproducibility
%Add forward and inverse model to path
addpath(genpath('dosiTrainingCode2019'))
addpath(genpath('utilities'))
%Load the OP set
load('../generatedData/OPSet10k.mat')


%Parameters to use for simulations
wt =0;
reImFlag = 1;
nind = 1.4; %Index of refraction
numTrials = length(trueMua);

%Range of OPs
muaAll = [0.001,0.05];
musAll = [.1,3];
%Grid of OPs to calculate LUT
allMua = linspace(muaAll(1),muaAll(2),500);
allMus = linspace(musAll(1),musAll(2),500);

%Allocate memory for OPs
recMuaLUT = zeros(1,numTrials);
recMusLUT = zeros(1,numTrials);

sepRange= [10,20,30]; %Range of source/detector separations
%Frequencies to explore
freqList = {[70],[110],[500],[50,500],[50:7:253], [50:499]};

%Allocate memory for MAT file variables
saveMuaTrue = zeros(numTrials,length(freqList),length(sepRange));
saveMusTrue = zeros(numTrials,length(freqList),length(sepRange));
saveMuaRec = zeros(numTrials,length(freqList),length(sepRange));
saveMusRec = zeros(numTrials,length(freqList),length(sepRange));

%Iterate over all separation ranges
for sepIdx = 1:length(sepRange)
    fprintf('working on SDsep %d of %d\n', sepIdx,length(sepRange));
    for f=1:length(freqList)
        fprintf("\tworking on num freqs %d of %d\n", f,length(freqList));              
        fa = freqList{f};
        LUT = generateLUT(allMua,allMus,sepRange(sepIdx),fa);
        parfor i=1:numTrials
            %Get this optical property pair
            randMua = trueMua(i);
            randMus = trueMus(i);
            %Generate forward model
            rawDat = p1seminf_mba([randMua,randMus],fa',nind,sepRange(sepIdx),wt,reImFlag);
            %extract real and imaginary parts
            rp = rawDat(1:length(fa));
            ip = rawDat(length(fa)+1:end);
            %Run inverse model
            [recMuaLUT(i),recMusLUT(i)] = getLUTOPs(LUT,[rp,ip]',0);
        end
        %Save values in array
        saveMuaTrue(:,f,sepIdx) = trueMua;
        saveMuaRec(:,f,sepIdx) = recMuaLUT;
        saveMusTrue(:,f,sepIdx) = trueMus;
        saveMusRec(:,f,sepIdx) = recMusLUT;

    end 
end
%%
%Save MAT FILES
save('../generatedData/noNoise_LUT.mat', 'saveMuaTrue','saveMusTrue','saveMuaRec','saveMusRec','freqList');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../generatedData/noNoise_LUT.mat')
eucError = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
strFreqList = {'70','110','500','50, 500', '50:7:253','50:499'};

sepRange = [10,20,30];
%%%Euclidean error
f=figure;
set(f,'Position',[230,260,1470,900])
sgtitle('No noise error single freq LUT')
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
print('../plots/ErrorLandscapeSingleFreq_noNoise_LUT.png','-dpng')


f=figure;
set(f,'Position',[230,260,1470,900])
sgtitle('No noise error multi-freq default')
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
print('../plots/ErrorLandscapeMultiFreq_noNoise_LUT.png','-dpng')


