% %%
% %Get 1 frequency no noise accuracy vs LUT density with linear scale
% 
% numPts = logspace(2,4,10);
% addpath(genpath('utilities')) 
% load('../generatedData/OPSet10k.mat')
% 
% freq = 110;
% rho = 30;
% nind = 1.4;
% numOPs = 10000;
% muaErr = zeros(1,length(numPts));
% muaStd = muaErr;
% musErr = zeros(1,length(numPts));
% musStd = musErr;
% timingLin = zeros(1,length(numPts));
% muaDistLin = zeros(numOPs,length(numPts));
% musDistLin = zeros(numOPs,length(numPts));
% for i = 1:length(numPts)
%     fprintf("Working on density %d of %d\n", i, length(numPts))
%     muaAll = linspace(0.0001,0.05,numPts(i));
%     musAll = linspace(0.1,3,numPts(i));
%     
%     LUT = generateLUT(muaAll,musAll,rho,freq);
%     recMua = zeros(1,numOPs);
%     recMus = zeros(1,numOPs);
%     tic
%     for j = 1:numOPs
%         rawDat = p1seminf([trueMua(j),trueMus(j)],freq,0,nind,rho,0,0,1,1);
%         
%         [recMua(j),recMus(j)] = getLUTOPs(LUT,rawDat,0);
%     end
%     timingLin(i) =toc;
%     muaDistLin(:,i) = sqrt((trueMua-recMua).^2);
%     musDistLin(:,i) = sqrt((trueMus-recMus).^2);
%     
%     muaErr(i) = mean(muaDistLin(:,i));
%     musErr(i) = mean(musDistLin(:,i));
%     muaStd(i) = std(muaDistLin(:,i));
%     musStd(i) = std(musDistLin(:,i));
%     
% end
% eucError = sqrt(muaDistLin.^2 + musDistLin.^2);
% %%
% %Plot the linspace plots
% figure
% loglog(numPts,numOPs./timingLin,'o-')
% xlabel('LUT Density')
% ylabel('Inversions/s')
% title('LUT speed')
% print('../plots/LUT_speed_vs_density_noNoise.png','-dpng')
% 
% figure
% [h,L,MX,MED]=violin(log10(eucError));
% %mxMuaLin = MX;
% xticks([1:length(numPts)])
% xticklabels(numPts);
% set(gcf,'Position', [600,800,900,550]);
% xtickangle(45);
% xlabel('LUT density')
% ylabel('log_{10}(Error)')
% title('Noise-free LUT accuracy')
% print('../plots/LUT_error_vs_density_noNoise.png','-dpng')

%%
%From the previous section I have learned that a 2000 x 2000 LUT is still
%3-4x faster than the default method, so let's use that as our HD LUT
%Start with no-noise
freq = 110;
rho = 30;
nind = 1.4;
numOPs = 10000;
numPts = 2000;

fprintf("Working on noise-free noise model\n")
muaAll = linspace(0.0001,0.05,numPts);
musAll = linspace(0.1,3,numPts);
LUT = generateLUT(muaAll,musAll,rho,freq);
recMua = zeros(1,numOPs);
recMus = zeros(1,numOPs);
tic
for j = 1:numOPs
    rawDat = p1seminf_mba([trueMua(j),trueMus(j)],freq,nind,rho,wt,reImFlag);

    [recMua(j),recMus(j)] = getLUTOPs(LUT,rawDat,0);
end
timingNoNoise =toc;
    
save('../generatedData/noNoise_110MHz_HDLUT_2k.mat','trueMua','trueMus','recMua','recMus','timingNoNoise')

%%
%From the previous section I have learned that a 2000 x 2000 LUT is still
%3-4x faster than the default method, so let's use that as our HD LUT
%Start with no-noise
freq = 500;
rho = 30;
nind = 1.4;
numOPs = 10000;
numPts = 2000;

fprintf("Working on noise-free noise model\n")
muaAll = linspace(0.0001,0.05,numPts);
musAll = linspace(0.1,3,numPts);
LUT = generateLUT(muaAll,musAll,rho,freq);
recMua = zeros(1,numOPs);
recMus = zeros(1,numOPs);
tic
for j = 1:numOPs
    rawDat = p1seminf_mba([trueMua(j),trueMus(j)],freq,nind,rho,wt,reImFlag);

    [recMua(j),recMus(j)] = getLUTOPs(LUT,rawDat,0);
end
timingNoNoise =toc;
    
save('../generatedData/noNoise_500MHz_HDLUT_2k.mat','trueMua','trueMus','recMua','recMus','timingNoNoise')

%%
%Repeat this adding noise similar to the all digital system
addpath(genpath('utilities')) 
load('../generatedData/OPSet10k.mat')

%Data directories containing raw data for building noise model
dataDir30 = '../experimentalData/allDigitalSystem/Gen2_30mm/8192/200614';
%Samples to use for noise models
sampName =  'bpav4*.asc';
waveIdx = 6;
%Acquisition parameters for noise model samples
numDiodes = 6;
numMeasFreqs = 70;

freq = 110;
rho = 30;
nind = 1.4;
numOPs = 10000;
numReps = 450;
numPts = 2000;

fprintf("Working on digital noise model\n")
muaAll = linspace(0.0001,0.05,numPts);
musAll = linspace(0.1,3,numPts);

LUT = generateLUT(muaAll,musAll,rho,freq);
%Calculate noise model
[~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir30,sampName,numMeasFreqs,numDiodes,freq, 5);
recMua = zeros(1,numOPs);
recMus = zeros(1,numOPs);
tic
for j = 1:numOPs
    rawDat = p1seminf_mba([trueMua(j),trueMus(j)],freq,nind,rho,wt,reImFlag);
    %extract real and imaginary parts
    rp = rawDat(1:length(freq));
    ip = rawDat(length(freq)+1:end);
    %Convert to amp/phase
    amp = sqrt(rp.^2 + ip.^2);
    phase = unwrap(atan2(ip,rp));
    noiseAmp = zeros(size(amp));
    noisePhase = zeros(size(phase));
    %Add noise to each frequency
    runSumAmp = zeros(size(amp));
    runSumPhase = zeros(size(phase));
    for r = 1:numReps
        for thisFreq = 1:size(amp)
            noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
            noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
        end
       runSumAmp = runSumAmp + noiseAmp;
       runSumPhase = runSumPhase + noisePhase;
    end
    avgAmp = runSumAmp /numReps;
    avgPhase = runSumPhase / numReps;
    %Convert back to real/imaginary
    noiseReal = avgAmp .* cos(avgPhase);
    noiseImag = avgAmp .* sin(avgPhase);
    %Run inverse model
    [recMua(j),recMus(j)] = getLUTOPs(LUT,[noiseReal;noiseImag],0);
end
timingDigNoise =toc;

save('../generatedData/digitalNoise_110MHz_HDLUT_2k.mat','trueMua','trueMus','recMua','recMus','timingDigNoise')

%%
%Repeat this adding noise similar to the network analyzer system
addpath(genpath('utilities')) 
load('../generatedData/OPSet10k.mat')
%Data directories containing raw data for building noise model
dataDir = '../experimentalData/networkAnalyzerSystem/mDOSI_Probe2020_Drift/200304';
%Samples to use for noise models
sampName =  'INO9*miniLBS.asc';
%Acquisition parameters for noise model samples
numDiodes = 6;
numMeasFreqs = 401;

freq = 110;
rho = 28;
nind = 1.4;
numOPs = 10000;
numReps = 450;
numPts = 2000;

fprintf("Working on network analyzer noise model\n")
muaAll = linspace(0.0001,0.05,numPts);
musAll = linspace(0.1,3,numPts);

LUT = generateLUT(muaAll,musAll,rho,freq);
%Calculate noise model
[~, ~, ampfit,phasefit,~]=getAmpPhaseSD(dataDir,sampName,numMeasFreqs,numDiodes,freq, 15);
recMua = zeros(1,numOPs);
recMus = zeros(1,numOPs);
tic
for j = 1:numOPs
    rawDat = p1seminf_mba([trueMua(j),trueMus(j)],freq,nind,rho,wt,reImFlag);
    %extract real and imaginary parts
    rp = rawDat(1:length(freq));
    ip = rawDat(length(freq)+1:end);
    %Convert to amp/phase
    amp = sqrt(rp.^2 + ip.^2);
    phase = unwrap(atan2(ip,rp));
    noiseAmp = zeros(size(amp));
    noisePhase = zeros(size(phase));
    %Add noise to each frequency
    runSumAmp = zeros(size(amp));
    runSumPhase = zeros(size(phase));
    for r = 1:numReps
        for thisFreq = 1:size(amp)
            noiseAmp(thisFreq) = amp(thisFreq)+amp(thisFreq)*randn(1)*ampfit(thisFreq,waveIdx);
            noisePhase(thisFreq) = phase(thisFreq) + randn(1)*phasefit(thisFreq,waveIdx);        
        end
       runSumAmp = runSumAmp + noiseAmp;
       runSumPhase = runSumPhase + noisePhase;
    end
    avgAmp = runSumAmp /numReps;
    avgPhase = runSumPhase / numReps;
    %Convert back to real/imaginary
    noiseReal = avgAmp .* cos(avgPhase);
    noiseImag = avgAmp .* sin(avgPhase);
    %Run inverse model
    [recMua(j),recMus(j)] = getLUTOPs(LUT,[noiseReal;noiseImag],0);
end
timingNANoise =toc;

save('../generatedData/NANoise_110MHz_HDLUT_2k.mat','trueMua','trueMus','recMua','recMus','timingDigNoise')
%%
%Plot the linspace plots
load('../generatedData/noNoise_110MHz_HDLUT_2k.mat')
eucErrorNoNoise = sqrt(((trueMua-recMua)./trueMua).^2 + (((trueMus-recMus)./trueMus)).^2);
load('../generatedData/digitalNoise_110MHz_HDLUT_2k.mat')
eucErrorDigNoise = sqrt(((trueMua-recMua)./trueMua).^2 + (((trueMus-recMus)./trueMus)).^2);
load('../generatedData/NANoise_110MHz_HDLUT_2k.mat')
eucErrorNANoise = sqrt(((trueMua-recMua)./trueMua).^2 + (((trueMus-recMus)./trueMus)).^2);

figure
histogram(log10(eucErrorNoNoise))





