function [LUT] = generateLUT(allMua,allMus, rho, freqs)
%Generate LUTs, returns structure of look up tables with values at allMua
%and allMus
%INPUTS:
%   allMua -- list of all absorption values to calculate
%   allMus -- list of all scattering values to calculate
%   rhos -- list of all source/detector separations to calculate
%   freqList -- list of modulation frequencies to calculate
%OUTPUTS:
%   LUT  -- Structure of LUT with following fields
        %rp -- real part of the data (size length(allMua)xlength(allMus)xnumFreqs)
        %ip -- imaginary part of the data (same size as above)
        %rho -- source-detector separation
        %freqs -- modulation frequencies
        %muaGrid -- Grid of mua values in the LUT
        %musGrid -- Grid of mus values in the LUT
%Also saves a copy of the LUT in the LUTs folder (if it's not there)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('dosiTrainingCode2019'))
wt = 0; %don't weight forward data
reImFlag = 1; %calculate real/imaginary part
nind = 1.4; %IOR
numEls = length(allMua); %Number of muas (for saving purposes)
if length(freqs) == 1
    saveName = sprintf('DOSILUT_%dmm_%dMHz_n%d.mat',rho,freqs,numEls);
elseif length(freqs) == 2
    saveName = sprintf('DOSILUT_%dmm_%dMHz_and_%dMHz_n%d.mat',rho,min(freqs),max(freqs),numEls);
else
    saveName = sprintf('DOSILUT_%dmm_%dMHz_to_%dMHz_by_%dMHz_n%d.mat',rho,min(freqs),max(freqs),round(mean(diff(freqs))),numEls);
end

%If the LUT already exists don't bother to regenerate it
if exist(fullfile('../LUTs',saveName),'file')
    load(fullfile('../LUTs',saveName),'LUT');
    disp('LUT already exists (move/rename it to regenerate)')
%Otherwise, make the LUT
else
    disp('Generating new LUT')
    %Turn vector into grid
    [X,Y ]= meshgrid(allMua, allMus);
    %Separate the real and imaginary parts
    rp = zeros(size(X,1),size(X,2),length(freqs));
    ip = zeros(size(X,1),size(X,2),length(freqs));
    %Generate the forward data and save it
    for i = 1:size(X,1)
        if mod(i,10) == 0
            fprintf('Finished %d of %d\n',i,size(X,1));
        end
        parfor j = 1:size(X,2) %Use a parallel for loop for speed
            rawDat = p1seminf_mba([X(i,j),Y(i,j)],freqs,nind,rho,wt,reImFlag);
            realPart = rawDat(1,:);
            imagPart = rawDat(2,:);
            rp(i,j,:) = realPart;
            ip(i,j,:) = imagPart;
        end
    end
    %Save LUT the structure and the MAT file
    LUT.rp = rp;
    LUT.ip = ip;
    LUT.freqs = freqs;
    LUT.rho = rho;
    LUT.muaGrid = X;
    LUT.musGrid = Y;
    save(fullfile('../LUTs',saveName),'LUT','-v7.3');
%end
end

