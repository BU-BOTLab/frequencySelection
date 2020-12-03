%getAmpPhaseSD.m
%
%PURPOSE: Return the amplitude and phase standard deviation based on the
%raw data found in a dataDir with the specified sampName
%
%INPUTS:
%  dataDir:   Directory with the raw data
%  sampName:  Name of the sample to calculate the Std for
%  numFreqs:  Number of modulation frequencies in the actual measurement
%  numDiodes: Number of diodes in the actual measurement
%  fitFreqs:  Frequencies to calculate the std deviation
%  smoothing: Number of data points to smooth over (default=0)
%
%OUTPUTS:
%  ampSD:     Standard deviation of amplitude at the measured freqs
%  phaseSD:   Same as above but for phase
%  ampFit:    Standard deviation of amplitude at fitFreqs
%  phaseFit:  Same as above for phase
%  f:         Actually measured frequencies

function [ampSD, phaseSD, ampEqns, phaseEqns, f] = getAmpPhaseSD_v2(dataDir, sampNames,sampMuas, numFreqs, numDiodes,fitFreqs,smoothing)
    %Set smoothing to 0 if not specified
    if nargin < 6
        smoothing = 0;
    end
     %Run the interpolation
    ampFit = zeros(length(fitFreqs),numDiodes,length(sampNames));
    phaseFit = zeros(length(fitFreqs),numDiodes,length(sampNames));
    ampSD = zeros(numFreqs,numDiodes,length(sampNames));
    phaseSD = zeros(numFreqs,numDiodes,length(sampNames));
    ampEqns = zeros(length(fitFreqs),numDiodes,2);
    phaseEqns = zeros(length(fitFreqs),numDiodes,2);
    for p = 1:length(sampNames)
        %List of files
        flist = dir(fullfile(dataDir,sampNames{p}));
        fnames = cell(1,length(flist));
        if isempty(fnames)
            error('No files on the path with the specified sample name')
        end
        %Preallocate for amplitude and phases
        amp = zeros(numFreqs,numDiodes,length(flist));
        phase = zeros(numFreqs,numDiodes,length(flist));
        %Fill the amplitude and phase vectors
        for i = 1:length(flist)
            fnames{i} = fullfile(dataDir,flist(i).name);
            dat = getASCData(fnames{i});
            amp(:,:,i) = dat.AC;
            phase(:,:,i) = deg2rad(dat.phase);
        end
        %Calculate the normalized amplitude and phase standard deviation
        normAmp = (amp./repmat(mean(amp,3),[1,1,length(flist)]));
        ampSD(:,:,p) = std(normAmp,0,3);
        %Same for phase
        %normPhase = (phase)./repmat(mean(phase,3),[1,1,length(flist)]);
        phaseSD(:,:,p) = std(phase,0,3);
        %Smooth the data
        smoothAmpSD =ampSD;
        smoothPhaseSD = phaseSD;
        if smoothing ~= 0
          for q = 1:size(ampSD,2)
              smoothAmpSD(:,q) = smooth(ampSD(:,q,p),smoothing);
              smoothPhaseSD(:,q) = smooth(phaseSD(:,q,p),smoothing);
          end
        end

        for t = 1:numDiodes
            phaseFit(:,t,p) = InterpWithClipExtrap(dat.freq,smoothPhaseSD(:,t),fitFreqs);
            ampFit(:,t,p) = InterpWithClipExtrap(dat.freq,smoothAmpSD(:,t),fitFreqs);
        end
    end
    f=dat.freq;
    for k = 1:length(fitFreqs)
        for j = 1:numDiodes
            ampEqns(k,j,:) = polyfit(sampMuas',squeeze(ampFit(k,j,:)),1);
            phaseEqns(k,j,:) = polyfit(sampMuas',squeeze(phaseFit(k,j,:)),1);
        end
    end
end


