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

function [ampMu, phaseMu, ampFitMu, phaseFitMu, f] = getAmpPhaseMu(dataDir, sampName, numFreqs, numDiodes,fitFreqs,smoothing)
    %Set smoothing to 0 if not specified
    if nargin < 6
        smoothing = 0;
    end
    %List of files
    flist = dir(fullfile(dataDir,sampName));
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
    ampMu = mean(amp,3);
    %Same for phase
    %normPhase = (phase)./repmat(mean(phase,3),[1,1,length(flist)]);
    phaseMu = mean(phase,3);
    %Smooth the data
    smoothAmpMu =ampMu;
    smoothPhaseMu = phaseMu;
    if smoothing ~= 0
      for q = 1:size(ampMu,2)
          smoothAmpMu(:,q) = smooth(ampMu(:,q),smoothing);
          smoothPhaseMu(:,q) = smooth(phaseMu(:,q),smoothing);
      end
    end
    
    %Run the interpolation
    ampFitMu = zeros(length(fitFreqs),numDiodes);
    phaseFitMu = zeros(length(fitFreqs),numDiodes);
    
    for t = 1:numDiodes
        phaseFitMu(:,t) = InterpWithClipExtrap(dat.freq,smoothPhaseMu(:,t),fitFreqs);
        ampFitMu(:,t) = InterpWithClipExtrap(dat.freq,smoothAmpMu(:,t),fitFreqs);
    end
    f=dat.freq;
end


