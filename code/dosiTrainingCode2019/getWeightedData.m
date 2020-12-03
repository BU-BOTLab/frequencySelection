function [YDATA,WT_REIM1] = getWeightedData(cal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Transform from amplitude and phase to real and imaginary parts
realAtEnd = cal.AC .* cos(cal.phase);
imagAtEnd = cal.AC .* sin(cal.phase);
%Concatenate the data
DATA = [realAtEnd;imagAtEnd];
%Calculate the weighting function
ERR_REAL1=sqrt((cal.damp.*cos(cal.phase)).^2 + (cal.dphi.*cal.AC.*sin(cal.phase)).^2);
ERR_IMAG1=sqrt((cal.damp.*sin(cal.phase)).^2 + (cal.dphi.*cal.AC.*cos(cal.phase)).^2);
WT_REIM1 = abs([1./ERR_REAL1;1./ERR_IMAG1]);

%ERR_AMP = abs(ERR_REAL1 + sqrt(-1) * ERR_IMAG1);
%ERR_PHASE = angle(ERR_REAL1+ sqrt(-1)*ERR_IMAG1);

% WT_PAA = abs([1./ERR_AMP(1:51);1./ERR_PHASE(1:51)]);
%WT_PAA = ones(size(DATA));
%Replace any NaNs in the weighting function with 1s
if sum(isnan(WT_REIM1(:))) ~= 0
    warning('NaN in the weighting function replaced with 1')
    idxs = find(isnan(WT_REIM1(:))==1);
    WT_REIM1(idxs) = 1;
end
%Apply the weight
YDATA = DATA.*WT_REIM1;

end

