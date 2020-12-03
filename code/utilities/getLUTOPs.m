%getLUTOPs.m
%
%PURPOSE: Get the Optical properties estimated using the LUT method
%
%INPUTS:
%  LUT:      Structure of the lookup table
%  ydata:    Raw input data
%  plotFlag: Do you want to plot the chi^2 space
%
%OUTPUTS:
%  recMua: Estimated Mua
%  recMus: Estimated Scattering
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [recMua,recMus] = getLUTOPs(LUT,ydata, plotFlag)
    %Extract real and imaginary part of the input data
    datReal(1,1,:) = ydata(1,:);
    datImag(1,1,:) = ydata(2,:);
    if size(datReal,3) ~= size(LUT.rp,3)
        error('Input data is the wrong shape. Check your transposes')
    end
    %Calculate the chi^2 statistic
    chisq = sum((datReal - LUT.rp).^2 + (datImag - LUT.ip).^2,3);
    %Find the index of the minimum chi^2
    [~,minIdx] = min(chisq(:));
    %Return the recovered OPs
    recMua = LUT.muaGrid(minIdx);
    recMus = LUT.musGrid(minIdx);
    %If you want to plot, make a plot
    if plotFlag == 1
        figure
        imagesc(LUT.muaGrid(1,:),LUT.musGrid(:,1),log10(chisq))
        pause
    end
end

