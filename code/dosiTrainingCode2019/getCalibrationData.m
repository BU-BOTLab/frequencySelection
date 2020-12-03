function cal=getCalibrationData(avgCalData, cal_file, n,dist)
%%%byh Calculates instrument response using measurement on phantom of known optical properties 

% calibration files need to be in path, will probably need to specify
% location at some point. . 
% phantom name must not have dash
% OUTPUT:
% cal.error
% cal.dist
% cal.AC
% cal.ACsd_AC_sqd
% cal.phase
% cal.phsd_sqd

    if avgCalData.error~=0
        cal.error=-1;
        return;
    end
    
    cal.dist = dist;
    nFreq=length(avgCalData.freq);         	%find number of data points
    noWt = ones(nFreq*2,1); %No weighting for calibration
    
    %%%byh For phantom files, standardized names are used which match up to
    %%%calibration files in the phantom directory.  the calibration files
    %%%consist of the phantom optical properties at a number of wavelengths.  The
    %%%wavelengths and optical properties are loaded here below 
    cfile=load(cal_file);
    %%%byh The wavelengths in the phantom files do not always match up to the fdpm wavelengths, 
    %%%so we do an interpolation of that dataset to get mua and mus at each
    %%%of the fdpm diode wavelengths
    pmua = interp1(cfile(:,1), cfile(:,2), avgCalData.wavelengths); %mua
    pmus = interp1(cfile(:,1), cfile(:,3), avgCalData.wavelengths); %mus
    
    pdmua = pmua.*0.05; %assume 5% error in phantom OPs
    pdmus = pmus.*0.05; %assume 5% error
    
    
%%%byh Instrument response is calculated for each diode separately.  Using
%%%the known phantom mua and mus at the diode wavelength, a forward calculation is done with the model
%%%and a theoretical phase and ampltidude is calculated covering all the
%%%measured frequencies
    for a=1:avgCalData.nDiodes
        %%%byh Forward calculation - for matlab purposes, the model is set up to return the
        %%%amplitude and phase as a single array, ampltiude at every
        %%%frequency followed by phase at each frequency
        theory=feval('p1seminf',[pmua(a),pmus(a)], avgCalData.freq, 0, n, avgCalData.dist, 0, noWt, 0);
        AC_phan = theory(1:nFreq);
        PHI_phan = theory(1+nFreq:2*nFreq);  %in radians
        real_phan = AC_phan .* cos(PHI_phan);
        imag_phan = AC_phan .* sin(PHI_phan);
        %%%byh Once the theoretical phase and amplitude are calculated,
        %%%this phase is subtracted from the measured phantom phase and the
        %%%measured amplitude is divided by the theoretical amplitude.  The
        %%%result is the instrument response that can then be used to
        %%%calibrate the measurement files
        cal.phase(:,a) =avgCalData.phase(:,a) - PHI_phan;    		
        cal.AC(:,a)   =avgCalData.AC(:,a) ./ AC_phan;
        cal.real(:,a) = avgCalData.real(:,a) ./ real_phan;
        cal.imag(:,a) = avgCalData.imag(:,a) ./ imag_phan;
        
        %mba The weighting is calculated from two sources: the error in
        %calibration phantom OPs and the sensitivity of the amplitude and
        %phase to that error.
        
        %dfdp does a partial deriviative. By making small changes to mua
        %and musp you can calculate how much that would change the
        %amplitude and phase. If the calibration phantom doesn't have
        %exactly the specified phase it can cause problems with the
        %theoretical values. How big the problems will be is based on the
        %derivative
        deriv_phan=dfdp('p1seminf',avgCalData.freq,theory, [pmua(a),pmus(a)],.0001*ones(2,1),0,n,cal.dist,0,noWt,0);
        dACdmua=deriv_phan(1:nFreq,1);		    %derivative of amp with respect to mua
        dPHIdmua=deriv_phan(1+nFreq:2*nFreq,1); %radians derivative of phase wrt mua
        dACdmus=deriv_phan(1:nFreq,2); %Same as above for musp
        dPHIdmus=deriv_phan(1+nFreq:2*nFreq,2); %radians
        %Amplitude and phase weighting based on the error model and the
        %sensitivity
        cal.ACsd_AC_sqd(:,a)=  (avgCalData.ACsd(:,a) ./avgCalData.AC(:,a)).^2 + ((pdmua(a)*dACdmua).^2 +(pdmus(a)*dACdmus).^2)./AC_phan.^2;
        cal.phsd_sqd(:,a)= avgCalData.phsd(:,a).^2 + (pdmua(a)*dPHIdmua).^2 + (pdmus(a)*dPHIdmus).^2 ;
        
        %Same as above for real and imaginary parts (not currently used)
%         deriv_phan_reim=dfdp('p1seminf',avgCalData.freq,[real_phan,imag_phan], [pmua(a),pmus(a)],.0001*ones(2,1),0,n,cal.dist,0,noWt,0);
%         drealdmua=deriv_phan_reim(1:nFreq,1);				%derivative of Re with respect to mua
%         dimagdmua=deriv_phan_reim(1+nFreq:2*nFreq,1); %radians
%         drealdmus=deriv_phan_reim(1:nFreq,2);
%         dimagdmus=deriv_phan_reim(1+nFreq:2*nFreq,2); %radians
%         cal.real_sd_sqd(:,a)=  (avgCalData.real(:,a) ./avgCalData.real(:,a)).^2 + ((pdmua(a)*drealdmua).^2 +(pdmus(a)*drealdmus).^2)./real_phan.^2;
%         cal.imag_sd_sqd(:,a) = (avgCalData.imagsd(:,a) ./ avgCalData.imag(:,a)).^2 + ((pdmua(a)*dimagdmua).^2 + (pdmus(a)*dimagdmus).^2)./imag_phan.^2;
    end
    
cal.error=0;