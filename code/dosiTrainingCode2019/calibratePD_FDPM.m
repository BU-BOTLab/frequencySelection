function calibrated=calibratePD_FDPM(cal,raw)
% send in calibration factors and raw data and return calibrated data
% Output:
% calibrated.AC
% calibrated.phase
% calibrated.damp
% calibrated.dphi
% calibrated.freq
% calibrated.dist


calibrated.AC=raw.AC./cal.AC;
calibrated.phase=raw.phase-cal.phase;
calibrated.real = calibrated.AC .* cos(calibrated.phase);
calibrated.imag = calibrated.AC .* sin(calibrated.phase);

calibrated.damp = calibrated.AC.*sqrt((raw.ACsd./raw.AC).^2 + cal.ACsd_AC_sqd);
calibrated.dphi = sqrt(cal.phsd_sqd + raw.phsd.^2);
% calibrated.dreal = calibrated.real.*sqrt((raw.realsd./raw.real).^2 + cal.real_sd_sqd);
% calibrated.dimag = calibrated.imag.*sqrt((raw.imagsd./raw.imag).^2 + cal.imag_sd_sqd);

calibrated.freq = raw.freq;

calibrated.dist = raw.dist;


