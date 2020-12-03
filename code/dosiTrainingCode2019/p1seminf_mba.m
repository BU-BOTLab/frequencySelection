%========================================================================================
%SEMI-INFINITE PHASE AND AMPLITUDE CALCULATOR FOR P1 MODEL, v 2.5
%========================================================================================
% Returns phase (radians) and amplitude of P1 seminfinite Photon density
% waves (PDW).  See Cerussi, A. and B. Tromberg (2003), Photon Migration Spectroscopy: 
% Frequency Domain Techniques. Biomedical Photonics Handbook. T. Vo-Dinh,
% CRC Press: 22.1- 22.17 for more details.
%
% USAGE
%    [Y, VER] = P1SEMINF(P,F,NIND,RHO,WT,REIM_FLAG);  
%
% Note on Units:    Accepts MHz for frequency, 1/mm for optical properties, and mm
%   for distances
%Modified from code developed at the Beckman Laser Institute.

%INPUTS: p -- 1x2 element vector of [mua, musp] in 1/mm
%        f -- frequencies to calculate in MHz. Be careful of shape as long
%             vs. tall vectors will give differently shaped outputs
%        nind -- Index of refraction of the semi-infinite medium
%        rho1 -- source-detector separation in mm
%        wt -- Value to weight the elements by which can be useful for
%        curve fitting. Use 0 for no weighting (1 will also work)
%        reim_flag -- 1 Returns the real/imaginary parts of the response, 0
%                     returns the amplitude and phase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = p1seminf_mba(p,f,nind,rho1, wt, reim_flag)

%*************************************************************************
%PRELIMINARY MATTERS
%*************************************************************************
mua = p(1);% Absorption coefficient in mm-1
mus = p(2);% Reduced scattering in mm-1
%********************************************************************
%BASIC CALCULATIONS
%*************************************************************************
%definitions 
c = 2.99792458e11/nind;				% speed of light in mm/s
mutr = mua+mus;                     %transport coefficient (1/mm)
ltr = 1/mutr;                       %transport length (mm)
D=1/3*ltr;                          %diffusion coefficient, uses the mua for kicks
I=sqrt(-1);                         %Imaginary unit
fbc = 1.0e6*2*pi*f./c;	% Angular frequency in MHZ for omega
alpha=3*fbc*D;		%Josh definition, such that alpha is 2pi*Tcoll/Tmod (page 109).

%boundary conditions
if nind==1.4
    reff=0.493;
elseif nind==1.33
    reff=0.431;
else
    %polynomial fit 6 order by Sophie and AEC
    reff = 2.1037.*nind.^6-19.8048.*nind.^5+76.8786.*nind.^4-156.9634.*nind.^3+176.4549.*nind.^2 -101.6004.*nind+22.9286;
end
%calcuate true distances   
zb = 2/3*(1+reff)/(1-reff)*ltr;					%extrapolated boundary
r01 = sqrt(ltr*ltr+rho1.*rho1);  	            %s-d separation for source, dist 1   
rb1 = sqrt((2*zb+ltr)*(2*zb+ltr)+rho1.*rho1);	%s-d separation for image, dist 1

%===kvector=====================================================================================
k_josh=sqrt((mua-fbc.*alpha-I.*(fbc+mua*alpha))./D);	%complete P1 wavevector
kr=abs(real(k_josh));	
ki=abs(imag(k_josh));		

%==photon density wave===========================================================================
er01 = exp(-kr.*r01);					%exponentials in form e(-kr)/r
erb1 = exp(-kr.*rb1);
Re1 = (+(er01./r01).*cos(ki.*r01) - (erb1./rb1).*cos(ki.*rb1))./D;
Im1 = (+(er01./r01).*sin(ki.*r01) - (erb1./rb1).*sin(ki.*rb1))./D;  

%*************************************************************************
%OUTPUT
%*************************************************************************
if(reim_flag == 0)
    fa=sqrt(Re1.^2+Im1.^2);
    fb=unwrap(atan2(Im1,Re1));		%in radians
else
    fa=Re1;
    fb=Im1;
end
%final outputs
y=[fa; fb];		%final return

if(wt~=0)
 y = y.*wt;  %Weight for the curvefitting
end