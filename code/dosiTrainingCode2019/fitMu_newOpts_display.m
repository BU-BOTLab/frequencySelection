function fit = fitMu_newOpts_display(YDATA, WT, freq, n, dist,name)
%%%byh Here the calibrated fdpm data is used along with the p1 model to
%%%recover optical properties at each of the fdpm wavelengths. Following
%%%that fit, the scattering is fit to a power law.


ndiodes=size(YDATA,2);
%fitoptions=optimset('display','off','TolFun',1e-16,'LargeScale','off');
fitoptions=optimset('display','iter','TolFun',1e-35,'LargeScale','off','TolX', 1e-14,'MaxIter', 10000,'maxfunevals',30000,'algorithm', 'levenberg-marquardt');
optimoptions(@lsqcurvefit,'OptimalityTolerance',1e-14,'StepTolerance', 1e-14);
%fitoptions=optimset('display','off');
%%%byh 5 sets of initial guesses to use for the fdpm fit
pinitial =	[.005 .8;.001 1.3;.01 1.0;.05 1.0;.005 0.6;.001 .1; 0.05 2];	%initial guess for non-linear least squares (2xnum_guesses)
%pinitial = [.001,.1; .001,3;.05,.5;.05,3;.025,1.55];
%    muaAll = [0.001,0.05];
%musAll = [.1,3];
%pinitial =	[.005 .8;.001 1.3;.01 1.0;.05 1.0;.02 0.6];	    %REI 3/29/2018
%pinitial =	[.01 1.0;.005 .8;.001 1.3];
%pinitial = [.01 1.0];
numStarts = size(pinitial,1);
PFIX = 0; %currently can't fix mua or mus
r2=0; %2 distance not coded
m = length(YDATA);
reimFlag = 1; %Fit real and imaginary parts (instead of amplitudes and phases)
%n=2;
flen = length(freq);
noWt = ones(size(WT(:,1)));
reff_option = 1;

%%%byh loop over each fdpm wavelength and do fdpm fit
for a = 1:ndiodes
    P1 = zeros(size(pinitial,1),2); RESID = zeros(m,size(pinitial,1)); CONVERGED = zeros(1,size(pinitial,1)); 
    CHI=zeros(size(pinitial,1),1); 
    %    OUTPUT=repmat(struct('iterations',0,'funcCount',0,'stepsize',0,'cgiterations','','firstorderopt',0,'algorithm',0,'message',''),1,size(pinitial,1));
    %%%byh loop over each initial guess
    for t = 1:size(pinitial,1)    
            %%%least squares fit to model using calibration.  P1 variable
            %%%is the set of recovered optical properties
        [P1(t,:), ss, RESID(:,t), CONVERGED(t)] = ...
            lsqcurvefit('p1seminf', pinitial(t,:), freq, YDATA(:,a), [], [], fitoptions, PFIX, ...
            n, dist, r2, WT(:,a), reimFlag, reff_option);
        CHI(t) =  ss;
        pause;
    end
    %     for cons scat
    
    %find lowest chi2 and retain the fit
    [fit.chi(a), best] = min(CHI);
    numex = 0;
    while 10*P1(best,1) > P1(best,2) && numex < numStarts
        numex = numex+1;
        CHI(best)= 1e30;
        [fit.chi(a),best] = min(CHI);
    end
 	disp(P1)
    %calculate good/bad fit index (AEC)
    fit.gbfi(a,1)=sum(abs(RESID(1:flen,best)))./(flen-1);
    fit.gbfi(a,2)=sum(abs(RESID(flen+1:2*flen,best)))./(flen-1);
    if numex < numStarts
        [fit.mua(a), fit.mus(a)] = ...
                deal(P1(best,1),P1(best,2));
    else
        fit.mua(a) = 0;
        fit.mus(a) = 0;
    end
    %Calculate the theoretical amplitudes and phases based on the fit
    theo = p1seminf([P1(best,1),P1(best,2)],freq,0,n,dist,0, noWt, reimFlag);
    fit.real(:,a) = theo(1:length(freq));
    fit.imag(:,a) = theo(length(freq)+1:end);
    fit.AC(:,a) = sqrt(fit.real(:,a).^2 + fit.imag(:,a).^2);
    fit.phase(:,a) = atan2(fit.imag(:,a),fit.real(:,a));
    fit.name = name;
end



