function fit = mufit(diodes, model_to_fit, YDATA, WT, freq, nofr, r1, options, verbose, reff_option)
%%%byh Here the calibrated fdpm data is used along with the p1 model to
%%%recover optical properties at each of the fdpm wavelengths. Following
%%%that fit, the scattering is fit to a power law.


ndiodes=length(diodes);
%fitoptions=optimset('display','off','TolFun',1e-10,'LargeScale','off','LevenbergMarquardt','on');
fitoptions=optimset('display','off');
%%%byh 5 sets of initial guesses to use for the fdpm fit
pinitial =	[.005 .8;.001 1.3;.01 1.0;.05 1.0;.005 0.6];	%initial guess for non-linear least squares (2xnum_guesses)
%pinitial =	[.005 .8;.001 1.3;.01 1.0;.05 1.0;.02 0.6];	    %REI 3/29/2018
%pinitial =	[.01 1.0;.005 .8;.001 1.3];
%pinitial = [.01 1.0];
PFIX = 0; %currently can't fix mua or mus
r2=0; %2 distance not coded
m = length(YDATA);
n=2;
flen = length(freq);
noWt = ones(size(WT(:,1)));
holdmusflag=0;
musavenum=0;

if(isfield(options,'holdmus'))
    if(options.holdmus>0)
        persistent musholdval;
        holdmusflag=1;
        musavenum=options.holdmus;
        if options.itr==(musavenum+1)
            musholdval=musholdval/musavenum;
        end
    end
end

%%%byh loop over each fdpm wavelength and do fdpm fit
for a = 1:ndiodes
    P1 = zeros(size(pinitial,1),2); RESID = zeros(m,size(pinitial,1)); CONVERGED = zeros(1,size(pinitial,1)); 
    JACOBIAN = zeros(m,2,size(pinitial,1)); CHI=zeros(size(pinitial,1),1); 
    %    OUTPUT=repmat(struct('iterations',0,'funcCount',0,'stepsize',0,'cgiterations','','firstorderopt',0,'algorithm',0,'message',''),1,size(pinitial,1));
    %%%byh loop over each initial guess
    for t = 1:size(pinitial,1)
%         [P1(t,:), ss, RESID(:,t), CONVERGED(t), OUTPUT(t), lambda, JACOBIAN(:, :, t)] = ...
%             lsqcurvefit(model_to_fit, pinitial(t,:), freq, YDATA(:,a), [], [], fitoptions, PFIX, ...
%             nofr, r1, r2, WT(:,a), options.imagfit, reff_option);
        if holdmusflag==0 || options.itr<=musavenum
            %%%least squares fit to model using calibration.  P1 variable
            %%%is the set of recovered optical properties
            [P1(t,:), ss, RESID(:,t), CONVERGED(t), OUTPUT(t), lambda, JACOBIAN(:, :, t)] = ...
                lsqcurvefit(model_to_fit, pinitial(t,:), freq, YDATA(:,a), [0 0.2], [1 10], fitoptions, PFIX, ...
                nofr, r1, r2, WT(:,a), options.imagfit, reff_option); %********REI 4/12/2018 constrain the fit in mua and musp***********
%                  [P1(t,:), ss, RESID(:,t), CONVERGED(t), OUTPUT(t), lambda, JACOBIAN(:, :, t)] = ...
%                 lsqcurvefit(model_to_fit, pinitial(t,:), freq, YDATA(:,a), [], [], fitoptions, PFIX, ...
%                 nofr, r1, r2, WT(:,a), options.imagfit, reff_option);
        else
%             [P1(t,:), ss, RESID(:,t), CONVERGED(t), OUTPUT(t), lambda, JACOBIAN(:, :, t)] = ... %     for cons scat
%                 lsqcurvefit('p1seminfsetmus', pinitial(t,:), freq, YDATA(:,a), [0 0.2], [1 10], fitoptions, PFIX, ...
%                 nofr, r1, r2, WT(:,a), options.imagfit, reff_option, musholdval(a)); %********REI 4/12/2018 constrain the fit in mua and musp***********
                [P1(t,:), ss, RESID(:,t), CONVERGED(t), OUTPUT(t), lambda, JACOBIAN(:, :, t)] = ... %     for cons scat
                lsqcurvefit('p1seminfsetmus', pinitial(t,:), freq, YDATA(:,a), [], [], fitoptions, PFIX, ...
                nofr, r1, r2, WT(:,a), options.imagfit, reff_option, musholdval(a));
        end
        CHI(t) =  ss./(m-2);
%         if P1(t,1)>0
%             break
%         end
    end
    %     for cons scat
    if(musavenum>0)
        if options.itr==1
            musholdval(a)=P1(t,2);
        elseif options.itr<=musavenum
            musholdval(a)=P1(t,2)+musholdval(a);
        else
            P1(t,2)=musholdval(a);
        end
    end
    
    %find lowest chi2 and retain the fit
 	[fit.chi(a), best] = min(CHI);	
 	[fit.alg{a}, fit.iter(a), fit.guessnum(a), fit.guesses(a), fit.converged(a)] = ...
 		deal (OUTPUT(best).algorithm, OUTPUT(best).iterations, best, t, CONVERGED(best));
	%fit.chi(a) = CHI(t);	
	%[fit.alg{a}, fit.iter(a), fit.guessnum(a), fit.guesses(a), fit.converged(a)] = ...
	%	deal (OUTPUT(t).algorithm, OUTPUT(t).iterations, t, t, CONVERGED(t));
    
    
    %calculate good/bad fit index (AEC)
    fit.gbfi(a,1)=sum(abs(RESID(1:flen,best)))./(flen-1);
    fit.gbfi(a,2)=sum(abs(RESID(flen+1:2*flen,best)))./(flen-1);
%   fit.gbfi(a,1)=sum(abs(RESID(1:flen,t)))./(flen-1);
%   fit.gbfi(a,2)=sum(abs(RESID(flen+1:2*flen,t)))./(flen-1);
    
	%Assign to correct values
	%which parameters to fit?fd
 	cov = estimated_cov(n, m, RESID(:,best), JACOBIAN(:,:,best));
%    cov = estimated_cov(n, m, RESID(:,t), JACOBIAN(:,:,t));
	err = sqrt(diag(cov));

 	[fit.mua(a), fit.mus(a), fit.dmua(a), fit.dmus(a), fit.scov] = ...
 			deal(P1(best,1),P1(best,2),err(1),err(2), sqrt(cov(1,2)));
%	[fit.mua(a), fit.mus(a), fit.dmua(a), fit.dmus(a), fit.scov] = ...
%			deal(P1(t,1),P1(t,2),err(1),err(2), sqrt(cov(1,2)));


	%confidence interval: 
	if cov(2,2) ~= 0
		DP1=[fit.dmua(a); fit.dmus(a)] ;CPINV=inv(cov');dchi=DP1'*CPINV*DP1;
		fit.conf(a)=1-gammainc(2/2,dchi/2);		%I think (albert)
	else
		fit.conf(a) = 0;  %temp fix .. if one variable is fixed, get d2 =0, bad matrix inversion
	end;
	
	%Calculate fitted amp and phi
 	OUTP = feval(model_to_fit, P1(best,:), freq, 0, nofr,r1,r2, noWt, options.imagfit, reff_option);
%    OUTP = feval(model_to_fit, P1(t,:), freq, 0, nofr,r1,r2, noWt, options.imagfit, reff_option);
% OUTP = feval(model_to_fit, P1(best,:), freq, 0, n,r1, noWt, options.imagfit);
	%split up OUTPut from function into correct components
	if(options.imagfit == 0)		%amp/phase fit
		fit.amp(:,a)=OUTP(1:flen);										
		fit.phi(:,a)=OUTP(1+flen:2*flen);					
	else     %function gave Re and Im ...convert to phase and amplitude
		FIT_R=OUTP(1:flen);
		FIT_I=OUTP(1+flen:2*flen);
		fit.amp(:,a)=sqrt(FIT_R.^2+FIT_I.^2);
		fit.phi(:,a)=unwrap(atan2(FIT_I,FIT_R));
	end;
	
	%print the information of the fit
    if verbose
        form = '%d\t%7.5f \t%7.6f \t%7.3f \t%7.4f \t%5.2g\t%d in %d itt\t %d';
        disp(sprintf(form, diodes(a),fit.mua(a),fit.dmua(a), ...
            fit.mus(a),fit.dmus(a),fit.chi(a),fit.guessnum(a), fit.iter(a), fit.guesses(a)));
    end
end

if ndiodes>=2
    %%%byh fdpmFitScat does power law fit to recovered scattering.  The
    %%%resulting power law coefficients are then passed on
    [fit.preft, fit.slope, fit.dpreft, fit.dslope] =  ...
        fdpmFitScat(diodes', fit.mus, fit.chi, options, verbose);
end

fit.phy=[];
fit.fitmua=[];
fit.ss=[];


