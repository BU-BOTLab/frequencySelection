%DFDP  Numerical derivative calculator (Matlab site function)
% prt= Jacobian Matrix prt(i,j)=df(i)/dp(j)
% x=vec or matrix of indep var(used as arg to func) x=[x0 x1 ....]
% f=func(x,p) vector initialsed by user before each call to dfdp
% p= vec of current parameter values
% dp= fractional increment of p for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
% vargin:  variables sent into function
%
%NOTES: Version 1.5
% modified by AEC
%===================================================================================================
function prt=dfdp_mba(func, x,f,p,dp,varargin)


m=2*length(x);n=length(p);      					%dimensions
ps=p; prt=zeros(m,n);del=zeros(n,1);         % initialize Jacobian to Zero
for j=1:n
	del(j)=dp(j) .*p(j);		%cal delx=fract(dp)*param value(p)
	if p(j)==0
		del(j)=dp(j);		%if param=0 delx=fraction
	end
	p(j)=ps(j) + del(j);
	if del(j)~=0, f1=reshape(feval(func,p,x,varargin{:}),[],1);
		if dp(j) < 0, prt(:,j)=(f1-f)./del(j);
		else
			p(j)=ps(j)- del(j);
			prt(:,j)=(f1-reshape(feval(func,p,x,varargin{:}),[],1))./(2 .*del(j));
		end
	end
	p(j)=ps(j);      %restore p(j)
end
return
