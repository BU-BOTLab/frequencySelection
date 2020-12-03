%%
addpath(genpath('utilities'))
addpath(genpath('dosiTrainingCode2019'))
%Load the data
dataDir = '../experimentalData/allDigitalSystem/allPhantoms30mm/200316';
%Load fitDat which is used for setting the "gold standard" OPs
load('../experimentalData/allDigitalSystem/allPhantomOPs.mat')

pName = cell(1,length(1:200:length(fitDat)));
muas = zeros(1,length(1:200:length(fitDat)));
muss = zeros(size(muas));
it = 1;
for i = 1:200:length(fitDat)
    thisSampIdx = i:i+199;
    pName{it} = fitDat(i).name;
    sampMua = zeros(1,200);
    sampMus = zeros(1,200);
    
    for j = 1:length(thisSampIdx)
        sampMua(j) = fitDat(thisSampIdx(j)).mua(5);
        sampMus(j) = fitDat(thisSampIdx(j)).mus(5);
    end
    muas(it) = mean(sampMua);
    muss(it) = mean(sampMus);
    it = it+1;
end
KMRef = 1 + muas./muss - sqrt((2*muas)./muss + muas.^2./muss.^2);
[sortedRef,idxRef] = sort(KMRef,'descend');
[sortedMua, idxMua] = sort(muas);
sortedMus = muss(idxMua);
%sortedMua = muas(idx);
%sortedName = pName(idx);

fitFreqs = 43:6:241;
ampVMua = zeros(length(sortedMua),length(fitFreqs));
phaseVMua = zeros(size(ampVMua));
ampmuVMua= zeros(size(ampVMua));
phasemuVMua=zeros(size(ampVMua));
rawDat = zeros(length(sortedMua),2);
for p = 1:length(sortedMua)
    thisSamp = strsplit(pName{p},'-');
    sampStr = [thisSamp{1},'*',thisSamp{3},'.asc'];
    [ampsd,phasesd,ampfit,phasefit,f] = getAmpPhaseSD(dataDir,sampStr,34,5,fitFreqs,0);
    [ampMu, phaseMu, ampFitMu, phaseFitMu, ~] = getAmpPhaseMu(dataDir,sampStr,34,5,fitFreqs,0);
    ampVMua(p,:)=ampfit(:,5);
    phaseVMua(p,:)=phasefit(:,5);
    ampmuVMua(p,:) = ampFitMu(:,5);
    phasemuVMua(p,:) = phaseFitMu(:,5);
    rawDat(p,:) = p1seminf_mba([muas(p),muss(p)],43,1.4,30,0,1);
end
%%
%It looks like there's a correlation between mua and amp/phase noise
%The trend is roughly linear and doesn't really depend on mus
%Let's try to fit a line to each frequency and compare slopes/intercepts
lamp = zeros(size(ampVMua,2),2);
lphase = zeros(size(lamp));
lamp0= zeros(size(ampVMua,2),1);
lphase0= zeros(size(lamp0));
% figure
rsq = zeros(1,size(ampVMua,2));
rsqPhase = zeros(size(rsq));
corrMuaAmp = zeros(size(rsq));
corrMusAmp = zeros(size(rsq));
corrMuaPhase = zeros(size(rsq));
corrMusPhase = zeros(size(rsq));
for i = 1:size(ampVMua,2)
   %DM = [muas',zeros(size(muas'))];
   la = muas'\ampVMua(:,i);
   lp = muas([1:6,8:17])'\phaseVMua([1:6,8:17],i);
   lamp0(i) = la(1);
   lphase0(i) = lp(1);
   lamp(i,:) = polyfit(muas',ampVMua(:,i),1);
   mdl = fitlm(muas',ampVMua(:,i));
   ca = corrcoef(muas,ampVMua(:,i));
   cs = corrcoef(muss,ampVMua(:,i));
   cap = corrcoef(muas,phaseVMua(:,i));
   csp = corrcoef(muss,phaseVMua(:,i));
   corrMuaAmp(i) = ca(1,2);
   corrMusAmp(i) = cs(1,2);
   corrMuaPhase(i) = cap(1,2);
   corrMusPhase(i) = csp(1,2);
   rsq(i)=mdl.Rsquared.Adjusted;
   mdl2 = fitlm(muas',phaseVMua(:,i),1);
   rsqPhase(i) = mdl2.Rsquared.Adjusted;
   lphase(i,:) = polyfit(muas([1:6,8:17])',phaseVMua([1:6,8:17],i),1);
   plot(muas,ampVMua(:,i),'o')
%    hold on
%    plot([0,.05],[0,.05]*la,'k--','linewidth',1.5)
%    pause()
   
end
Rsquared = mean(rsq)
RsquaredPhase = mean(rsqPhase)
corrMuaMean = mean(corrMuaAmp)
corrMusMean = mean(corrMusAmp)
corrMuaPhaseMean = mean(corrMuaPhase)
corrMusPhaseMean = mean(corrMusPhase)
mlamp = mean(lamp);
mphase = mean(lphase);

figure
plot(muas,phaseVMua,'o')
hold on
plot([.001,.05],[0.001,.05]*mphase(1)+mphase(2),'k--')

figure
plot(f,ampVMua)

l1 = polyfit(muas',ampVMua(:,12),1);
figure
subplot(121)
plot(muas,ampVMua(:,12),'o')
hold on
plot([0.001,0.035],l1(1)*[0.001,0.035]+l1(2),'k--')
title('Absorption vs. amplitude SD')
xlabel('\mu_a (1/mm)')
ylabel('\sigma_A')
subplot(122)
plot(muss,ampVMua(:,12),'o')
title('Scattering vs. ampltidue SD')
xlabel('\mu_s'' (1/mm)')
ylabel('\sigma_A')
print('../plots/amplitudeCorrelation.png','-dpng')

l2 = polyfit(muas([1:6,8:17])',phaseVMua([1:6,8:17],12),1);
figure
subplot(121)
plot(muas,phaseVMua(:,12),'o')
hold on
plot([0.001,0.035],l2(1)*[0.001,0.035]+l2(2),'k--')
title('Absorption vs. phase SD')
xlabel('\mu_a (1/mm)')
ylabel('\sigma_\phi')
subplot(122)
plot(muss,phaseVMua(:,12),'o')
title('Scattering vs. phase SD')
xlabel('\mu_s'' (1/mm)')
ylabel('\sigma_\phi')
print('../plots/phaseCorrelation.png','-dpng')

figure
plot(repmat(muas,1,length(f)),ampVMua(:),'o')

figure
plot(ampmuVMua(idxRef,1))
hold on
plot(ampmuVMua(idxMua,1))

figure
plot(sortedRef, ampmuVMua(:,1),'o')
figure
plot(sortedMua,ampmuVMua(:,1),'o')

figure
plot(ampmuVMua(:),ampVMua(:),'o')

figure
plot(phasemuVMua(:,1),phaseVMua(:,1),'o')

%%
%Let's fit the data
l1 = polyfit(log(ampmuVMua(:)), log(ampVMua(:)),1);
smX = linspace(min(log(ampmuVMua(:))), max(log(ampmuVMua(:))), 1000);
smY = l1(2) + l1(1)*smX;
figure
plot(log(ampmuVMua(:)),log(ampVMua(:)),'o')
hold on
plot(smX,smY,'k--')

figure
plot(f,ampVMua,'o')
