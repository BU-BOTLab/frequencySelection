%makePaperFigures.m
%
%Script for reproducing figures that appear in the multifrequency paper

%Add forward and inverse model to path
addpath(genpath('dosiTrainingCode2019'))
addpath(genpath('utilities'))

%Load any data
load('../generatedData/noNoise_LMA_optimized.mat')
noNoiseOpt = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/noNoise_LMA_default.mat')
noNoiseDef = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/noNoise_LUT.mat')
noNoiseLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

load('../generatedData/Digital_noise_850nm_equalData_optimized_v2.mat')
digitalEqualDatOpt = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/Digital_noise_850nm_unequalData_optimized_v2.mat')
digitalUnEqualDatOpt = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/Digital_noise_850nm_equalData_default_v2.mat')
digitalEqualDatDef = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/Digital_noise_850nm_unequalData_default_v2.mat')
digitalUnEqualDatDef = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/Digital_noise_850nm_equalData_LUT_v2.mat')
digitalEqualDatLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/Digital_noise_850nm_equalData_HDLUT_v2.mat')
digitalEqualDatHDLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/Digital_noise_850nm_unequalData_LUT_v2.mat')
digitalUnEqualDatLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

load('../generatedData/NA_noise_850nm_equalData_optimized.mat')
naEqualDatOpt = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_noise_850nm_unequalData_optimized.mat')
naUnEqualDatOpt = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_noise_850nm_equalData_default.mat')
naEqualDatDef = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_noise_850nm_unequalData_default.mat')
naUnEqualDatDef = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_noise_850nm_equalData_LUT.mat')
naEqualDatLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_noise_850nm_equalData_HDLUT.mat')
naEqualDatHDLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_noise_850nm_unequalData_LUT.mat')
naUnEqualDatLUT = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

%High density LUT 
load('../generatedData/noNoise_110MHz_HDLUT_2k.mat')
eucErrorNoNoise110 = sqrt(((trueMua-recMua)./trueMua).^2 + (((trueMus-recMus)./trueMus)).^2);

load('../generatedData/noNoise_500MHz_HDLUT_2k.mat')
eucErrorNoNoise500 = sqrt(((trueMua-recMua)./trueMua).^2 + (((trueMus-recMus)./trueMus)).^2);

load('../generatedData/OPSet10k.mat')
%%
%Figure 2
f=figure;
set(f,'Position',[660,300,900,840])
%sgtitle('No noise error single freq (tolFun = 1e-20)')

[ha,pa] = tight_subplot(2,2,[.15 .13],[.2 .16],[.1 .13]);

axes(ha(1))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),noNoiseOpt(:,2,3),4,log10(noNoiseOpt(:,2,3)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-15,0])
ylabel(c,'log_{10}(Error)')
title([{'Optimized parameters'},{'110 MHz, 30 mm'}])
set(gca,'fontsize',18)
view(0,90)

axes(ha(2))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),noNoiseDef(:,2,3),4,log10(noNoiseDef(:,2,3)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-15,0])
ylabel(c,'log_{10}(Error)')
title([{'Default parameters'},{'110 MHz, 30 mm'}])
set(gca,'fontsize',18)
view(0,90)

axes(ha(3))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),noNoiseOpt(:,6,3),4,log10(noNoiseOpt(:,6,3)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-15,0])
ylabel(c,'log_{10}(Error)')
title('50:500 MHz, 30 mm')
set(gca,'fontsize',18)
view(0,90)

axes(ha(4))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),noNoiseDef(:,6,3),4,log10(noNoiseDef(:,6,3)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-15,0])
ylabel(c,'log_{10}(Error)')
title('50:500 MHz, 30 mm')
set(gca,'fontsize',18)
view(0,90)
print('../plots/Fig2_defaultVsOptimized_v2.png','-dpng')

%%
%Figure 3
load('../generatedData/inversionTimingData.mat')
invPerSec = numTrials./[timingLUT; timingOpt; timingDefault];
invPerSecHD = numTrials./timingHDLUT;
figure
l1=semilogy(numFreqs, invPerSec(1,:)','-o','markerfacecolor',[77,175,77]/255,'color',[77,175,77]/255);
hold on
l4=semilogy(numFreqsHD,invPerSecHD,'-o','markerfacecolor',[216,164,155]/255,'color',[216,164,155]/255);
l2=semilogy(numFreqs, invPerSec(2,:)','-o','color',[55,126,184]/255,'markerfacecolor',[55,126,184]/255);
l3=semilogy(numFreqs, invPerSec(3,:)','-o','color',[228,26,28]/255,'markerfacecolor',[228,26,28]/255);
c1=get(l1,'color');
c2=get(l2,'color');
c3=get(l3,'color');
set(l1,'markerfacecolor',c1)
set(l2,'markerfacecolor',c2)
set(l3,'markerfacecolor',c3)
xlim([-10,500])
xticks([1,100:100:500])
xlabel('Number of frequencies')
ylabel('Inversions/second')
yticks([10,100,1000])
title('Inversion speed comparison')
legend('LUT','HD-LUT','Optimized','Default')
set(gca,'fontsize',18)
print('../plots/Fig3_inversionTimingPlot.png','-dpng')

figure
plot(numFreqs,invPerSec(3,:),'o')
lm = polyfit(numFreqs,invPerSec(3,:),1);
hold on
plot([1,1000], [1,1000].*lm(1) + lm(2))

figure
semilogy(numFreqs,invPerSec(3,:),'o')
hold on
plot(1:500, ([1:500].^lm(1))*exp(lm(2)))

%%
%Figure 4

%Data directories containing raw data for building noise model
dataDir10='../experimentalData/allDigitalSystem/Gen2_10mm/8192/200614';
dataDir20='../experimentalData/allDigitalSystem/Gen2_20mm/8192/200614';
dataDir30='../experimentalData/allDigitalSystem/Gen2_30mm/8192/200614';
dataDirM ='../experimentalData/networkAnalyzerSystem/SystemNoiseOP_Review/201027';
%Samples to use for noise models
sampNameLowMua = 'bpav4*.asc';
sampNameHighMua = 'b4h2*.asc';

%Acquisition parameters for digital system
numDiodes = 6;
numMeasFreqs = 70;
nind = 1.4;
allFa= 50:500;
mDosCol = [126,47,142]/255;
%Get noise model for digital system
[~, ~, ampfit10,phasefit10,~]=getAmpPhaseSD(dataDir10,sampNameLowMua,numMeasFreqs,numDiodes,allFa,5);
[~, ~, ampfit20,phasefit20,~]=getAmpPhaseSD(dataDir20,sampNameLowMua,numMeasFreqs,numDiodes,allFa,5);
[~, ~, ampfit30,phasefit30,~]=getAmpPhaseSD(dataDir30,sampNameLowMua,numMeasFreqs,numDiodes,allFa,5);

[~, ~, ampfit10high,phasefit10high,~]=getAmpPhaseSD(dataDir10,sampNameHighMua,numMeasFreqs,numDiodes,allFa,5);
[~, ~, ampfit20high,phasefit20high,~]=getAmpPhaseSD(dataDir20,sampNameHighMua,numMeasFreqs,numDiodes,allFa,5);
[~, ~, ampfit30high,phasefit30high,~]=getAmpPhaseSD(dataDir30,sampNameHighMua,numMeasFreqs,numDiodes,allFa,5);
%Get noise model for network analyzer system
[~, ~, ampfitM,phasefitM,~]=getAmpPhaseSD(dataDirM,sampNameLowMua,401,5,allFa,15);
[~, ~, ampfitMHigh,phasefitMHigh,~]=getAmpPhaseSD(dataDirM,sampNameHighMua,401,5,allFa,15);


f=figure;
set(f,'Position',[300,660,1740,640])
subplot(121)
semilogy(allFa, ampfit10(:,6)*100,'linewidth',2,'color',[16,78, 139]/255)
hold on
plot(allFa, ampfit10high(:,6)*100,'--','linewidth',2,'color',[16,78,139]/255)
plot(allFa, ampfit20(:,6)*100,'linewidth',2,'color', [178,34,34]/255)
plot(allFa, ampfit20high(:,6)*100,'--','linewidth',2,'color',[178,34,34]/255)
plot(allFa, ampfit30(:,6)*100,'linewidth',2,'color', [51,160,44]/255)
plot(allFa, ampfit30high(:,6)*100,'--','linewidth',2,'color',[51,160,44]/255)
plot(allFa, ampfitM(:,5)*100,'linewidth',2,'color',mDosCol)
plot(allFa, ampfitMHigh(:,5)*100,'--','linewidth',2,'color',mDosCol)
z(1)=plot(nan,nan,'linewidth',2,'color',[16,78, 139]/255);
z(2)=plot(nan,nan,'linewidth',2,'color',[178,34,34]/255);
z(3)=plot(nan,nan,'linewidth',2,'color',[51,160,44]/255);
z(4)=plot(nan,nan,'linewidth',2,'color',mDosCol);
z(5)=plot(nan,nan,'k','linewidth',2);
z(6)=plot(nan,nan,'k--','linewidth',2);
legend(z,'digital 10 mm', 'digital 20 mm', 'digital 30 mm', 'NA 28 mm','low \mu_a','high \mu_a','location','northwest')
xlabel('Frequency (MHz)')
ylabel('\sigma_{Amp} (%)')
title('Amplitude noise model')
set(gca,'fontsize',18)
subplot(122)
semilogy(allFa, phasefit10(:,6),'linewidth',2,'color',[16,78, 139]/255)
hold on
plot(allFa, phasefit10high(:,6),'--','linewidth',2,'color',[16,78,139]/255)
plot(allFa, phasefit20(:,6),'linewidth',2,'color', [178,34,34]/255)
plot(allFa, phasefit20high(:,6),'--','linewidth',2,'color',[178,34,34]/255)
plot(allFa, phasefit30(:,6),'linewidth',2,'color', [51,160,44]/255)
plot(allFa, phasefit30high(:,6),'--','linewidth',2,'color',[51,160,44]/255)
plot(allFa, phasefitM(:,5),'linewidth',2,'color',mDosCol)
plot(allFa, phasefitMHigh(:,5),'--','linewidth',2,'color',mDosCol)
%legend('digital 10 mm', 'digital 20 mm', 'digital 30 mm', 'NA 28 mm','location','northwest')
xlabel('Frequency (MHz)')
ylabel('\sigma_{Phase} (rad)')
title('Phase noise model')
set(gca,'fontsize',18)
print('../plots/Fig4_errorModelSDSep_v4.png','-dpng')

%%
%Figure 5
f=figure;
set(f,'Position',[762,805,800,533])
[ha,pa] = tight_subplot(1,2,[.07 .07],[.25 .1],[.07 .03]);
axes(ha(1));
[h,L,MX,MED,bw]=violin(log10(digitalEqualDatOpt(:,:,3)),'x',[1,4,7,10,13,16],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10(digitalUnEqualDatOpt(:,:,3)),'x',[2,5,8,11,14,17],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
xticks([1.5:3:16.5])
delete(L)
xticklabels({'70','110','500','50,500','50:7:253','50:499'})
xtickangle(-45)
ax=gca;
set(ax,'TickLength',[.01,.01])
ax.XAxis
xlim([0,18])
ylim([-6,2])
xlabel('Frequency Set (MHz)')
ylabel('log_{10}(Error)')
h = zeros(2, 1);
h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
legend(h,'Equal Data','Unequal Data','fontsize',16);
set(gca,'fontsize',18)
title('All Digital FD-DOS')

axes(ha(2));
set(f,'Position',[762,805,1600,533])
[h,L,MX,MED,bw]=violin(log10(naEqualDatOpt),'x',[1,4,7,10,13,16],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10(naUnEqualDatOpt),'x',[2,5,8,11,14,17],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
xticks([1.5:3:16.5])
delete(L)
xticklabels({'70','110','500','50,500','50:7:253','50:500'})
xtickangle(-45)
ax=gca;
set(ax,'TickLength',[.01,.01])
ax.XAxis
xlim([0,18])
ylim([-6,2])
xlabel('Frequency Set (MHz)')
ylabel('log_{10}(Error)')
h = zeros(2, 1);
h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
legend(h,'Equal Data','Unequal Data','fontsize',16);
set(gca,'fontsize',18)
title('Network Analyzer FD-DOS')
print('../plots/Fig5_violinPlots_systemComp_v3.png','-dpng')

%%
%Figure 6

f=figure;
noNoiseOpt(noNoiseOpt == 0) = eps;
noNoiseDef(noNoiseDef == 0) = eps;
noNoiseLUT(noNoiseLUT == 0) = eps;
set(f,'Position',[762,100,1200,1066])
[ha,pa] = tight_subplot(2,2,[.1 .1],[.15 .1],[.1 .03]);
%110 MHz
axes(ha(1));
[h,L,MX,MED,bw]=violin(log10([noNoiseOpt(:,2,3),digitalEqualDatOpt(:,2,3),naEqualDatOpt(:,2,1)]),'x',[1,6,11],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseDef(:,2,3),digitalEqualDatDef(:,2,3),naEqualDatDef(:,2,1)]),'x',[2,7,12],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseLUT(:,2,3),digitalEqualDatLUT(:,2,1),naEqualDatLUT(:,2,1)]),'x',[3,8,13],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([eucErrorNoNoise110',digitalEqualDatHDLUT(:,2,1),naEqualDatHDLUT(:,2,1)]),'x',[4,9,14],'facecolor',[216,164,155]/255,'facealpha',1,'medc','c','mc','k');
xticks([2.5,7.5,12.5])
delete(L)
xticklabels({'No Noise','Digital FD-DOS','NA FD-DOS'})
%xtickangle(-45)
ax=gca;
set(ax,'TickLength',[.01,.01])

xlim([0,15])
ylim([-18,4])
xlabel('Noise Model')
ylabel('log_{10}(Error)')

h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
h(4) = plot(nan,nan,'s','markerfacecolor',[216,164,155]/255,'markersize',15);
legend(h,'Optimized','Default','LUT','HD-LUT','fontsize',16,'Location','southeast');
set(gca,'fontsize',18)
title('110 MHz')
%50:500
axes(ha(2));
[h,L,MX,MED,bw]=violin(log10([noNoiseOpt(:,6,2),digitalEqualDatOpt(:,6,2),naEqualDatOpt(:,6,1)]),'x',[1,5,9],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseDef(:,6,2),digitalEqualDatDef(:,6,2),naEqualDatDef(:,6,1)]),'x',[2,6,10],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseLUT(:,6,2),digitalEqualDatLUT(:,6,1),naEqualDatLUT(:,6,1)]),'x',[3,7,11],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
xticks([2,6,10])
delete(L)
xticklabels({'No Noise','Digital FD-DOS','NA FD-DOS'})
%xtickangle(-45)
ax=gca;
set(ax,'TickLength',[.01,.01])
xlim([0,12])
ylim([-18,4])
xlabel('Noise Model')
ylabel('log_{10}(Error)')
h = zeros(2, 1);
h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
legend(h,'Optimized','Default','LUT','fontsize',16,'Location','southeast');
set(gca,'fontsize',18)
title('50 to 500 MHz by 1 MHz')
%500 MHz
axes(ha(3));
[h,L,MX,MED,bw]=violin(log10([noNoiseOpt(:,3,3),digitalEqualDatOpt(:,3,3),naEqualDatOpt(:,3,1)]),'x',[1,6,11],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseDef(:,3,3),digitalEqualDatDef(:,3,3),naEqualDatDef(:,3,1)]),'x',[2,7,12],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseLUT(:,3,3),digitalEqualDatLUT(:,3,1),naEqualDatLUT(:,3,1)]),'x',[3,8,13],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([eucErrorNoNoise500',digitalEqualDatHDLUT(:,3,1),naEqualDatHDLUT(:,3,1)]),'x',[4,9,14],'facecolor',[216,164,155]/255,'facealpha',1,'medc','c','mc','k');
xticks([2.5,7.5,12.5])
delete(L)
xticklabels({'No Noise','Digital FD-DOS','NA FD-DOS'})
%xtickangle(-45)
ax=gca;
set(ax,'TickLength',[.01,.01])

xlim([0,15])
ylim([-18,4])
xlabel('Noise Model')
ylabel('log_{10}(Error)')

h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
h(4) = plot(nan,nan,'s','markerfacecolor',[216,164,155]/255,'markersize',15);
legend(h,'Optimized','Default','LUT','HD-LUT','fontsize',16,'Location','southeast');
set(gca,'fontsize',18)
title('500 MHz')

%50:500
axes(ha(4));
[h,L,MX,MED,bw]=violin(log10([noNoiseOpt(:,5,2),digitalEqualDatOpt(:,5,2),naEqualDatOpt(:,5,1)]),'x',[1,5,9],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseDef(:,5,2),digitalEqualDatDef(:,5,2),naEqualDatDef(:,5,1)]),'x',[2,6,10],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseLUT(:,5,2),digitalEqualDatLUT(:,5,1),naEqualDatLUT(:,5,1)]),'x',[3,7,11],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
xticks([2,6,10])
delete(L)
xticklabels({'No Noise','Digital FD-DOS','NA FD-DOS'})
%xtickangle(-45)
ax=gca;
set(ax,'TickLength',[.01,.01])
xlim([0,12])
ylim([-18,4])
xlabel('Noise Model')
ylabel('log_{10}(Error)')
h = zeros(2, 1);
h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
legend(h,'Optimized','Default','LUT','fontsize',16,'Location','southeast');
set(gca,'fontsize',18)
title('50 to 253 MHz by 7 MHz')

print('../plots/Fig6_inversionMethodAccuracy_v4.png','-dpng')
%%
%Updated violin plot figure to more easily show the differences between
%methods as requested by the reviewer

f=figure;
noNoiseOpt(noNoiseOpt == 0) = eps;
noNoiseDef(noNoiseDef == 0) = eps;
noNoiseLUT(noNoiseLUT == 0) = eps;
set(f,'Position',[762,100,1500,533])
[ha,pa] = tight_subplot(1,3,[.0 .05],[.18 .1],[.06 .03]);
%110 MHz
axes(ha(1));
[h,L,MX,MED,bw]=violin(log10([noNoiseOpt(:,2,3),noNoiseOpt(:,3,3),noNoiseOpt(:,5,3),noNoiseOpt(:,6,3)]),'x',[1,7,13,17],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseDef(:,2,3),noNoiseDef(:,3,3),noNoiseDef(:,5,3),noNoiseDef(:,6,3)]),'x',[2,8,14,18],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([noNoiseLUT(:,2,3),noNoiseLUT(:,3,3),noNoiseLUT(:,5,3),noNoiseLUT(:,6,3)]),'x',[3,9,15,19],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([eucErrorNoNoise110',eucErrorNoNoise500']),'x',[4,10],'facecolor',[216,164,155]/255,'facealpha',1,'medc','c','mc','k');
xticks([2.5,8.5,14,18])
delete(L)
xticklabels({'110','500','50:7:253','50:1:500'})
xtickangle(-15)
ax=gca;
set(ax,'TickLength',[.01,.01])
xlim([0,20])
ylim([-18,4])
xlabel('Frequencies (MHz)')
ylabel('log_{10}(Error)')
% 
% h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
% h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
% h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
% h(4) = plot(nan,nan,'s','markerfacecolor',[216,164,155]/255,'markersize',15);
% legend(h,'Optimized','Default','LUT','HD-LUT','fontsize',16,'Location','southeast');
set(gca,'fontsize',18)
title('No added noise')
%50:500
axes(ha(2));
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatOpt(:,2,3),digitalEqualDatOpt(:,3,3),digitalEqualDatOpt(:,5,3),digitalEqualDatOpt(:,6,3)]),'x',[1,7,13,17],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatDef(:,2,3),digitalEqualDatDef(:,3,3),digitalEqualDatDef(:,5,3),digitalEqualDatDef(:,6,3)]),'x',[2,8,14,18],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatLUT(:,2),digitalEqualDatLUT(:,3),digitalEqualDatLUT(:,5),digitalEqualDatLUT(:,6)]),'x',[3,9,15,19],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatHDLUT(:,2),digitalEqualDatHDLUT(:,3)]),'x',[4,10],'facecolor',[216,164,155]/255,'facealpha',1,'medc','c','mc','k');
xticks([2.5,8.5,14,18])
delete(L)
xticklabels({'110','500','50:7:253','50:1:500'})
xtickangle(-15)
ax=gca;
set(ax,'TickLength',[.01,.01])
xlim([0,20])
ylim([-7,4])
xlabel('Frequencies (MHz)')
ylabel('log_{10}(Error)')
h = zeros(2, 1);
h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
h(4)=plot(nan,nan,'s','markerfacecolor',[216,164,155]/255,'markersize',15);
legend(h,'Optimized','Default','LUT','HD-LUT','fontsize',16,'Location','northeast');
set(gca,'fontsize',18)
title('All-digital system noise')
%500 MHz
axes(ha(3));

[h,L,MX,MED,bw]=violin(log10([naEqualDatOpt(:,2),naEqualDatOpt(:,6),naEqualDatOpt(:,3),naEqualDatOpt(:,5)]),'x',[1,7,13,17],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([naEqualDatDef(:,2),naEqualDatDef(:,6),naEqualDatDef(:,3),naEqualDatDef(:,5)]),'x',[2,8,14,18],'facecolor',[55,126,184]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([naEqualDatLUT(:,2),naEqualDatLUT(:,6),naEqualDatLUT(:,3),naEqualDatLUT(:,5)]),'x',[3,9,15,19],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([naEqualDatHDLUT(:,1),naEqualDatHDLUT(:,3)]),'x',[4,10],'facecolor',[216,164,155]/255,'facealpha',1,'medc','c','mc','k');
xticks([2.5,8.5,14,18])
delete(L)
xticklabels({'110','500','50:7:253','50:1:500'})
xtickangle(-15)
ax=gca;
set(ax,'TickLength',[.01,.01])
xlim([0,20])
ylim([-7,4])
xlabel('Frequencies (MHz)')
ylabel('log_{10}(Error)')
% h = zeros(2, 1);
% h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
% h(2)=plot(nan,nan,'s','markerfacecolor',[55,126,184]/255,'markersize',15);
% h(3)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
% legend(h,'Optimized','Default','LUT','fontsize',16,'Location','best');
set(gca,'fontsize',18)
title('Network analyzer system noise')

print('../plots/Fig6_inversionMethodAccuracy_bySystem_v2.png','-dpng')

%%
%Fig XX, plot of "optimal frequencies"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%850 high mua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../generatedData/digital_850nm_optimalFreq_30mm.mat')
wavelength = 850;
eucErrordDOS = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);
load('../generatedData/NA_850nm_optimalFreq.mat')
eucErrormDOS = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2 + (((saveMusTrue-saveMusRec)./saveMusTrue)).^2);

eucErrordDOS(1059,19) = median(eucErrordDOS(:,19)); %This point is a ridiculously extreme outlier, not sure why
eucErrordDOS(1795,19) = median(eucErrordDOS(:,19)); %This point is a ridiculously extreme outlier, not sure why

distMua = sqrt(((saveMuaTrue-saveMuaRec)./saveMuaTrue).^2);
distMus = sqrt(((saveMusTrue-saveMusRec)./saveMusTrue).^2);
avgErrDOS = mean(eucErrordDOS,1);
avgErrMDOS = mean(eucErrormDOS,1);
medErrDOS = median(eucErrordDOS,1);
medErrMDOS = median(eucErrormDOS,1);
stdErr = std(log10(eucErrordDOS),[],1);
stdErrM = std(log10(eucErrormDOS),[],1);
freqs = 50:20:490;

figure
e=plot(freqs,log10(squeeze(medErrDOS(1,1:length(freqs),1))),'o');
c=get(e,'color');
set(e,'markerfacecolor',c)
hold on
f=plot(freqs,log10(squeeze(medErrMDOS(1,1:length(freqs),1))),'o');
c2 = get(f,'color');
set(f,'markerfacecolor',c2);
xlabel('Modulation Frequency (MHz)')
ylabel('log_{10}(Median Error)')
title('Comparison of single frequencies')
xlim([40,425])
legend('all-digital','network analyzer')
set(gca,'fontsize',18)
print('../plots/optimalFrequencies.png','-dpng')


f=figure;
e=errorbar(freqs,log10(medErrDOS),stdErr,'o-','linewidth',1.5);
c=get(e,'color');
set(e,'markerfacecolor',c)
hold on
f=errorbar(freqs+6,log10(medErrMDOS),stdErrM,'o-','linewidth',1.5);
c2 = get(f,'color');
set(f,'markerfacecolor',c2);
xlabel('Modulation Frequency (MHz)')
ylabel('log_{10}(Median Error)')
title('Comparison of single frequencies')
xlim([40,510])
legend('all-digital','network analyzer','location','northwest')
set(gca,'fontsize',18)
print('../plots/optimalFrequencies_errBar.png','-dpng')

f=figure;
e=boxplot(log10(eucErrordDOS),'symbol','');
hold on
f=boxplot(log10(eucErrormDOS),'symbol','');
xlabel('Modulation Frequency (MHz)')
ylabel('log_{10}(Error)')
title('Comparison of single frequencies')
set(gca,'fontsize',18)
print('../plots/optimalFrequencies_boxPlot.png','-dpng')

%%
figure
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatOpt(:,3,3),naEqualDatOpt(:,3,1)]),'x',[1,5],'facecolor',[228,26,28]/255,'facealpha',1,'medc','c','mc','k');
hold on
delete(L)
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatLUT(:,3,1),naEqualDatLUT(:,3,1)]),'x',[2,6],'facecolor',[77,175,77]/255,'facealpha',1,'medc','c','mc','k');
delete(L)
[h,L,MX,MED,bw]=violin(log10([digitalEqualDatHDLUT(:,3,1),naEqualDatHDLUT(:,3,1)]),'x',[3,7],'facecolor',[216,164,155]/255,'facealpha',1,'medc','c','mc','k');
xticks([2,6])
delete(L)
xlim([0,8])
ylim([-6,1])
xticklabels({'High noise','Low noise'})
ylabel('log_{10} (OP Error)')
title('Look-up table accuracy')
set(gca,'fontsize',18)
h(1)=plot(nan,nan,'s','markerfacecolor',[228,26,28]/255,'markersize',15);
h(2)=plot(nan,nan,'s','markerfacecolor',[77,175,77]/255,'markersize',15);
h(3) = plot(nan,nan,'s','markerfacecolor',[216,164,155]/255,'markersize',15);
legend(h,'Iterative','LUT','HD-LUT','fontsize',16,'Location','southwest');


%%
%
load('../generatedData/Digital_noise_850nm_unEqualData_optimized_v2.mat')
errorDigOpt = sqrt((saveMuaTrue-saveMuaRec).^2 + (saveMusTrue-saveMusRec).^2);
recMuaDigOpt = saveMuaRec;
recMusDigOpt = saveMusRec;
%load('../generatedData/NA_noise_850nm_equalData_optimized.mat')

h=figure;
set(h,'Position',[415,840,1156,512])
subplot(121)
plot(saveMuaTrue(:,2),log10(errorDigOpt(:,2,3)),'.')
xlabel('True \mu_a (mm^{-1})')
ylabel('log_{10} (OP Error)')
title('Error vs. Absorption')
set(gca,'fontsize',18)
subplot(122)
plot(saveMusTrue(:,2),log10(errorDigOpt(:,2,3)),'.')
xlabel('True \mu_s'' (mm^{-1})')
ylabel('log_{10} (OP Error)')
title('Error vs. Scattering')
set(gca,'fontsize',18)
print('../plots/OPErrorAnalysis.png','-dpng')

figure
subplot(121)
plot(saveMuaTrue(:,1),saveMuaRec(:,1,3),'.')
hold on
plot([0,0.05],[0,0.05],'k--','linewidth',1.5)
ylim([0,.05])
subplot(122)
plot(saveMusTrue(:,1),saveMusRec(:,1,3),'.')

%%
f=figure;
set(f,'Position',[660,300,1550,450])
%sgtitle('No noise error single freq (tolFun = 1e-20)')

[ha,pa] = tight_subplot(1,3,[.1 .1],[.18 .16],[.05 .05]);

axes(ha(1))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),digitalEqualDatOpt(:,2,3),4,log10(digitalEqualDatOpt(:,2,3)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-4,-1])
ylabel(c,'log_{10}(Error)')
title('Optimized Iterative')
set(gca,'fontsize',18)
xlim([0,0.05])
view(0,90)

axes(ha(2))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),digitalEqualDatLUT(:,2),4,log10(digitalEqualDatLUT(:,2)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-4,-1])
xlim([0,0.05])
ylabel(c,'log_{10}(Error)')
title('LUT')
set(gca,'fontsize',18)
view(0,90)

axes(ha(3))
scatter3(saveMuaTrue(:,1,1),saveMusTrue(:,1,1),digitalEqualDatHDLUT(:,2),4,log10(digitalEqualDatHDLUT(:,2)),'filled')
xlabel('\mu_a (1/mm)')
ylabel('\mu_s'' (1/mm)')
c=colorbar;
caxis([-4,-1])
xlim([0,0.05])
ylabel(c,'log_{10}(Error)')
title('HD-LUT')
set(gca,'fontsize',18)
view(0,90)
print('../plots/LUT_landscape.png','-dpng')

%%
%Easier speed plot
%Figure 3
load('../generatedData/inversionTimingData.mat')
invPerSec = numTrials./[timingLUT; timingOpt; timingDefault];
invPerSecHD = numTrials./timingHDLUT;
figure
l1=semilogy(numFreqs, invPerSec(1,:)','-o','markerfacecolor',[77,175,77]/255,'color',[77,175,77]/255);
hold on
l4=semilogy(numFreqsHD,invPerSecHD,'-o','markerfacecolor',[216,164,155]/255,'color',[216,164,155]/255);
%l2=semilogy(numFreqs, invPerSec(2,:)','-o','color',[55,126,184]/255,'markerfacecolor',[55,126,184]/255);
l3=semilogy(numFreqs, invPerSec(3,:)','-o','color',[228,26,28]/255,'markerfacecolor',[228,26,28]/255);
c1=get(l1,'color');
%c2=get(l2,'color');
c3=get(l3,'color');
set(l1,'markerfacecolor',c1)
%set(l2,'markerfacecolor',c2)
set(l3,'markerfacecolor',c3)
xlim([-10,500])
xticks([1,100:100:500])
xlabel('Number of frequencies')
ylabel('Inversions/second')
yticks([10,100,1000])
title('Inversion speed comparison')
legend('LUT','HD-LUT','Iterative')
set(gca,'fontsize',18)
print('../plots/Fig3_inversionTimingPlot_simple.png','-dpng')


