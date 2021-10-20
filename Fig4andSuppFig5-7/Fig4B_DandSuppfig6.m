%% load data
tumormarker_data=readtable('tumormarkers_paper.csv');

%% calculate time to resistance (TTR)
% for TTR, estimate time to return to starting size for purposes of model

%patient 468 
TTR_Let468=(tumormarker_data.relativedays2(212)+tumormarker_data.relativedays2(213))/2-tumormarker_data.relativedays2(205);
TTR_NabPac468=tumormarker_data.relativedays2(222)-tumormarker_data.relativedays2(216);

%patient 288
TTR_Let288=tumormarker_data.relativedays2(103)-tumormarker_data.relativedays2(51);
TTR_Everol288=(tumormarker_data.relativedays2(111)+tumormarker_data.relativedays2(112))/2-tumormarker_data.relativedays2(104);

%patient 15
TTR_Gem15=tumormarker_data.relativedays2(28)-tumormarker_data.relativedays2(11);

%patient 422
TTR_Cap422=tumormarker_data.relativedays2(167)-tumormarker_data.relativedays2(152);
TTR_Erib422=tumormarker_data.relativedays2(190)-tumormarker_data.relativedays2(180);
TTR_Pac422=tumormarker_data.relativedays2(199)-tumormarker_data.relativedays2(190);



% plot patient data overlay with predictions
%% 468 letrozole
M=[-0.277601468039394,0;0.000100000000000000,0.0810090694623333];
Moff=[0.426907753811401,0;0.000100000000000000,0.0810090694623333];
n0=[55.8000000000000;6.20000000000000];
x=1.6731;
timeon=3.9669;
dayplot=(tumormarker_data.relativedays2(205:213)-19939);
%[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
[breakadvantage,time,dataplot,timeplot,trelapse,singledosen] = plotOptimalBreak_delayedstart(M,Moff,n0,x,timeon,dayplot(2));
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'k');
hold on
%timescale2=linspace(0,time,length(dataplot))*10;
%plot(timescale2,dataplot(:,3),'m')
plot(timeplot*10,dataplot(:,3),'m')
scatter(dayplot,tumormarker_data.CA_15_3num(205:213),'r','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') %284
plot([TTR_Let468,TTR_Let468],[0,maxy(2)],'k--') %277
plot([timeplot(end)*10,timeplot(end)*10],[0,maxy(2)],'k--') %469 454
hold off
xlabel('days')
ylabel('tumor burden')
maxad_486Let=454/284;
%% 468 NabPac
M=[-21.0155074411907,0;0.00100000000000000,0.0881259919218638];
Moff=[0.482768329871426,0;0.00100000000000000,0.0881259919218638];
n0=[388.842038724327;60.1579612756727];
x_in=round(21/2);
x=1.01;
timeon=1.04;
%[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,trelapse] = plotOptimalBreak_pulsevspulse(M,Moff,n0,x_in,x,1,20,timeon,1);
[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,timeplot,trelapse] = plotOptimalBreak_pulsevspulse_delayedstart(M,Moff,n0,x_in,x,1,20,timeon,1,4);

figure
ind=find(dataplot_in(:,3)<252);
timescale=linspace(0,trelapse,ind(end));
plot(timescale,dataplot_in(1:ind(end),3),'m');
hold on
%ind2=find(dataplot(:,3)<252);
% timescale2=linspace(0,time,ind2(end));
% plot(timescale2,dataplot(1:ind2(end),3),'m--')
plot(timeplot(1:length(dataplot))*10,dataplot(:,3),'m--')
ind3=find(singledosen<252);
timescale3=linspace(0,singledose_trelapse,ind3(end));
plot(timescale3,singledosen(1:ind3(end)),'k')
dayplot=(tumormarker_data.relativedays2(216:222)-20337);
scatter(dayplot,tumormarker_data.CA_15_3num(216:222),'r','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') %171.5
plot([TTR_NabPac468,TTR_NabPac468],[0,maxy(2)],'k--') %172
plot([timeplot(length(dataplot))*10,timeplot(length(dataplot))*10],[0,maxy(2)],'k--') %179 176
plot([timescale3(end),timescale3(end)],[0,maxy(2)],'k--') %164
hold off
ylim([0,450])
xlabel('days')
ylabel('tumor burden')
maxad_486NabPac=176/164;
%% 288 letrozole
M=[-0.199966609866911,0;0.000100000000000000,0.00490377420494712];
Moff=[0.0436,0;0.000100000000000000,0.00490377420494712];
n0=[117.820000000000;19.1800000000000];
x=1.6;
timeon=8.8;
dayplot=(tumormarker_data.relativedays2(51:103)-14975);
%[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
[breakadvantage,time,dataplot,timeplot,trelapse,singledosen] = plotOptimalBreak_delayedstart(M,Moff,n0,x,timeon,dayplot(2));
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'k');
hold on
%timescale2=linspace(0,time,length(dataplot))*10;
%plot(timescale2,dataplot(:,3),'m')
plot(timeplot*10,dataplot(:,3),'m')
scatter(dayplot,tumormarker_data.CA_15_3num(51:103),'r','filled')
maxy=ylim;
ind=timescale(find(singledosen<121));
plot([ind(end),ind(end)],[0,maxy(2)],'k--') %3750
plot([TTR_Let288,TTR_Let288],[0,maxy(2)],'k--') %4188
%ind2=timescale2(find(dataplot(:,3)<121));
ind2=timeplot(find(dataplot(:,3)<121));
plot([ind2(end)*10,ind2(end)*10],[0,maxy(2)],'k--') %4254 %4425
hold off
xlabel('days')
ylabel('tumor burden')
maxad_288Let=4425/3750;
%% 288 everolimus-exmestesane
M=[-0.0906726335969283,0;0.000100000000000000,0.0569641460063333];
Moff=[0.0334414764411500,0;0.000100000000000000,0.0569641460063333];
n0=[143.178963893250;41.8210361067504];
x=1.9925;
timeon=3.7275;
dayplot=(tumormarker_data.relativedays2([104,106:112])-19210);
%[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
[breakadvantage,time,dataplot,timeplot,trelapse,singledosen] = plotOptimalBreak_delayedstart(M,Moff,n0,x,timeon,dayplot(2));
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'k');
hold on
%timescale2=linspace(0,time,length(dataplot))*10;
%plot(timescale2,dataplot(:,3),'m')
plot(timeplot*10,dataplot(:,3),'m')
scatter(dayplot,tumormarker_data.CA_15_3num([104,106:112]),'r','filled')
ylim([50,250])
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') %245.5
plot([TTR_Everol288,TTR_Everol288],[0,maxy(2)],'k--') %252
plot([timeplot(end)*10,timeplot(end)*10],[0,maxy(2)],'k--') %333 %335
hold off
xlabel('days')
ylabel('tumor burden')
maxad_288Everol=335/245;
%% 15 Gemcitabine
M=[-9.57769093958680,0;0.000100000000000000,0.0510110906458706];
Moff=[0.703220893447389,0;0.000100000000000000,0.0510110906458706];
n0=[235.294117647059;164.705882352941];
x_in=round(6/2);
x=1.14;
timeon=1.03;
%[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,trelapse] = plotOptimalBreak_pulsevspulse(M,Moff,n0,x_in,x,1,6,timeon,0);
[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,timeplot,trelapse] = plotOptimalBreak_pulsevspulse_delayedstart(M,Moff,n0,x_in,x,1,6,timeon,0,2);
figure
timescale=linspace(0,trelapse,length(dataplot_in));
plot(timescale,dataplot_in(:,3),'m');
hold on
%timescale2=linspace(0,time,length(dataplot));
%plot(timescale2,dataplot(:,3),'m--')
plot(timeplot(1:length(dataplot))*10,dataplot(:,3),'m--')
timescale3=linspace(0,singledose_trelapse,length(singledosen));
plot(timescale3,singledosen,'k')
dayplot=(tumormarker_data.relativedays2(11:28)-16313);
scatter(dayplot,tumormarker_data.CA_15_3num(11:28),'r','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') %198
plot([TTR_Gem15,TTR_Gem15],[0,maxy(2)],'k--') %197
plot([timeplot(length(dataplot))*10,timeplot(length(dataplot))*10],[0,maxy(2)],'k--') %231 230
plot([timescale3(end),timescale3(end)],[0,maxy(2)],'k--') %174
hold off
ylim([0,450])
xlabel('days')
ylabel('tumor burden')
maxad_15Gem=230/174;
%% 422 eribulin
M=[-10.1776168548742,0;0.00100000000000000,0.0641910572531241];
Moff=[0.193072530736878,0;0.00100000000000000,0.0641910572531241];
n0=[200.489795918367;106.510204081633];
x_in=round(14/2);
x=1.1;
timeon=1;
%[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,trelapse] = plotOptimalBreak_pulsevspulse(M,Moff,n0,x_in,x,1,13,timeon,2);
[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,timeplot,trelapse] = plotOptimalBreak_pulsevspulse_delayedstart(M,Moff,n0,x_in,x,1,13,timeon,2,2);

figure
ind=find(dataplot_in(:,3)<252);
timescale=linspace(0,trelapse,ind(end));
plot(timescale,dataplot_in(1:ind(end),3),'m');
hold on
ind2=find(dataplot(:,3)<252);
%timescale2=linspace(0,time,ind2(end));
%plot(timescale2,dataplot(1:ind2(end),3),'m--')
plot(timeplot(1:length(dataplot))*10,dataplot(:,3),'m--')
ind3=find(singledosen<252);
timescale3=linspace(0,singledose_trelapse,ind3(end));
plot(timescale3,singledosen(1:ind3(end)),'k')
dayplot=(tumormarker_data.relativedays2(180:190)-19214);
scatter(dayplot,tumormarker_data.CA_15_3num(180:190),'r','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') %143
plot([TTR_Erib422,TTR_Erib422],[0,maxy(2)],'k--') %139
plot([timeplot(length(dataplot))*10,timeplot(length(dataplot))*10],[0,maxy(2)],'k--') %152 150
plot([timescale3(end),timescale3(end)],[0,maxy(2)],'k--') %134
hold off
xlabel('days')
ylabel('tumor burden')
maxad_422Erib=150/134;
%% 422 paclitaxel
M=[-20.2543174771295,0;0.00100000000000000,0.101082790799119];
Moff=[0.652969934663618,0;0.00100000000000000,0.101082790799119];
n0=[190.488188976378;61.5118110236220];
x_in=round(21/2);
x=1.08;
timeon=1.1;
%[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,trelapse] = plotOptimalBreak_pulsevspulse(M,Moff,n0,x_in,x,1,20,timeon,0);
[breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,timeplot,trelapse] = plotOptimalBreak_pulsevspulse_delayedstart(M,Moff,n0,x_in,x,1,20,timeon,0,2);

figure
timescale=linspace(0,trelapse,length(dataplot_in));
plot(timescale,dataplot_in(:,3),'m');
hold on
%timescale2=linspace(0,time,length(dataplot));
%plot(timescale2,dataplot(:,3),'m--')
plot(timeplot(1:length(dataplot))*10,dataplot(:,3),'m--')
timescale3=linspace(0,singledose_trelapse,length(singledosen));
plot(timescale3,singledosen,'k')
dayplot=(tumormarker_data.relativedays2(190:199)-19353);
scatter(dayplot,tumormarker_data.CA_15_3num(190:199),'r','filled')
ylim([50,250])
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') %145.5
plot([TTR_Pac422,TTR_Pac422],[0,maxy(2)],'k--') %137
plot([timeplot(length(dataplot))*10,timeplot(length(dataplot))*10],[0,maxy(2)],'k--') %151 149
plot([timescale3(end),timescale3(end)],[0,maxy(2)],'k--') %139.5
hold off
xlabel('days')
ylabel('tumor burden')
maxad_422pac=149/139.5;
%% 422 both
n0=[128.072230435481
85.2081632653061
59.9496062992125
33.7700000000000];
u=10^-3;

MA=zeros(4,4);
MA(4,4)=.0856; %M(4,4) from both Pac and erib, average
MA(2,2)=.064; %cells resistant to erib
MA(3,3)=-10.1776168548742; %M from Erib422
MA(1,1)=-15.2160; %average both Ms
MA(2,1)=u;
MA(3,1)=u;
MA(4,2)=u;
MA(4,3)=u;
MB=MA;
MB(2,2)=-20.2543174771295; %M from Pac422
MB(3,3)=.1; %cells resistant to pac
Moff=MA;
Moff(1,1)=.4; %Moff from both Pac and erib
Moff(2,2)=.06;
Moff(3,3)=.11;

Atime=5.896;
Btime=5.89;
x=9.1206;
timeon=2.08;
%[trelapse_in,dataplot_in_all,timeplot_in_all,breakadvantage_cont,dataplot_continuous,timeplot_continuous,breakadvantage_alternating,dataplot_alternating,timeplot_alternating,breakadvantage_combobreak,dataplot_combobreak,timeplot_combobreak] = plotOptimalBreak2drugs_pulsevspulse( MA,MB,Moff,n0,1,13,1,20,7,11,Atime,Btime,x,timeon);
[trelapse_in,dataplot_in_all,timeplot_in_all,breakadvantage_cont,dataplot_continuous,timeplot_continuous,breakadvantage_alternating,dataplot_alternating,timeplot_alternating,breakadvantage_combobreak,dataplot_combobreak,timeplot_combobreak] = plotOptimalBreak2drugs_pulsevspulse_delayedstart( MA,MB,Moff,n0,1,13,1,20,7,11,Atime,Btime,x,timeon,1);

figure
plot(timeplot_in_all*10,dataplot_in_all(:,5),'m');
hold on
plot(timeplot_continuous*10,dataplot_continuous(:,5),'k')
plot(timeplot_alternating*10,dataplot_alternating(:,5),'r')
plot(timeplot_combobreak*10,dataplot_combobreak(:,5),'b')
dayplot=(tumormarker_data.relativedays2(180:199)-19214);
scatter(dayplot,tumormarker_data.CA_15_3num(180:199),'r','filled')
maxy=ylim;
plot([timeplot_in_all(end)*10,timeplot_in_all(end)*10],[0,maxy(2)],'k--') % 270
plot([TTR_Pac422+TTR_Erib422,TTR_Pac422+TTR_Erib422],[0,maxy(2)],'k--') %276
plot([timeplot_continuous(end)*10,timeplot_continuous(end)*10],[0,maxy(2)],'k--') %255 
plot([timeplot_alternating(end)*10,timeplot_alternating(end)*10],[0,maxy(2)],'k--') %321 318
plot([timeplot_combobreak(end)*10,timeplot_combobreak(end)*10],[0,maxy(2)],'k--') %315 311
hold off
xlabel('days')
ylabel('tumor burden')
%% all stats
Ratio_468Let=1;
Ratio_468NabPac=1.3;
Ratio_288Let=.93;
Ratio_288Exem=.9;
Ratio_15Gem=1.13;
Ratio_422Erib=.9;
Ratio_422Pac=.99;

maxad_468Let=1.6;
maxad_468NabPac=1.07;
maxad_288Let=1.18;
maxad_288Exem=1.4;
maxad_15Gem=1.32;
maxad_422Erib=1.12;
maxad_422Pac=1.07;

relapse_ratioall=[Ratio_468Let,Ratio_468NabPac,Ratio_288Let,Ratio_288Exem,Ratio_15Gem,Ratio_422Erib,Ratio_422Pac];
relapse_ratioall_mean=mean(relapse_ratioall); %1.02
relapse_ratioall_meandiff=mean(abs(relapse_ratioall-1)); 
maxadvantage_all=[maxad_468Let,maxad_468NabPac,maxad_288Let,maxad_288Exem,maxad_15Gem,maxad_422Erib,maxad_422Pac];
maxadvantage_all_mean=mean(maxadvantage_all); %1.25

%plot as overlapping bars
figure
bar([1:7],maxadvantage_all,.5,'FaceColor',[0.2 0.2 0.5])
hold on
bar([1:7],relapse_ratioall,.25,'FaceColor',[0 0.7 0.7])
plot([0,8],[relapse_ratioall_mean,relapse_ratioall_mean],'--','Color',[0 0.7 0.7])
plot([0,8],[maxadvantage_all_mean,maxadvantage_all_mean],'--','Color',[.2 0.2 0.5])
hold off
