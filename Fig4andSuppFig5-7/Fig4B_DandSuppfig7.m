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
timeon=(3.9669-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2(205:213)-19939);
scatter(dayplot,tumormarker_data.CA_15_3num(205:213),'k','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') 
plot([TTR_Let468,TTR_Let468],[0,maxy(2)],'k--') 
plot([timescale2(end),timescale2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 468 NabPac
M=[-0.697257232685537,0;0.00100000000000000,0.0613309909961578];
Moff=[0.431081420637815,0;0.00100000000000000,0.0613309909961578];
n0=[388.842038724327;60.1579612756727];
x=2.3;
timeon=(2.16-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2(216:222)-20337);
scatter(dayplot,tumormarker_data.CA_15_3num(216:222),'k','filled')
maxy=ylim;
ind=timescale(find(singledosen<256));
plot([ind(end),ind(end)],[0,maxy(2)],'k--') 
plot([TTR_NabPac468,TTR_NabPac468],[0,maxy(2)],'k--') 
ind2=timescale2(find(dataplot(:,3)<256));
plot([ind2(end),ind2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 288 letrozole
M=[-0.156842898604616,0;0.000100000000000000,0.0188237392714033];
Moff=[0.0397952017353970,0;0.000100000000000000,0.0188237392714033];
M=[-0.199966609866911,0;0.000100000000000000,0.00490377420494712];
Moff=[0.0436,0;0.000100000000000000,0.00490377420494712];
n0=[117.820000000000;19.1800000000000];
x=2.97;
x=1.6;
timeon=(100-1)/20;
timeon=(8.8-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2(51:103)-14975);
scatter(dayplot,tumormarker_data.CA_15_3num(51:103),'k','filled')
maxy=ylim;
ind=timescale(find(singledosen<121));
plot([ind(end),ind(end)],[0,maxy(2)],'k--') 
plot([TTR_Let288,TTR_Let288],[0,maxy(2)],'k--') 
ind2=timescale2(find(dataplot(:,3)<121));
plot([ind2(end),ind2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 288 everolimus-exmestesane
M=[-0.0906726335969283,0;0.000100000000000000,0.0569641460063333];
Moff=[0.0334414764411500,0;0.000100000000000000,0.0569641460063333];
n0=[143.178963893250;41.8210361067504];
x=1.9925;
timeon=(3.7275-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2([104,106:112])-19210);
scatter(dayplot,tumormarker_data.CA_15_3num([104,106:112]),'k','filled')
ylim([50,250])
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') 
plot([TTR_Everol288,TTR_Everol288],[0,maxy(2)],'k--') 
plot([timescale2(end),timescale2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 15 Gemcitabine
M=[-0.0506759862375532,0;0.000100000000000000,0.0328193057075717];
Moff=[0.778846717125000,0;0.000100000000000000,0.0328193057075717];
n0=[235.294117647059;164.705882352941];
x=2.19;
timeon=(3.19-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2(11:28)-16313);
scatter(dayplot,tumormarker_data.CA_15_3num(11:28),'k','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') 
plot([TTR_Gem15,TTR_Gem15],[0,maxy(2)],'k--') 
plot([timescale2(end),timescale2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 422 eribulin
M=[-0.854375066836716,0;0.00100000000000000,0.0406733029715053];
Moff=[0.325414311189373,0;0.00100000000000000,0.0406733029715053];
n0=[232.062992125984;74.9370078740156];
n0=[200.489795918367;106.510204081633];
x=2.2131;
x=1.5;
timeon=(2.3414-1)/20;
timeon=(3.78-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2(180:190)-19214);
scatter(dayplot,tumormarker_data.CA_15_3num(180:190),'k','filled')
maxy=ylim;
ind=timescale(find(singledosen<252));
plot([ind(end),ind(end)],[0,maxy(2)],'k--') 
plot([TTR_Erib422,TTR_Erib422],[0,maxy(2)],'k--') 
ind2=timescale2(find(dataplot(:,3)<252));
plot([ind2(end),ind2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 422 paclitaxel
M=[-0.877367588218591,0;0.00100000000000000,0.105170755980731];
Moff=[0.276057896874889,0;0.00100000000000000,0.105170755980731];
n0=[190.488188976378;61.5118110236220];
x=2.9;
timeon=(2.8-1)/20;
[breakadvantage,time,dataplot,trelapse,singledosen] = plotOptimalBreak(M,Moff,n0,x,timeon);
timescale=linspace(0,trelapse,length(singledosen))*10;
figure,plot(timescale,singledosen,'b');
hold on
timescale2=linspace(0,time,length(dataplot))*10;
plot(timescale2,dataplot(:,3),'r')
dayplot=(tumormarker_data.relativedays2(190:199)-19353);
scatter(dayplot,tumormarker_data.CA_15_3num(190:199),'k','filled')
ylim([50,250])
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') 
plot([TTR_Pac422,TTR_Pac422],[0,maxy(2)],'k--') 
plot([timescale2(end),timescale2(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% 422 both
n0=[126.787596014784
35.5034013605442
24.9790026246719
119.730000000000];
u=10^-3;

MA=zeros(4,4);
MA(4,4)=.07; %M(4,4) from both Pac and erib, average
MA(2,2)=.16; %average of (cells sensitive to pac but off pac treatment, and cells resistant to erib)
MA(3,3)=-0.85; %M from Erib422
MA(1,1)=-0.86; %average both Ms
MA(2,1)=u;
MA(3,1)=u;
MA(4,2)=u;
MA(4,3)=u;
MB=MA;
MB(2,2)=-0.88; %M from Pac422
MB(3,3)=.2;%average of (cells sensitive to erib but off erib treatment, and cells resistant to pac)
Moff=MA;
Moff(1,1)=.3; %Moff from both Pac and erib
Moff(2,2)=.2; %assuming single resistant same growth rate as double
Moff(3,3)=.3;

Atime=(6.8-1)/20;
Btime=(4.7-1)/20;
x=2.57;
timeon=(2-1)/20;
[trelapseAthenB,AthenBn,advantagealternating2,trelapsealternating2,dataplot_alt,advantagecombobreak,dataplot_combo,trelapsecombobreak] = plotOptimalBreak2drugs( MA,MB,Moff,n0,Atime,Btime,x,timeon);

timescale=linspace(0,trelapseAthenB,length(AthenBn))*20;
figure,plot(timescale,AthenBn,'k');
hold on
timescale2=linspace(0,trelapsealternating2,length(dataplot_alt))*20;
plot(timescale2,dataplot_alt(:,5),'m')
timescale3=linspace(0,trelapsecombobreak,length(dataplot_combo))*20;
plot(timescale3,dataplot_combo(:,5),'b')
dayplot=(tumormarker_data.relativedays2(180:199)-19214);
scatter(dayplot,tumormarker_data.CA_15_3num(180:199),'k','filled')
maxy=ylim;
plot([timescale(end),timescale(end)],[0,maxy(2)],'k--') 
plot([TTR_Pac422+TTR_Erib422,TTR_Pac422+TTR_Erib422],[0,maxy(2)],'k--') 
plot([timescale2(end),timescale2(end)],[0,maxy(2)],'k--') 
plot([timescale3(end),timescale3(end)],[0,maxy(2)],'k--') 
hold off
xlabel('days')
ylabel('tumor burden')
%% all stats
Ratio_468Let=1;
Ratio_468NabPac=1.7;
Ratio_288Let=.93;
Ratio_288Everol=.98;
Ratio_288Exem=.9;
Ratio_288EverolExem=1.01;
Ratio_15Gem=.7;
Ratio_422Erib=.9;
Ratio_422Pac=1.4;

maxad_468Let=1.6;
maxad_468NabPac=1.3;
maxad_288Let=1.12;
maxad_288Everol=1.2;
maxad_288Exem=1.4;
maxad_288EverolExem=1.4;
maxad_15Gem=1.6;
maxad_422Erib=1.5;
maxad_422Pac=1.3;

relapse_ratioall=[Ratio_468Let,Ratio_468NabPac,Ratio_288Let,Ratio_288Exem,Ratio_15Gem,Ratio_422Erib,Ratio_422Pac];
relapse_ratioall_mean=mean(relapse_ratioall); 
relapse_ratioall_meandiff=mean(abs(relapse_ratioall-1)); 
maxadvantage_all=[maxad_468Let,maxad_468NabPac,maxad_288Let,maxad_288Exem,maxad_15Gem,maxad_422Erib,maxad_422Pac];
maxadvantage_all_mean=mean(maxadvantage_all); 

%plot as overlapping bars
figure
bar([1:7],maxadvantage_all,.5,'FaceColor',[0.2 0.2 0.5])
hold on
bar([1:7],relapse_ratioall,.25,'FaceColor',[0 0.7 0.7])
plot([0,8],[relapse_ratioall_mean,relapse_ratioall_mean],'--','Color',[0 0.7 0.7])
plot([0,8],[maxadvantage_all_mean,maxadvantage_all_mean],'--','Color',[.2 0.2 0.5])
hold off