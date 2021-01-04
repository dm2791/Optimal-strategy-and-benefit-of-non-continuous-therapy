function [trelapse_in,dataplot_in_all,timeplot_in_all,breakadvantage_cont,dataplot_continuous,timeplot_continuous,breakadvantage_alternating,dataplot_alternating,timeplot_alternating,breakadvantage_combobreak,dataplot_combobreak,timeplot_combobreak] = plotOptimalBreak2drugs_pulsevspulse( MA,MB,Moff,n0,timeon_in_erib,timeoff_in_erib,timeon_in_pac,timeoff_in_pac,x_in_erib,x_in_pac,Atime,Btime,x,timeon)
%% baseline: pulse MA and then pulse MB
tspanon=[0:.05:timeon_in_erib/10];
tspanoff=[0:.05:timeoff_in_erib/10];
n0new=n0;
dataplot_in=[];
timeplot_in=[];
for j=1:9
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA );
    [t,n] = ode45(dn, tspanon, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    dataplot_in=[dataplot_in;n];
    if j==1
        timeplot_in=[timeplot_in;t];
    else
        timeplot_in=[timeplot_in;t+.05+timeplot_in(end)];
    end
    n0new=[n(end,1:4)'];
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff );
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    dataplot_in=[dataplot_in;n];
    timeplot_in=[timeplot_in;t+.05+timeplot_in(end)];
    n0new=[n(end,1:4)'];
end
dataplot_in_erib=dataplot_in;
timeplot_in_erib=timeplot_in;

tspanon=[0:.05:timeon_in_pac/10];
tspanoff=[0:.05:timeoff_in_pac/10];
dataplot_in=[];
timeplot_in=[];
for j=1:8
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB );
    [t,n] = ode45(dn, tspanon, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    dataplot_in=[dataplot_in;n];
    if j==1
        timeplot_in=[timeplot_in;t];
    else
        timeplot_in=[timeplot_in;t+.05+timeplot_in(end)];
    end
    n0new=[n(end,1:4)'];
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff );
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    dataplot_in=[dataplot_in;n];
    timeplot_in=[timeplot_in;t+.05+timeplot_in(end)];
    n0new=[n(end,1:4)'];
end
dataplot_in_pac=dataplot_in;
timeplot_in_pac=timeplot_in;
dataplot_in_all=[dataplot_in_erib;dataplot_in_pac];
timeplot_in_all=[timeplot_in_erib;timeplot_in_erib(end)+timeplot_in_pac];
dataplot_in_all=dataplot_in_all(1:max(find(dataplot_in_all(:,5)<sum(n0)))+1,:);
timeplot_in_all=timeplot_in_all(1:max(find(dataplot_in_all(:,5)<sum(n0)))+1);
trelapse_in=timeplot_in_all(end)*10;

%% continuous
fA=Moff(1,1)-MA(1,1); %drug effect
f_baseA=fA/x_in_erib;

tspan=[0:.05:500];
M_baseA=Moff;
M_baseA(1,1)=Moff(1,1)-f_baseA;
M_baseA(3,3)=Moff(3,3)-f_baseA;
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,M_baseA );
[t,n] = ode45(dn, tspan, n0);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
[c ind]=min(n(:,5));
trelapse=t(min(find(n(ind:end,5)>252))+ind-1)*10;
singledosenA = n(1:round(trelapse*2),5);
singledose_trelapseA=trelapse;
n0new=[n(min(find(n(ind:end,5)>252))+ind-1,1:4)'];
dataplot_contA=n(1:min(find(n(ind:end,5)>252))+ind-1,:);
timeplot_contA=t(1:min(find(n(ind:end,5)>252))+ind-1,:);

fB=Moff(1,1)-MB(1,1); %drug effect
f_baseB=fB/x_in_pac;

tspan=[0:.05:500];
M_baseB=Moff;
M_baseB(1,1)=Moff(1,1)-f_baseB;
M_baseB(2,2)=Moff(2,2)-f_baseB;
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,M_baseB );
[t,n] = ode45(dn, tspan, n0new);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
[c ind]=min(n(:,5));
trelapse=t(min(find(n(ind:end,5)>sum(n0)))+ind-1)*10;
singledosenB = n(1:round(trelapse*2),5);
singledose_trelapseB=trelapse;
dataplot_contB=n(1:min(find(n(ind:end,5)>sum(n0)))+ind-1,:);
timeplot_contB=t(1:min(find(n(ind:end,5)>sum(n0)))+ind-1,:);

singledose_trelapse=singledose_trelapseA+singledose_trelapseB;
dataplot_continuous=[dataplot_contA;dataplot_contB];
timeplot_continuous=[timeplot_contA;timeplot_contB+.05+timeplot_contA(end)];
breakadvantage_cont=trelapse_in/singledose_trelapse;

%% alternating
n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeonA=0;
totaltimeonB=0;
trelapse=0;
dataplot=[];
timeplot=[];
Atime=(Atime-1)/20;
Btime=(Btime-1)/20;
tspanonA=[0:.05:Atime];
tspanonB=[0:.05:Btime];
iteration=0;
while 1
    if sum(n0new)>sum(n0)  %if tumor exceeds starting size 
        break
    end
    
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,M_baseA );
    [t,n] = ode45(dn, tspanonA, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    n0new=[n(end,1:4)'];  %calculate n at end of on interval; will be new n0 for start of off interval

    if n(end,5)>sum(n0)
                stoppoint=min(t(n(:,5)>sum(n0)));
        if round(stoppoint*20)<=1
            break
        end
        n0new=[n(round(stoppoint*20),1:4)'];
        n=n(1:round(stoppoint*20),:);
        t=t(1:round(stoppoint*20),:);
        totaltimeonA=totaltimeonA+stoppoint;
        totaltimeon=totaltimeon+stoppoint;
    else
    totaltimeonA=totaltimeonA+Atime;
    totaltimeon=totaltimeon+Atime;
    end
    
    dataplot=[dataplot;n];
    if iteration==0
    timeplot=[timeplot;t];
    else
         timeplot=[timeplot;t+.05+timeplot(end)];
    end
        iteration=iteration+1;
        
     
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,M_baseB );
    [t,n] = ode45(dn, tspanonB, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    n0new=[n(end,1:4)'];  %calculate n at end of on interval; will be new n0 for start of off interval
    if n(end,5)>sum(n0)
                stoppoint=min(t(n(:,5)>sum(n0)));
        if round(stoppoint*20)<=1
            break
        end
        n0new=[n(round(stoppoint*20),1:4)'];
        n=n(1:round(stoppoint*20),:);
        t=t(1:round(stoppoint*20),:);
        totaltimeonB=totaltimeonB+stoppoint;
        totaltimeon=totaltimeon+stoppoint;
    else
    totaltimeonB=totaltimeonB+Btime;
    totaltimeon=totaltimeon+Btime;
end
    
        totaln=n(:,5);
    dataplot=[dataplot;n];
        timeplot=[timeplot;t+.05+timeplot(end)];
    iteration=iteration+1;
end
dataplot_alternating=dataplot;
timeplot_alternating=timeplot;
trelapse_alternating=timeplot(end)*10;
breakadvantage_alternating=trelapse_alternating/trelapse_in;

%% Combination with breaks
dosecombo=f_baseA;
doseB=(f_baseB/(f_baseA+f_baseB))*dosecombo; %weighted proportional to drug effect 
doseA=f_baseA-doseB;
timeon=(timeon-1)/20;
timeoff=timeon*(x-1);

Mnew=Moff;
Mnew(1,1)=Moff(1,1)-x*doseA-x*doseB;  %make new drug matrix with higher dose
Mnew(2,2)=Moff(2,2)-x*doseB;
Mnew(3,3)=Moff(3,3)-x*doseA;

n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeoff=0;
iteration=0;
dataplot=[];
timeplot=[];
tspanon=[0:.05:timeon];
tspanoff=[0:.05:timeoff];
if length(tspanon)<2 || length(tspanoff)<2
    advantagecombobreak=0;
    totaltimeon=0;
    totaltimeoff=0;
    return
end

while 1
    if sum(n0new)>sum(n0)  %if tumor exceeds starting size
        break
    end
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Mnew);
    [t,n] = ode45(dn, tspanon, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    n0new=[n(end,1:4)'];  %calculate n at end of on interval; will be new n0 for start of off interval
    
    if n(end,5)>sum(n0)
        stoppoint=min(t(n(:,5)>sum(n0)));
        if round(stoppoint*20)<=1
            break
        end
        n0new=[n(round(stoppoint*20),1:4)'];
        n=n(1:round(stoppoint*20),:);
        t=t(1:round(stoppoint*20),:);
        add=stoppoint;
        totaltimeon=totaltimeon+stoppoint;
    else
        add=timeon;
        totaltimeon=totaltimeon+timeon;
    end
    
    dataplot=[dataplot;n];
    if iteration==0
        timeplot=[timeplot;t];
    else
        timeplot=[timeplot;t+.05+timeplot(end)];
    end
    iteration=iteration+1;
    
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff);
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
    n0new=[n(end,1:4)'];
    if n(end,5)>sum(n0)
        stoppoint=min(t(n(:,5)>sum(n0)));
        if round(stoppoint*20)<=1 || stoppoint<add*(x-1)
            totaltimeoff=totaltimeoff+stoppoint; %THIS IS NEW
            n=n(1:round(stoppoint*20),:);
            dataplot=[dataplot;n];
            t=t(1:round(stoppoint*20),:);
            timeplot=[timeplot;t+.05+timeplot(end)];
            iteration=iteration+1;
            break
        end
        n0new=[n(round(stoppoint*20),1:4)'];
        n=n(1:round(stoppoint*20),:);
        t=t(1:round(stoppoint*20),:);
        totaltimeoff=totaltimeoff+stoppoint;
    else
        totaltimeoff=totaltimeoff+timeoff;
    end
    
    dataplot=[dataplot;n];
    timeplot=[timeplot;t+.05+timeplot(end)];
    iteration=iteration+1;
    
end
trelapsecombobreak=timeplot(end)*10;
dataplot_combobreak=dataplot;
timeplot_combobreak=timeplot;
breakadvantage_combobreak=trelapsecombobreak/trelapse_in; %if greater than 1, break strategy is better

