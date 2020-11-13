%% combo pulse example
Moff =[0.7699         0         0         0
    0.0007    0.6374         0         0
    0.0009         0    0.6486         0
         0    0.0001    0.0003    0.1148];
MA =[ -0.7623         0         0         0
    0.0010    0.6374         0         0
    0.0009         0   -0.8837         0
         0    0.0001    0.0008    0.1148];
MB =[-0.5935         0         0         0
    0.0007   -0.7261         0         0
    0.0009         0    0.6486         0
         0    0.0009    0.0003    0.1148];
n0=[58;29;10;3];
tspan=[0:.05:500];
fA=Moff(1,1)-MA(1,1); %drug effect of each drug
fB=Moff(1,1)-MB(1,1);

%standard of care conditions: each drug individually or sequential treatment
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA);
[t,n] = ode45(dn, tspan, n0);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseA=min(t(n(:,5)>sum(n0))); %time to relapse if drugA was given continuously
singledosenA = n(1:round(trelapseA*20),5);
singledose_trelapseA=trelapseA;
n0Anew=n(round(trelapseA*20),1:4)';

dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB);
[t,n] = ode45(dn, tspan, n0);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseB=min(t(n(:,5)>sum(n0))); %time to relapse if drugB was given continuously (would be equal to time to relapse under drug A in the simple case that both drugs are equally effective)
singledosenB = n(1:round(trelapseB*20),5);
singledose_trelapseB=trelapseB;
n0Bnew=n(round(trelapseB*20),1:4)';

%time to relapse if drug A was maxed out and then drug B was given until
%relapse
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB);
[t,n] = ode45(dn, tspan, n0Anew);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseB2=min(t(n(:,5)>sum(n0)));
if isempty(trelapseB2)
    trelapseB2=0;
end
trelapseAthenB=trelapseB2+singledose_trelapseA;
AthenBn=[singledosenA;n(1:round(trelapseB2*20),5)];

%time to relapse if drug B was maxed out and then drug A was given until

dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA);
[t,n] = ode45(dn, tspan, n0Bnew);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseA2=min(t(n(:,5)>sum(n0)));
if isempty(trelapseA2)
    trelapseA2=0;
end
trelapseBthenA=trelapseA2+singledose_trelapseB;
BthenAn=[singledosenB;n(1:round(trelapseA2*20),5)];
BthenAt=t(1:length(BthenAn));

bestSOC=max(trelapseAthenB,trelapseBthenA);
if trelapseAthenB>trelapseBthenA
    Afirst=1;
else
    Afirst=0;
end

x = 5.0842;
%maxtimeon=trelapsecombo/x; 
%timeon=.05 + (maxtimeon)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeon=0.5475;
timeoff=timeon*(x-1);

if fA>fB
    dosecombo=fA;
    doseB=(fB/(fA+fB))*dosecombo; %weighted proportional to drug effect 
    doseA=fA-doseB;
else
    dosecombo=fB;
    doseA=(fA/(fA+fB))*dosecombo;
    doseB=fB-doseA;
end


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
        timeplot=[timeplot;t+timeplot(end)];
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
            timeplot=[timeplot;t+timeplot(end)];
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
    timeplot=[timeplot;t+timeplot(end)];
    iteration=iteration+1;
    
end
trelapsecombobreak=totaltimeon+totaltimeoff;
advantagecombobreak=trelapsecombobreak/bestSOC; %if greater than 1, break strategy is better
figure, plot(timeplot,dataplot(:,5))
hold on
plot(BthenAt,BthenAn)
hold off

%% alternating example

Moff =[0.8164         0         0         0
    0.0008    0.4898         0         0
    0.0010         0    0.3799         0
         0    0.0002    0.0006    0.0417];
MA =[-0.4585         0         0         0
    0.0005    0.4898         0         0
    0.0004         0   -0.8950         0
         0    0.0001    0.0007    0.0417];
MB =[-0.6172         0         0         0
    0.0004   -0.9439         0         0
    0.0010         0    0.3799         0
         0    0.0007    0.0005    0.0417];
n0= [48
    17
    20
    15];

tspan=[0:.05:500];

%standard of care conditions: each drug individually or sequential treatment
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA);
[t,n] = ode45(dn, tspan, n0);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseA=min(t(n(:,5)>sum(n0))); %time to relapse if drugA was given continuously
singledosenA = n(1:round(trelapseA*20),5);
singledose_trelapseA=trelapseA;
n0Anew=n(round(trelapseA*20),1:4)';

dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB);
[t,n] = ode45(dn, tspan, n0);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseB=min(t(n(:,5)>sum(n0))); %time to relapse if drugB was given continuously (would be equal to time to relapse under drug A in the simple case that both drugs are equally effective)
singledosenB = n(1:round(trelapseB*20),5);
singledose_trelapseB=trelapseB;
n0Bnew=n(round(trelapseB*20),1:4)';

%time to relapse if drug A was maxed out and then drug B was given until
%relapse
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB);
[t,n] = ode45(dn, tspan, n0Anew);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseB2=min(t(n(:,5)>sum(n0)));
if isempty(trelapseB2)
    trelapseB2=0;
end
trelapseAthenB=trelapseB2+singledose_trelapseA;
AthenBn=[singledosenA;n(1:round(trelapseB2*20),5)];

%time to relapse if drug B was maxed out and then drug A was given until

dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA);
[t,n] = ode45(dn, tspan, n0Bnew);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseA2=min(t(n(:,5)>sum(n0)));
if isempty(trelapseA2)
    trelapseA2=0;
end
trelapseBthenA=trelapseA2+singledose_trelapseB;
BthenAn=[singledosenB;n(1:round(trelapseA2*20),5)];
BthenAt=t(1:length(BthenAn));

bestSOC=max(trelapseAthenB,trelapseBthenA);
if trelapseAthenB>trelapseBthenA
    Afirst=1;
else
    Afirst=0;
end

% Alternating: no breaks 
Atime=1.7;
Btime=1.5;
%Atime=.05 + (trelapseA)*rand(1); %random amount of time to give drug A, as long as it's less than the time it was given continuously
%Btime=.05 + (trelapseB)*rand(1); 

n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeonA=0;
totaltimeonB=0;
trelapse=0;
dataplot=[];
timeplot=[];
tspanonA=[0:.05:Atime];
tspanonB=[0:.05:Btime];
iteration=0;

if Afirst
while 1
    if sum(n0new)>sum(n0)  %if tumor exceeds starting size 
        break
    end
    
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA );
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
         timeplot=[timeplot;t+timeplot(end)];
    end
        iteration=iteration+1;
        
     
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB );
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
        timeplot=[timeplot;t+timeplot(end)];
    iteration=iteration+1;
end

else
    while 1
    if sum(n0new)>sum(n0)  %if tumor exceeds starting size 
        break
    end
    
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MB );
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
    
    dataplot=[dataplot;n];
    if iteration==0
    timeplot=[timeplot;t];
    else
         timeplot=[timeplot;t+timeplot(end)];
    end
        iteration=iteration+1;
        
     
    dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA );
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
        timeplot=[timeplot;t+timeplot(end)];

    iteration=iteration+1;
    end
end
trelapsealternating2=totaltimeon;
figure, plot(timeplot,dataplot(:,5))
hold on
plot(BthenAt,BthenAn)
hold off

advantagealternating2=trelapsealternating2/bestSOC;

