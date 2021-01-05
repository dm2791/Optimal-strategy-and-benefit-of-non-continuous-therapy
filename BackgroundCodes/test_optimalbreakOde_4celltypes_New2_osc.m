function [advantagecombo,advantagealternating2,advantageoscillating,advantagealternatingbreak,advantagecombobreak,Afirst,Atime,Btime,xA,timeonA,timeonB,x,timeon] = test_optimalbreakOde_4celltypes_New2_osc( MA,MB,Moff,n0)
%%  Ode model
% 4 cell types: sensitive to both drugs, resistant to drug A, resistant to
% drug B, double resistant
%  Optimize duration of each drug to extend time to relapse. No constraints
%  on total amount of dose, only ratio of on/off
tspan=[0:.05:500];

fA=Moff(1,1)-MA(1,1); %drug effect of each drug
fB=Moff(1,1)-MB(1,1);
%% standard of care conditions: each drug individually or sequential treatment
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
%relapse
dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MA);
[t,n] = ode45(dn, tspan, n0Bnew);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapseA2=min(t(n(:,5)>sum(n0)));
if isempty(trelapseA2)
    trelapseA2=0;
end
trelapseBthenA=trelapseA2+singledose_trelapseB;
BthenAn=[singledosenB;n(1:round(trelapseA2*20),5)];

bestSOC=max(trelapseAthenB,trelapseBthenA);
if trelapseAthenB>trelapseBthenA
    Afirst=1;
else
    Afirst=0;
end

%% Combo therapy
%assumption: drug must be equal to the greater value between fA and fB, and it needs to be weighted so
%that the more potent drug is given proportionally more

if fA>fB
    dosecombo=fA;
    doseB=(fB/(fA+fB))*dosecombo; %weighted proportional to drug effect
    doseA=fA-doseB;
else
    dosecombo=fB;
    doseA=(fA/(fA+fB))*dosecombo;
    doseB=fB-doseA;
end


Mon=Moff;
Mon(1,1)=Moff(1,1)-doseA-doseB;
Mon(2,2)=Moff(2,2)-doseB;
Mon(3,3)=Moff(3,3)-doseA;

dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Mon);
[t,n] = ode45(dn, tspan, n0);
n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
trelapsecombo=min(t(n(:,5)>sum(n0))); %time to relapse if combo was given continuously
combon = n(1:round(trelapsecombo*20),5);
combo_trelapse=trelapsecombo;
n0combonew=n(round(trelapsecombo*20),1:4)';

advantagecombo=combo_trelapse/bestSOC;

%% Alternating 2: no breaks
%keep alternating until relapse is reached, no limit on time since continuous treatment

Atime=.05 + (trelapseA)*rand(1); %random amount of time to give drug A, as long as it's less than the time it was given continuously
Btime=.05 + (trelapseB)*rand(1);

n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeonA=0;
totaltimeonB=0;
trelapse=0;
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
        
        
        
        iteration=iteration+1;
    end
end
trelapsealternating2=totaltimeon;

advantagealternating2=trelapsealternating2/bestSOC;

%% Oscillation (not used in data analysis)
n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeonA=0;
totaltimeonB=0;
totaltimeoff=0;
trelapse=0;
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
        
        
        iteration=iteration+1;
        
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff);
        [t,n] = ode45(dn, tspan, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        trelapse=min(t(n(:,5)>sum(n0)));  %calculate time to reach relapse size
        if trelapse==0
            break
        end
        totaltimeoff=totaltimeoff+trelapse;

        
        n0new=n(round(trelapse*20),1:4)'; %calculate n at end of off interval, will be new n0 for start of next on interval
        
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
        
  
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff);
        [t,n] = ode45(dn, tspan, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        trelapse=min(t(n(:,5)>sum(n0)));  %calculate time to reach relapse size
        if trelapse==0
            break
        end
        totaltimeoff=totaltimeoff+trelapse;
        
        n0new=n(round(trelapse*20),1:4)'; %calculate n at end of off interval, will be new n0 for start of next on interval
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
        
      
        iteration=iteration+1;
        
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff);
        [t,n] = ode45(dn, tspan, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        trelapse=min(t(n(:,5)>sum(n0)));  %calculate time to reach relapse size
        if trelapse==0
            break
        end
        totaltimeoff=totaltimeoff+trelapse;

        
        n0new=n(round(trelapse*20),1:4)'; %calculate n at end of off interval, will be new n0 for start of next on interval
        
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

        
        iteration=iteration+1;
        
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff);
        [t,n] = ode45(dn, tspan, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        trelapse=min(t(n(:,5)>sum(n0)));  %calculate time to reach relapse size
        if trelapse==0
            break
        end
        totaltimeoff=totaltimeoff+trelapse;

        
        n0new=n(round(trelapse*20),1:4)'; %calculate n at end of off interval, will be new n0 for start of next on interval
    end
end
trelapseoscillation=totaltimeon+totaltimeoff;
advantageoscillating=trelapseoscillation/bestSOC;
p=1;

%% Alternating with breaks

xA = 1 + (9)*rand(1); %random dose escalation up to 10 times original dose
maxtimeonA=trelapseA/xA; %maximum time this dose can be given is 1/x amount of time
timeonA=.05 + (maxtimeonA)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoffA=timeonA*(xA-1);

xB=xA;
maxtimeonB=trelapseB/xB; %maximum time this dose can be given is 1/x amount of time
timeonB=.05 + (maxtimeonB)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoffB=timeonB*(xB-1);

MnewA=MA;
MnewA(1,1)=Moff(1,1)-xA*fA;  %make new drug matrix with higher dose, affecting only those lines that are sensitive to A
MnewA(3,3)=Moff(3,3)-xA*fA;
MnewB=MB;
MnewB(1,1)=Moff(1,1)-xB*fB;  %make new drug matrix with higher dose, affecting only those lines that are sensitive to B
MnewB(2,2)=Moff(2,2)-xB*fB;

n0new=n0; %initialize n0new
totaltimeonA=0;
totaltimeonB=0;
totaltimeoffA=0;
totaltimeoffB=0;
iteration=0;
tspanonA=[0:.05:timeonA];
tspanonB=[0:.05:timeonB];
tspanoffA=[0:.05:timeoffA];
tspanoffB=[0:.05:timeoffB];
if length(tspanonA)<2 || length(tspanoffA)<2 || length(tspanonB)<2 || length(tspanoffB)<2
    advantagealternatingbreak=0;
    totaltimeonA=0;
    totaltimeoffA=0;
    totaltimeonB=0;
    totaltimeoffB=0;
    advantagecombobreak=0;
    return
end

if Afirst
    while 1
        if sum(n0new)>sum(n0)  %if tumor exceeds starting size
            break
        end
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MnewA);
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
            add=stoppoint;
            totaltimeonA=totaltimeonA+stoppoint;
        else
            add=timeonA;
            totaltimeonA=totaltimeonA+timeonA;
        end
        
        
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff );
        [t,n] = ode45(dn, tspanoffA, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        n0new=[n(end,1:4)'];
        if n(end,5)>sum(n0)
            stoppoint=min(t(n(:,5)>sum(n0)));
            if round(stoppoint*20)<=1 || stoppoint<add*(xA-1)
                totaltimeoffA=totaltimeoffA+stoppoint; %THIS IS NEW
                n=n(1:round(stoppoint*20),:);
                t=t(1:round(stoppoint*20),:);
                iteration=iteration+1;
                break
            end
            n0new=[n(round(stoppoint*20),1:2)'];
            n=n(1:round(stoppoint*20),:);
            t=t(1:round(stoppoint*20),:);
            totaltimeoffA=totaltimeoffA+stoppoint;
        else
            totaltimeoffA=totaltimeoffA+timeoffA;
        end
        
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MnewB);
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
            add=stoppoint;
            totaltimeonB=totaltimeonB+stoppoint;
        else
            add=timeonB;
            totaltimeonB=totaltimeonB+timeonB;
        end
        
        
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff );
        [t,n] = ode45(dn, tspanoffB, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        n0new=[n(end,1:4)'];
        if n(end,5)>sum(n0)
            stoppoint=min(t(n(:,5)>sum(n0)));
            if round(stoppoint*20)<=1 || stoppoint<add*(xB-1)
                totaltimeoffB=totaltimeoffB+stoppoint; %THIS IS NEW
                n=n(1:round(stoppoint*20),:);
                t=t(1:round(stoppoint*20),:);
                iteration=iteration+1;
                break
            end
            n0new=[n(round(stoppoint*20),1:2)'];
            n=n(1:round(stoppoint*20),:);
            t=t(1:round(stoppoint*20),:);
            totaltimeoffB=totaltimeoffB+stoppoint;
        else
            totaltimeoffB=totaltimeoffB+timeoffB;
        end
        iteration=iteration+1;
        
    end
else
    while 1
        if sum(n0new)>sum(n0)  %if tumor exceeds starting size
            break
        end
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MnewB);
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
            add=stoppoint;
            totaltimeonB=totaltimeonB+stoppoint;
        else
            add=timeonB;
            totaltimeonB=totaltimeonB+timeonB;
        end
        
        
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff );
        [t,n] = ode45(dn, tspanoffB, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        n0new=[n(end,1:4)'];
        if n(end,5)>sum(n0)
            stoppoint=min(t(n(:,5)>sum(n0)));
            if round(stoppoint*20)<=1 || stoppoint<add*(xB-1)
                totaltimeoffB=totaltimeoffB+stoppoint; %THIS IS NEW
                n=n(1:round(stoppoint*20),:);
                t=t(1:round(stoppoint*20),:);
                iteration=iteration+1;
                break
            end
            n0new=[n(round(stoppoint*20),1:2)'];
            n=n(1:round(stoppoint*20),:);
            t=t(1:round(stoppoint*20),:);
            totaltimeoffB=totaltimeoffB+stoppoint;
        else
            totaltimeoffB=totaltimeoffB+timeoffB;
        end
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,MnewA);
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
            add=stoppoint;
            totaltimeonA=totaltimeonA+stoppoint;
        else
            add=timeonA;
            totaltimeonA=totaltimeonA+timeonA;
        end
        
        iteration=iteration+1;
        
        dn = @(t,n) tumorgrowthmatrixV1_4celltypes( t,n,Moff );
        [t,n] = ode45(dn, tspanoffA, n0new);
        n(:,5)=n(:,1)+n(:,2)+n(:,3)+n(:,4);
        n0new=[n(end,1:4)'];
        if n(end,5)>sum(n0)
            stoppoint=min(t(n(:,5)>sum(n0)));
            if round(stoppoint*20)<=1 || stoppoint<add*(xA-1)
                totaltimeoffA=totaltimeoffA+stoppoint; %THIS IS NEW
                n=n(1:round(stoppoint*20),:);
                t=t(1:round(stoppoint*20),:);
                iteration=iteration+1;
                break
            end
            n0new=[n(round(stoppoint*20),1:2)'];
            n=n(1:round(stoppoint*20),:);
            t=t(1:round(stoppoint*20),:);
            totaltimeoffA=totaltimeoffA+stoppoint;
        else
            totaltimeoffA=totaltimeoffA+timeoffA;
        end
        iteration=iteration+1;
        
    end
    
end
trelapsealternatingbreak=totaltimeonA+totaltimeonB+totaltimeoffA+totaltimeoffB;
advantagealternatingbreak=trelapsealternatingbreak/bestSOC;

%% Combination with breaks
x = max(xA,xB);
maxtimeon=trelapsecombo/x;
timeon=.05 + (maxtimeon)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoff=timeon*(x-1);

Mnew=Moff;
Mnew(1,1)=Moff(1,1)-x*doseA-x*doseB;  %make new drug matrix with higher dose
Mnew(2,2)=Moff(2,2)-x*doseB;
Mnew(3,3)=Moff(3,3)-x*doseA;

n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeoff=0;
iteration=0;
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
            t=t(1:round(stoppoint*20),:);
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
    
    iteration=iteration+1;
    
end
trelapsecombobreak=totaltimeon+totaltimeoff;
advantagecombobreak=trelapsecombobreak/bestSOC; %if greater than 1, break strategy is better

end

