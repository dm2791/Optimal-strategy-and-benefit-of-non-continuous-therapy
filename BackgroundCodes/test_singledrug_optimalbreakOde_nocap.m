function [ breakadvantage,x,timeon,timeoff,totaltimeon,totaltimeoff] = test_singledrug_optimalbreakOde_nocap( M,Moff,n0)

f=Moff(1,1)-M(1,1); %drug effect

tspan=[0:.05:500];

dn = @(t,n) tumorgrowthmatrixV1( t,n,M );
[t,n] = ode45(dn, tspan, n0);
n(:,3)=n(:,1)+n(:,2);
trelapse=min(t(n(:,3)>sum(n0))); %time to relapse if drug was given continuously
singledosen = n(1:round(trelapse*20),3);
singledose_trelapse=trelapse;


x = 1 + (9)*rand(1); %random dose escalation up to 10 times original dose
maxtimeon=trelapse/x; %maximum time this dose can be given is 1/x amount of time
timeon=.05 + (maxtimeon)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoff=timeon*(x-1);

Mnew=M;
Mnew(1,1)=Moff(1,1)-x*f;  %make new drug matrix with higher dose

n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeoff=0;
iteration=0;
tspanon=[0:.05:timeon];
tspanoff=[0:.05:timeoff];
if length(tspanon)<2 || length(tspanoff)<2
    breakadvantage=0;
    totaltimeon=0;
    totaltimeoff=0;
    return
end

while 1
    if sum(n0new)>sum(n0)  %if tumor exceeds starting size
        break
    end
    
    dn = @(t,n) tumorgrowthmatrixV1( t,n,Mnew );
    [t,n] = ode45(dn, tspanon, n0new);
    n(:,3)=n(:,1)+n(:,2);
    n0new=[n(end,1:2)'];  %calculate n at end of on interval; will be new n0 for start of off interval
    if n(end,3)>sum(n0)
        stoppoint=min(t(n(:,3)>sum(n0)));
        if round(stoppoint*20)<=1
            break
        end
        n0new=[n(round(stoppoint*20),1:2)'];
        n=n(1:round(stoppoint*20),:);
        t=t(1:round(stoppoint*20),:);
        add=stoppoint;
        totaltimeon=totaltimeon+stoppoint;
    else
        add=timeon;
        totaltimeon=totaltimeon+timeon;
    end
    
    iteration=iteration+1;
    
    dn = @(t,n) tumorgrowthmatrixV1( t,n,Moff );
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,3)=n(:,1)+n(:,2);
    n0new=[n(end,1:2)'];
    if n(end,3)>sum(n0)
        stoppoint=min(t(n(:,3)>sum(n0)));
        if round(stoppoint*20)<=1 || stoppoint<add*(x-1)
            totaltimeoff=totaltimeoff+stoppoint;
            n=n(1:round(stoppoint*20),:);
            t=t(1:round(stoppoint*20),:);
            iteration=iteration+1;
            break
        end
        n0new=[n(round(stoppoint*20),1:2)'];
        n=n(1:round(stoppoint*20),:);
        t=t(1:round(stoppoint*20),:);
        totaltimeoff=totaltimeoff+stoppoint;
    else
        totaltimeoff=totaltimeoff+timeoff;
    end

    iteration=iteration+1;
    
end

time=totaltimeon+totaltimeoff;
breakadvantage=time/trelapse;




end

