function [breakadvantage,time,dataplot,timeplot,trelapse,singledosen] = plotOptimalBreak_delayedstart(M,Moff,n0,x,timeon,startpulse,halflife)

%incorporate half life

f=Moff(1,1)-M(1,1); %drug effect
timeon=(timeon-1)/20;

lambda=log(2)/halflife;

tspan=[0:.05:500];

dn = @(t,n) tumorgrowthmatrixV1( t,n,M );
[t,n] = ode45(dn, tspan, n0);
n(:,3)=n(:,1)+n(:,2);
trelapse=min(t(n(:,3)>sum(n0))); %time to relapse if drug was given continuously
singledosen = n(1:round(trelapse*20),3);
singledose_trelapse=trelapse;

timeoff=timeon*(x-1);
% timeoff=round(timeoff,2);
% timeon=round(timeon,2);

Mnew=M;
Mnew(1,1)=Moff(1,1)-x*f;  %make new drug matrix with higher dose
Moffnew=Moff;
fnew=Mnew(1,1)-Moff(1,1); %drug effect

% find start point
starttime=find(t*10<startpulse,1,'last');
n0new=n(starttime+1,1:2);
%n0new=n0; %initialize n0new
totaltimeon=0;
totaltimeoff=0;
iteration=0;
%dataplot=[];
dataplot=n(1:starttime,:);
%timeplot=[];
timeplot=t(1:starttime);
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
  
    dataplot=[dataplot;n];
    timeplot=[timeplot;t+.05+timeplot(end)];
    iteration=iteration+1;
    
    %dn = @(t,n) tumorgrowthmatrixV1( t,n,Moff );
    %dn = @(t,n) tumorgrowthmatrixV1_nonconstant( t,n,(Moff(1,1)+(fnew)*exp(-lambda*t*20)),Moff(2,2),Moff(2,1) );
    dn = @(t,n) tumorgrowthmatrixV1_nonconstant( t,n,Mnew(1,1),Moff(2,2),Moff(2,1) ,Moff(1,1), lambda );
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,3)=n(:,1)+n(:,2);
n0new=[n(end,1:2)'];
        if n(end,3)>sum(n0)
                stoppoint=min(t(n(:,3)>sum(n0)));
        if round(stoppoint*20)<=1 || stoppoint<add*(x-1)
             totaltimeoff=totaltimeoff+stoppoint; 
             n=n(1:round(stoppoint*20),:);
              t=t(1:round(stoppoint*20),:);
             dataplot=[dataplot;n];
             timeplot=[timeplot;t+.05+timeplot(end)];
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
    
    dataplot=[dataplot;n];
    timeplot=[timeplot;t+.05+timeplot(end)];
      iteration=iteration+1;

end

time=totaltimeon+totaltimeoff+startpulse/10;
breakadvantage=time/trelapse;




end



