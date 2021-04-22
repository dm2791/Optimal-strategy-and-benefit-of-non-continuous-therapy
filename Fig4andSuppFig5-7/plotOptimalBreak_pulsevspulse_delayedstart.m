function [breakadvantage,breakadvantage_base,singledose_trelapse,singledosen,time,dataplot_in,dataplot,timeplot,trelapse] = plotOptimalBreak_pulsevspulse_delayedstart(M,Moff,n0,x_in,x,timeon_in,timeoff_in,timeon,token,startpulse)

f=Moff(1,1)-M(1,1); %drug effect
f_base=f/x_in; %baseline dose if given continuously

tspan=[0:.05:500];
M_base=Moff;
M_base(1,1)=Moff(1,1)-f_base;
dn = @(t,n) tumorgrowthmatrixV1( t,n,M_base );
[t,n] = ode45(dn, tspan, n0);
n(:,3)=n(:,1)+n(:,2);
[c ind]=min(n(:,3));
if token==0
    trelapse=min(t(n(:,3)>sum(n0)))*10; %time to relapse if drug was given continuously
elseif token==1
    trelapse=t(min(find(n(ind:end,3)>256))+ind-1)*10;
elseif token==2
    trelapse=t(min(find(n(ind:end,3)>252))+ind-1)*10;
end
singledosen = n(1:round(trelapse*2),3);
singledose_trelapse=trelapse;


tspanon=[0:.05:timeon_in/10];
tspanoff=[0:.05:timeoff_in/10];
n0new=n0;
dataplot_in=[];
timeplot_in=[];
iteration=0;
while sum(n0new)<sum(n0)+1
    dn = @(t,n) tumorgrowthmatrixV1( t,n,M );
    [t,n] = ode45(dn, tspanon, n0new);
    n(:,3)=n(:,1)+n(:,2);
    dataplot_in=[dataplot_in;n];
    if iteration==0
        timeplot_in=[timeplot_in;t];
    else
        timeplot_in=[timeplot_in;t+.05+timeplot_in(end)];
    end
    n0new=[n(end,1:2)'];
    iteration=iteration+1;
    dn = @(t,n) tumorgrowthmatrixV1( t,n,Moff );
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,3)=n(:,1)+n(:,2);
    dataplot_in=[dataplot_in;n];
    timeplot_in=[timeplot_in;t+.05+timeplot_in(end)];
    
    n0new=[n(end,1:2)'];
    iteration=iteration+1;
    if iteration==startpulse
        dataplot_save=dataplot_in(1:end-1,:);
        timeplot_save=timeplot_in(1:end-1);
        n0_save=n0new;
    end
end
[c ind]=min(dataplot_in(:,3));
if token==0
    trelapse=min(find(dataplot_in(:,3)>sum(n0)));
elseif token==1
    trelapse=min(find(dataplot_in(ind:end,3)>256))+ind-1; %%
elseif token==2
    trelapse=min(find(dataplot_in(ind:end,3)>252))+ind-1; %%
end
trelapse=timeplot_in(trelapse)*10;
breakadvantage_base=trelapse/singledose_trelapse;
dataplot_in = dataplot_in(1:round(trelapse*2),:);

timeoff=timeon*(x*x_in-1);
Mnew=M;
Mnew(1,1)=Moff(1,1)-x*f;  %make new drug matrix with higher dose

%n0new=n0; %initialize n0new
n0new=n0_save;
%dataplot=[];
dataplot=dataplot_save;
timeplot=timeplot_save;
tspanon=[0:.05:timeon/10];
tspanoff=[0:.05:timeoff/10];
if length(tspanon)<2 || length(tspanoff)<2
    breakadvantage=0;
    totaltimeon=0;
    totaltimeoff=0;
    return
end

iteration=startpulse*2;
while 1
    if sum(n0new)>sum(n0)  %if tumor exceeds starting size
        break
    end
    
    dn = @(t,n) tumorgrowthmatrixV1( t,n,Mnew );
    [t,n] = ode45(dn, tspanon, n0new);
    n(:,3)=n(:,1)+n(:,2);
    n0new=[n(end,1:2)'];  %calculate n at end of on interval; will be new n0 for start of off interval
    
    dataplot=[dataplot;n];
    if iteration==0
        timeplot=[timeplot;t];
    else
        timeplot=[timeplot;t+.05+timeplot(end)];
    end
    iteration=iteration+1;
    dn = @(t,n) tumorgrowthmatrixV1( t,n,Moff );
    [t,n] = ode45(dn, tspanoff, n0new);
    n(:,3)=n(:,1)+n(:,2);
    n0new=[n(end,1:2)'];
    
    dataplot=[dataplot;n];
    timeplot=[timeplot;t+.05+timeplot(end)];
    iteration=iteration+1;
    
end
[c ind]=min(dataplot(:,3));
if token==0
    time=min(find(dataplot(:,3)>sum(n0)));
elseif token==1
    time=min(find(dataplot(ind:end,3)>256))+ind-1;
elseif token==2
    time=min(find(dataplot(ind:end,3)>252))+ind-1;
end
time=timeplot(time)*10;
breakadvantage=time/trelapse;
dataplot = dataplot(1:round(time*2),:);



end



