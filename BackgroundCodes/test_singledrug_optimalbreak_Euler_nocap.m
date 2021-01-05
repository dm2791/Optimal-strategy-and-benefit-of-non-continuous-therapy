function [ breakadvantage,x,timeon,totaltimeon,totaltimeoff,iteration] = test_singledrug_optimalbreak_Euler_nocap( M,Moff,N0)
%test defined on and off periods
%  Euler method to add noise at each step

f=Moff(1,1)-M(1,1); %drug effect
%Noise parameter
noisevariance = .1*ones(size(N0));

%Simulation parameters
numberTimesteps = 100000;
numberTimesteps = ceil(numberTimesteps);
stepsize = .01;

%% simulate single dose
N = zeros(numberTimesteps,length(N0));
T = stepsize*(1:numberTimesteps)';
N(1,:) = N0;
currN = N0;
for jj =2:numberTimesteps
    currdNdT = myderivative(currN,M);
    currN = currN + stepsize*currdNdT;
    currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
    N(jj,:) = currN.*double(currN>0);
end
N(:,3)=sum(N,2);
threshold=sum(N0);
abovethresh=find(N(:,3)>threshold);
minval=min(N(:,3));
singledose_minval=minval;
mintime=find(N(:,3)==minval,1,'last');
singledose_trelapse=T(abovethresh(find(abovethresh>mintime,1,'first')));
if isempty(singledose_trelapse) || minval==0
    breakadvantage=-100;
    x=0;
    timeon=0;
    totaltimeon=0;
    totaltimeoff=0;
    iteration=0;
    return
end
singleT=T;
singleN=N;

%% setup
x = 1 + (9)*rand(1); %random dose escalation up to 10 times original dose
maxtimeon=singledose_trelapse/x; %maximum time this dose can be given is 1/x amount of time
timeon=.05 + (maxtimeon)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoff=timeon*(x-1);

Mnew=M;
Mnew(1,1)=Moff(1,1)-x*f;  %make new drug matrix with higher dose

N0new=N0; %initialize n0new
totaltimeon=0;
totaltimeoff=0;
dataplot=[];
timeplot=[];
numberTimestepson = ceil(timeon/stepsize);
numberTimestepsoff = ceil(timeoff/stepsize);
endtime=0;
iteration=0;

%%
while 1
    if sum(N0new)>sum(N0)  %if tumor exceeds starting size
        break
    end
    
    %on drug
    clear T
    clear N
    T = stepsize*(1:numberTimestepson)';
    N=[];
    N(1,:) = N0new;
    currN = N0new;
    for jj =2:numberTimestepson
        currdNdT = myderivative(currN,Mnew);
        currN = currN + stepsize*currdNdT;
        currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
        N(jj,:) = currN.*double(currN>0);
    end
    N(:,3)=sum(N,2);
    
    N0new=(N(end,1:2));
    if N(end,3)>sum(N0)
        stoppoint=min(T(N(:,3)>sum(N0)));
        if round(stoppoint*100)<=1
            break
        end
        N0new=[N(round(stoppoint*100),1:2)];
        N=N(1:round(stoppoint*100),:);
        T=T(1:round(stoppoint*100),:);
        add=stoppoint;
        totaltimeon=totaltimeon+stoppoint;
    else
        add=timeon;
        totaltimeon=totaltimeon+timeon;
    end
    
    iteration=iteration+1;
    
    dataplot=[dataplot;N];
    newT=T+endtime;
    timeplot=[timeplot;newT];
    endtime=timeplot(end);
    minval=min(N(:,3));
    if minval==0
        break
    end
    
    %off drug
    clear T
    clear N
    T = stepsize*(1:numberTimestepsoff)';
    N=[];
    N(1,:) = N0new;
    currN = N0new;
    for jj =2:numberTimestepsoff
        currdNdT = myderivative(currN,Moff);
        currN = currN + stepsize*currdNdT;
        currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
        N(jj,:) = currN.*double(currN>0);
    end
    N(:,3)=sum(N,2);
    
    N0new=round(N(end,1:2));
    if N(end,3)>sum(N0)
        stoppoint=min(T(N(:,3)>sum(N0)));
        if round(stoppoint*100)<=1 || stoppoint<add*(x-1)
            totaltimeoff=totaltimeoff+stoppoint;
            N=N(1:round(stoppoint*100),:);
            dataplot=[dataplot;N];
            iteration=iteration+1;
            break
        end
        N0new=[N(round(stoppoint*100),1:2)];
        N=N(1:round(stoppoint*100),:);
        T=T(1:round(stoppoint*100),:);
        totaltimeoff=totaltimeoff+stoppoint;
    else
        totaltimeoff=totaltimeoff+timeoff;
    end
    
    iteration=iteration+1;
    dataplot=[dataplot;N];
    newT=T+endtime;
    timeplot=[timeplot;newT];
    endtime=timeplot(end);
    minval=min(N(:,3));
    if minval==0
        break
    end
    
end
time=totaltimeon+totaltimeoff;

if singledose_trelapse==0
    breakadvantage=100;
elseif minval==0
    breakadvantage=200;
else
    breakadvantage=time/singledose_trelapse;
end

end