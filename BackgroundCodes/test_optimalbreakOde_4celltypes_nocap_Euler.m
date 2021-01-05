function [advantagecombo,advantagealternating2,advantageoscillating,advantagealternatingbreak,advantagecombobreak] = test_optimalbreakOde_4celltypes_nocap_Euler( MA,MB,Moff,N0)

%  Euler method to add noise at each step

fA=Moff(1,1)-MA(1,1); %drug effect of each drug
fB=Moff(1,1)-MB(1,1);
%Noise parameter
noisevariance = .01*ones(size(N0));
%Simulation parameters
numberTimesteps = 100000;
numberTimesteps = ceil(numberTimesteps);
stepsize = .01;


%% standard of care conditions: each drug indivitually or sequential treatment
%continuous A only
N = zeros(numberTimesteps,length(N0));
T = stepsize*(1:numberTimesteps)';
N(1,:) = N0;
currN = N0;
for jj =2:numberTimesteps
    currdNdT = myderivative4celltypes(currN,MA);
    currN = currN + stepsize*currdNdT;
    currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
    N(jj,:) = currN.*double(currN>0);
end
N(:,5)=sum(N,2);
threshold=sum(N0);
abovethresh=find(N(:,5)>threshold);
minval=min(N(:,5));
singledose_minvalA=minval;
mintime=find(N(:,5)==minval,1,'last');
trelapse_singledoseA=T(abovethresh(find(abovethresh>mintime,1,'first')));
singleTA=T;
singleNA=N;
N0Anew=N(round(trelapse_singledoseA/stepsize),1:4);
if isempty(trelapse_singledoseA)
    advantagecombo=-100;
    advantagealternating2=-100;
    advantageoscillating=-100;
    advantagealternatingbreak=-100;
    advantagecombobreak=-100;
    return
end

%continuous B only
N = zeros(numberTimesteps,length(N0));
T = stepsize*(1:numberTimesteps)';
N(1,:) = N0;
currN = N0;
for jj =2:numberTimesteps
    currdNdT = myderivative4celltypes(currN,MB);
    currN = currN + stepsize*currdNdT;
    currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
    N(jj,:) = currN.*double(currN>0);
end
N(:,5)=sum(N,2);
threshold=sum(N0);
abovethresh=find(N(:,5)>threshold);
minval=min(N(:,5));
singledose_minvalB=minval;
mintime=find(N(:,5)==minval,1,'last');
trelapse_singledoseB=T(abovethresh(find(abovethresh>mintime,1,'first')));
singleTB=T;
singleNB=N;
N0Bnew=N(round(trelapse_singledoseB/stepsize),1:4);
if isempty(trelapse_singledoseB)
    advantagecombo=-100;
    advantagealternating2=-100;
    advantageoscillating=-100;
    advantagealternatingbreak=-100;
    advantagecombobreak=-100;
    return
end

% continuous A then B
N = zeros(numberTimesteps,length(N0));
N(1,:) = N0Anew;
currN = N0Anew;
for jj =2:numberTimesteps
    currdNdT = myderivative4celltypes(currN,MB);
    currN = currN + stepsize*currdNdT;
    currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
    N(jj,:) = currN.*double(currN>0);
end
N(:,5)=sum(N,2);
threshold=sum(N0);
abovethresh=find(N(:,5)>threshold);
minval=min(N(:,5));
singledoseB_minval=minval;
mintime=find(N(:,5)==minval,1,'last');
trelapse_B2=T(abovethresh(find(abovethresh>mintime,1,'first')));
trelapse_AthenB=trelapse_singledoseA+trelapse_B2;  %time to relapse if drug A was maxed out and then drug B was given until relapse
if isempty(trelapse_B2)
    trelapse_AthenB=trelapse_singledoseA;
end
AthenBN=[singleNA(1:round(trelapse_singledoseA/stepsize),:);N];
if isempty(trelapse_AthenB)
    advantagecombo=-100;
    advantagealternating2=-100;
    advantageoscillating=-100;
    advantagealternatingbreak=-100;
    advantagecombobreak=-100;
    return
end

% continuous B then A
N = zeros(numberTimesteps,length(N0));
N(1,:) = N0Bnew;
currN = N0Bnew;
for jj =2:numberTimesteps
    currdNdT = myderivative4celltypes(currN,MA);
    currN = currN + stepsize*currdNdT;
    currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
    N(jj,:) = currN.*double(currN>0);
end
N(:,5)=sum(N,2);
threshold=sum(N0);
abovethresh=find(N(:,5)>threshold);
minval=min(N(:,5));
singledoseA_minval=minval;
mintime=find(N(:,5)==minval,1,'last');
trelapse_A2=T(abovethresh(find(abovethresh>mintime,1,'first')));
trelapse_BthenA=trelapse_singledoseB+trelapse_A2;  %time to relapse if drug B was maxed out and then drug A was given until relapse
if isempty(trelapse_A2)
    trelapse_BthenA=trelapse_singledoseB;
end
BthenAN=[singleNB(1:round(trelapse_singledoseB/stepsize),:);N];

if isempty(trelapse_BthenA)
    advantagecombo=-100;
    advantagealternating2=-100;
    advantageoscillating=-100;
    advantagealternatingbreak=-100;
    advantagecombobreak=-100;
    return
end

bestSOC=max(trelapse_AthenB,trelapse_BthenA);
if trelapse_AthenB>trelapse_BthenA
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

N = zeros(numberTimesteps,length(N0));
N(1,:) = N0;
currN = N0;
for jj =2:numberTimesteps
    currdNdT = myderivative4celltypes(currN,Mon);
    currN = currN + stepsize*currdNdT;
    currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
    N(jj,:) = currN.*double(currN>0);
end
N(:,5)=sum(N,2);
threshold=sum(N0);
abovethresh=find(N(:,5)>threshold);
minval=min(N(:,5));
singledoseB_minval=minval;
mintime=find(N(:,5)==minval,1,'last');
trelapse_combo=T(abovethresh(find(abovethresh>mintime,1,'first')));  %time to relapse if combo was given continuously
comboN=N;

if minval==0
    advantagecombo=200;
else
    advantagecombo=trelapse_combo/bestSOC;
end

%% Alternating 2: no breaks

Atime=.05 + (trelapse_singledoseA)*rand(1); %random amount of time to give drug A, as long as it's less than the time it was given continuously
Btime=.05 + (trelapse_singledoseB)*rand(1);

N0new=N0; %initialize n0new
totaltimeonA=0;
totaltimeonB=0;
dataplot=[];
timeplot=[];
iteration=0;
endtime=0;
numberTimestepsonA = ceil(Atime/stepsize);
numberTimestepsonB = ceil(Btime/stepsize);

if Afirst
    while 1
        if sum(N0new)>sum(N0)  %if tumor exceeds starting size
            break
        end
        
        %on drugA
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonA
            currdNdT = myderivative4celltypes(currN,MA);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonA=totaltimeonA+stoppoint;
        else
            totaltimeonA=totaltimeonA+Atime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonB
            currdNdT = myderivative4celltypes(currN,MB);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonB=totaltimeonB+stoppoint;
        else
            totaltimeonB=totaltimeonB+Btime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
    end
else
    while 1
        if sum(N0new)>sum(N0)  %if tumor exceeds starting size
            break
        end
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonB
            currdNdT = myderivative4celltypes(currN,MB);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonB=totaltimeonB+stoppoint;
        else
            totaltimeonB=totaltimeonB+Btime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %on drugA
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonA
            currdNdT = myderivative4celltypes(currN,MA);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonA=totaltimeonA+stoppoint;
        else
            totaltimeonA=totaltimeonA+Atime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
    end
end
trelapsealternating2=totaltimeonA+totaltimeonB;
if minval==0
    advantagealternating2=200;
else
    advantagealternating2=trelapsealternating2/bestSOC;
end

%% Oscillation: not used in data analysis

N0new=N0; %initialize n0new
totaltimeonA=0;
totaltimeonB=0;
totaltimeoff=0;
trelapse=0;
dataplot=[];
timeplot=[];
iteration=0;
endtime=0;
numberTimestepsonA = ceil(Atime/stepsize);
numberTimestepsonB = ceil(Btime/stepsize);

if Afirst
    while 1
        if sum(N0new)>sum(N0)  %if tumor exceeds starting size
            break
        end
        
        %on drugA
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonA
            currdNdT = myderivative4celltypes(currN,MA);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonA=totaltimeonA+stoppoint;
        else
            totaltimeonA=totaltimeonA+Atime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        N = zeros(numberTimesteps,length(N0));
        T = stepsize*(1:numberTimesteps)';
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimesteps
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        threshold=sum(N0);
        abovethresh=find(N(:,5)>threshold);
        minval=min(N(:,5));
        if minval==0
            break
        end
        mintime=find(N(:,5)==minval,1,'last');
        trelapse=T(abovethresh(find(abovethresh>mintime,1,'first')));
        N0new=N(round(trelapse/stepsize),1:4);
        
        totaltimeoff=totaltimeoff+trelapse;
        iteration=iteration+1;
        dataplot=[dataplot;N(1:round(trelapse/stepsize),:)];
        newT=T(1:round(trelapse/stepsize))+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonB
            currdNdT = myderivative4celltypes(currN,MB);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonB=totaltimeonB+stoppoint;
        else
            totaltimeonB=totaltimeonB+Btime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        N = zeros(numberTimesteps,length(N0));
        T = stepsize*(1:numberTimesteps)';
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimesteps
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        threshold=sum(N0);
        abovethresh=find(N(:,5)>threshold);
        minval=min(N(:,5));
        if minval==0
            break
        end
        mintime=find(N(:,5)==minval,1,'last');
        trelapse=T(abovethresh(find(abovethresh>mintime,1,'first')));
        N0new=N(round(trelapse/stepsize),1:4);
        
        totaltimeoff=totaltimeoff+trelapse;
        iteration=iteration+1;
        dataplot=[dataplot;N(1:round(trelapse/stepsize),:)];
        newT=T(1:round(trelapse/stepsize))+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        
    end
else
    while 1
        if sum(N0new)>sum(N0)  %if tumor exceeds starting size
            break
        end
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonB
            currdNdT = myderivative4celltypes(currN,MB);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonB=totaltimeonB+stoppoint;
        else
            totaltimeonB=totaltimeonB+Btime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        N = zeros(numberTimesteps,length(N0));
        T = stepsize*(1:numberTimesteps)';
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimesteps
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        threshold=sum(N0);
        abovethresh=find(N(:,5)>threshold);
        minval=min(N(:,5));
        if minval==0
            break
        end
        mintime=find(N(:,5)==minval,1,'last');
        trelapse=T(abovethresh(find(abovethresh>mintime,1,'first')));
        N0new=N(round(trelapse/stepsize),1:4);
        
        totaltimeoff=totaltimeoff+trelapse;
        iteration=iteration+1;
        dataplot=[dataplot;N(1:round(trelapse/stepsize),:)];
        newT=T(1:round(trelapse/stepsize))+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        
        %on drugA
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonA
            currdNdT = myderivative4celltypes(currN,MA);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeonA=totaltimeonA+stoppoint;
        else
            totaltimeonA=totaltimeonA+Atime;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        N = zeros(numberTimesteps,length(N0));
        T = stepsize*(1:numberTimesteps)';
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimesteps
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        threshold=sum(N0);
        abovethresh=find(N(:,5)>threshold);
        minval=min(N(:,5));
        if minval==0
            break
        end
        mintime=find(N(:,5)==minval,1,'last');
        trelapse=T(abovethresh(find(abovethresh>mintime,1,'first')));
        N0new=N(round(trelapse/stepsize),1:4);
        
        totaltimeoff=totaltimeoff+trelapse;
        iteration=iteration+1;
        dataplot=[dataplot;N(1:round(trelapse/stepsize),:)];
        newT=T(1:round(trelapse/stepsize))+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
    end
    
end
totaltimeon=totaltimeonA+totaltimeonB;
trelapseoscillation=totaltimeon+totaltimeoff;
if minval==0
    advantageoscillating=200;
else
    advantageoscillating=trelapseoscillation/bestSOC;
end

%% Alternating with breaks

xA = 1 + (9)*rand(1); %random dose escalation up to 10 times original dose
maxtimeonA=trelapse_singledoseA/xA; %maximum time this dose can be given is 1/x amount of time
timeonA=.05 + (maxtimeonA)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoffA=timeonA*(xA-1);

xB=xA;
maxtimeonB=trelapse_singledoseB/xB; %maximum time this dose can be given is 1/x amount of time
timeonB=.05 + (maxtimeonB)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoffB=timeonB*(xB-1);

MnewA=MA;
MnewA(1,1)=Moff(1,1)-xA*fA;  %make new drug matrix with higher dose, affecting only those lines that are sensitive to A
MnewA(3,3)=Moff(3,3)-xA*fA;
MnewB=MB;
MnewB(1,1)=Moff(1,1)-xB*fB;  %make new drug matrix with higher dose, affecting only those lines that are sensitive to B
MnewB(2,2)=Moff(2,2)-xB*fB;

N0new=N0; %initialize n0new
totaltimeonA=0;
totaltimeonB=0;
totaltimeoffA=0;
totaltimeoffB=0;
iteration=0;
dataplot=[];
timeplot=[];
endtime=0;
numberTimestepsonA = ceil(timeonA/stepsize);
numberTimestepsoffA = ceil(timeoffA/stepsize);
numberTimestepsonB = ceil(timeonB/stepsize);
numberTimestepsoffB = ceil(timeoffB/stepsize);


if Afirst
    while 1
        if sum(N0new)>sum(N0)  %if tumor exceeds starting size
            break
        end
        
        %on drugA
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonA
            currdNdT = myderivative4celltypes(currN,MnewA);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            add=stoppoint;
            totaltimeonA=totaltimeonA+stoppoint;
        else
            add=timeonA;
            totaltimeonA=totaltimeonA+timeonA;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        T = stepsize*(1:numberTimestepsoffA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsoffA
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1 || stoppoint<add*(xA-1)
                totaltimeoffA=totaltimeoffA+stoppoint;
                N=N(1:round(stoppoint*100),:);
                T=T(1:round(stoppoint*100),:);
                dataplot=[dataplot;N];
                newT=T+endtime;
                timeplot=[timeplot;newT];
                iteration=iteration+1;
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeoffA=totaltimeoffA+stoppoint;
        else
            totaltimeoffA=totaltimeoffA+timeoffA;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonB
            currdNdT = myderivative4celltypes(currN,MnewB);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            add=stoppoint;
            totaltimeonB=totaltimeonB+stoppoint;
        else
            add=timeonB;
            totaltimeonB=totaltimeonB+timeonB;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        T = stepsize*(1:numberTimestepsoffB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsoffB
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1 || stoppoint<add*(xB-1)
                totaltimeoffB=totaltimeoffB+stoppoint;
                N=N(1:round(stoppoint*100),:);
                T=T(1:round(stoppoint*100),:);
                dataplot=[dataplot;N];
                newT=T+endtime;
                timeplot=[timeplot;newT];
                iteration=iteration+1;
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeoffB=totaltimeoffB+stoppoint;
        else
            totaltimeoffB=totaltimeoffB+timeoffB;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
    end
else
    while 1
        if sum(N0new)>sum(N0)  %if tumor exceeds starting size
            break
        end
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonB
            currdNdT = myderivative4celltypes(currN,MnewB);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            add=stoppoint;
            totaltimeonB=totaltimeonB+stoppoint;
        else
            add=timeonB;
            totaltimeonB=totaltimeonB+timeonB;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        T = stepsize*(1:numberTimestepsoffB)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsoffB
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1 || stoppoint<add*(xB-1)
                totaltimeoffB=totaltimeoffB+stoppoint;
                N=N(1:round(stoppoint*100),:);
                T=T(1:round(stoppoint*100),:);
                dataplot=[dataplot;N];
                newT=T+endtime;
                timeplot=[timeplot;newT];
                iteration=iteration+1;
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeoffB=totaltimeoffB+stoppoint;
        else
            totaltimeoffB=totaltimeoffB+timeoffB;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %on drugB
        clear T
        clear N
        T = stepsize*(1:numberTimestepsonA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsonA
            currdNdT = myderivative4celltypes(currN,MnewA);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            add=stoppoint;
            totaltimeonA=totaltimeonA+stoppoint;
        else
            add=timeonA;
            totaltimeonA=totaltimeonA+timeonA;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
        
        %off drug
        clear T
        clear N
        T = stepsize*(1:numberTimestepsoffA)';
        N=[];
        N(1,:) = N0new;
        currN = N0new;
        for jj =2:numberTimestepsoffA
            currdNdT = myderivative4celltypes(currN,Moff);
            currN = currN + stepsize*currdNdT;
            currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
            currN=currN(1,:);
            N(jj,:) = currN.*double(currN>0);
        end
        N(:,5)=sum(N,2);
        N0new=N(end,1:4);
        if N(end,5)>sum(N0)
            stoppoint=min(T(N(:,5)>sum(N0)));
            if round(stoppoint*100)<=1 || stoppoint<add*(xA-1)
                totaltimeoffA=totaltimeoffA+stoppoint;
                N=N(1:round(stoppoint*100),:);
                T=T(1:round(stoppoint*100),:);
                dataplot=[dataplot;N];
                newT=T+endtime;
                timeplot=[timeplot;newT];
                iteration=iteration+1;
                break
            end
            N0new=[N(round(stoppoint*100),1:4)];
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            totaltimeoffA=totaltimeoffA+stoppoint;
        else
            totaltimeoffA=totaltimeoffA+timeoffA;
        end
        
        iteration=iteration+1;
        dataplot=[dataplot;N];
        newT=T+endtime;
        timeplot=[timeplot;newT];
        endtime=timeplot(end);
        minval=min(N(:,5));
        if minval==0
            break
        end
    end
end
trelapsealternatingbreak=totaltimeonA+totaltimeonB+totaltimeoffA+totaltimeoffB;
if minval==0
    advantagealternatingbreak=200;
else
    advantagealternatingbreak=trelapsealternatingbreak/bestSOC;
end
p=1;

%% Combination with breaks
x = max(xA,xB);
maxtimeon=trelapse_combo/x;
timeon=.05 + (maxtimeon)*rand(1); %random amount of time to give this higher dose, as long as this time is less than maximum time dose can be given
timeoff=timeon*(x-1);

Mnew=Moff;
Mnew(1,1)=Moff(1,1)-x*doseA-x*doseB;  %make new drug matrix with higher dose
Mnew(2,2)=Moff(2,2)-x*doseB;
Mnew(3,3)=Moff(3,3)-x*doseA;

N0new=N0; %initialize n0new
totaltimeon=0;
totaltimeoff=0;
iteration=0;
dataplot=[];
timeplot=[];
endtime=0;
numberTimestepson = ceil(timeon/stepsize);
numberTimestepsoff = ceil(timeoff/stepsize);
if length(numberTimestepson)<1 || length(numberTimestepsoff)<1
    advantagecombobreak=0;
    return
end

while 1
    if sum(N0new)>sum(N0)  %if tumor exceeds starting size
        break
    end
    
    %on drugA
    clear T
    clear N
    T = stepsize*(1:numberTimestepson)';
    N=[];
    N(1,:) = N0new;
    currN = N0new;
    for jj =2:numberTimestepson
        currdNdT = myderivative4celltypes(currN,Mnew);
        currN = currN + stepsize*currdNdT;
        currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
        currN=currN(1,:);
        N(jj,:) = currN.*double(currN>0);
    end
    N(:,5)=sum(N,2);
    N0new=N(end,1:4);
    if N(end,5)>sum(N0)
        stoppoint=min(T(N(:,5)>sum(N0)));
        if round(stoppoint*100)<=1
            break
        end
        N0new=[N(round(stoppoint*100),1:4)];
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
    minval=min(N(:,5));
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
        currdNdT = myderivative4celltypes(currN,Moff);
        currN = currN + stepsize*currdNdT;
        currN = currN + sqrt(abs(currN)).*noisevariance.*randn(size(N0));
        currN=currN(1,:);
        N(jj,:) = currN.*double(currN>0);
    end
    N(:,5)=sum(N,2);
    N0new=N(end,1:4);
    if N(end,5)>sum(N0)
        stoppoint=min(T(N(:,5)>sum(N0)));
        if round(stoppoint*100)<=1 || stoppoint<add*(x-1)
            totaltimeoff=totaltimeoff+stoppoint;
            N=N(1:round(stoppoint*100),:);
            T=T(1:round(stoppoint*100),:);
            dataplot=[dataplot;N];
            newT=T+endtime;
            timeplot=[timeplot;newT];
            iteration=iteration+1;
            break
        end
        N0new=[N(round(stoppoint*100),1:4)];
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
    minval=min(N(:,5));
    if minval==0
        break
    end
end
trelapsecombobreak=totaltimeon+totaltimeoff;
if minval==0
    advantagecombobreak=200;
else
    advantagecombobreak=trelapsecombobreak/bestSOC;
end

end