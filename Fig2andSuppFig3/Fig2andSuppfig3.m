%% load data
alldata1=readmatrix('2drugsimulationdata_paper.csv');
alldata2=readmatrix('2drugsimulation_paper_revpart2.csv');
alldata=[alldata1;alldata2];

%% Unsupervised analysis
advantagescut=alldata(:,53:57);

matrices=alldata(:,1:52);
[coeff,score,latent,tsquared,explained,mu] = pca(matrices);
conditions=ones(size(advantagescut,1),1);
conditionsadd=advantagescut(:,1:5);
conditions=[conditions,conditionsadd];

conditions(:,4)=[];
color=[];
combo=0;
combobreaks=0;
alternatingbreaks=0;
alternating2=0;
soc=0;
for i=1:size(conditions,1)
    if find(conditions(i,:)==max(conditions(i,:)))==2;
        color(i,:)=[0.5843,0.8157,0.9882]; %light blue if combo best
        combo=combo+1;
    elseif find(conditions(i,:)==max(conditions(i,:)))==5;
        color(i,:)=[0,0,1]; %blue if combobreaksbest
        combobreaks=combobreaks+1;
    elseif find(conditions(i,:)==max(conditions(i,:)))==4;
        color(i,:)=[1,0,1]; %magenta if alternating breaks best
        alternatingbreaks=alternatingbreaks+1;
    elseif find(conditions(i,:)==max(conditions(i,:)))==3
        color(i,:)=[1,0,0]; %red if alternating2 best
        alternating2=alternating2+1;
    elseif find(conditions(i,:)==max(conditions(i,:)))==1;
        color(i,:)=[0,0,0]; %black if sequential best
        soc=soc+1;
        i;
    end
end


figure
scatter3(score(:,1),score(:,2),score(:,3),50,color,'filled')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
figure,pie([soc,combo,alternating2,combobreaks,alternatingbreaks])

%% find correlations between outputs
[coeff,score,latent,tsquared,explained,mu] = pca(conditions);
figure
biplot(coeff(:,1:2),'Scores',score(:,1:2),'Varlabels',{'1','2','3','4','5'}');
xlabel('component1')
ylabel('component2')
%% consider only conditions that don't need dose escalation
conditions_cont=conditions(:,1:3);
color=[];
combo_cont=0;
alternating_cont=0;
soc_cont=0;
for i=1:size(conditions,1)
    if find(conditions_cont(i,:)==max(conditions_cont(i,:)))==2;
        color(i,:)=[0.5843,0.8157,0.9882]; %light blue if combo best
        combo_cont=combo_cont+1;
    elseif find(conditions_cont(i,:)==max(conditions_cont(i,:)))==3
        color(i,:)=[1,0,0]; %red if alternating2 best
        alternating_cont=alternating_cont+1;
    elseif find(conditions_cont(i,:)==max(conditions_cont(i,:)))==1;
        color(i,:)=[0,0,0]; %black if sequential best
        soc_cont=soc_cont+1;
    end
end

best_cont=[alternating_cont,combo_cont,soc_cont]./length(conditions);
figure,pie(best_cont)

%%
%compare percentage of success of each winning strategy
successpercent=alldata(:,58:62);
successpercent(:,3)=[];
meansuccess=mean(successpercent);
figure,bar(meansuccess)
ylabel('%success')
title('mean success rate for a given tumor')

meanad=mean(conditions);
figure,bar(meanad(2:end))
ylabel('log(pulse advantage)')
title('Mean advantage')

maxad=max(conditions);
figure,bar(maxad(2:end))
ylabel('log(pulse advantage)')
title('Maximum advantage')
%% Make Matrix
logc=log(conditions);

matrices(:,53)=alldata(:,63); %sum(Moff*n0)
matrices(:,54)=alldata(:,64); %sum(MA*n0)
matrices(:,55)=alldata(:,65); %sum(MB*n0)
matrices(:,56)=matrices(:,49)/100; %prop double sensitive
matrices(:,57)=(matrices(:,50)+matrices(:,52))/100; %prop resistant to A
matrices(:,58)=(matrices(:,51)+matrices(:,52))/100; %prop resistant to B
matrices(:,59)=matrices(:,52)/100; %prop double resistant
matrices(:,60)=matrices(:,1)./(matrices(:,6)+matrices(:,11)+matrices(:,16)); %Moff; growthrate double sensitive/others
matrices(:,61)=(matrices(:,1)+matrices(:,11))./(matrices(:,6)+matrices(:,16)); %Moff; growthrate sensitiveA/resistantA
matrices(:,62)=matrices(:,1)./(matrices(:,11)+matrices(:,6)+matrices(:,16)); %Moff; growthrate sensitiveB/resistantB
matrices(:,63)=matrices(:,16)./(matrices(:,1)+matrices(:,6)+matrices(:,11)); %Moff; growthrate double resistant/others
matrices(:,64)=matrices(:,17)./(matrices(:,22)+matrices(:,27)+matrices(:,32)); %MA; growthrate double sensitive/others
matrices(:,65)=(matrices(:,17)+matrices(:,27))./(matrices(:,22)+matrices(:,32)); %MA; growthrate sensitiveA/resistantA
matrices(:,66)=(matrices(:,17)+matrices(:,22))./(matrices(:,27)+matrices(:,32)); %MA; growthrate sensitiveB/resistantB
matrices(:,67)=matrices(:,32)./(matrices(:,17)+matrices(:,22)+matrices(:,27)); %MA; growthrate double resistant/others
matrices(:,68)=matrices(:,33)./(matrices(:,38)+matrices(:,43)+matrices(:,48)); %MB; growthrate double sensitive/others
matrices(:,69)=(matrices(:,33)+matrices(:,43))./(matrices(:,38)+matrices(:,48)); %MB; growthrate sensitiveA/resistantA
matrices(:,70)=(matrices(:,33)+matrices(:,38))./(matrices(:,43)+matrices(:,48)); %MB; growthrate sensitiveB/resistantB
matrices(:,71)=matrices(:,48)./(matrices(:,43)+matrices(:,38)+matrices(:,33)); %MB; growthrate double resistant/others
matrices(:,72)=matrices(:,1)-matrices(:,17); %fA
matrices(:,73)=matrices(:,1)-matrices(:,33); %fB
matrices(:,74)=(matrices(:,18)./matrices(:,2))+(matrices(:,28)./matrices(:,12)); %uA/uoff
matrices(:,75)=(matrices(:,35)./matrices(:,3))+(matrices(:,40)./matrices(:,8)); %uB/uoff
%%
split=round(length(alldata)*.6);
testsetmatrices=matrices(split+1:end,:);
matrices=matrices(1:split,:);
testsetlogc=logc(split+1:end,:);
logc=logc(1:split,:);

%% predict strength of accuracy for top  strategies (alternating breaks)
table=array2table(matrices);
tableselect=array2table(matrices(:,53:75));
t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); 
Mdl = fitrensemble(tableselect,logc(:,4),'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2

impOOB = oobPermutedPredictorImportance(Mdl);
figure
bar(sort(impOOB))
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = table.Properties.VariableNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none'; 

%% combo breaks
t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); 
Mdl = fitrensemble(tableselect,logc(:,5),'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2 

impOOB = oobPermutedPredictorImportance(Mdl);
figure
bar(sort(impOOB))
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = table.Properties.VariableNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';   

%% Alternating
t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); 
%remove rows where alternating advantage =0 because mdl can't handle that
indzero=find(logc(:,3)==0);
logc2test=logc(:,3);
logc2test(indzero,:)=[];
tableselect2test=tableselect;
tableselect2test(indzero,:)=[];
indinf=find(logc2test==-Inf);
logc2test(indinf,:)=[];
tableselect2test(indinf,:)=[];


Mdl = fitrensemble(tableselect2test,logc2test,'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2 

impOOB = oobPermutedPredictorImportance(Mdl);
figure
bar(sort(impOOB))
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = table.Properties.VariableNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';  

%% binarize method: predict winner
winner=[];
for i=1:size(conditions,1)
    winner(i)=find(conditions(i,:)==max(conditions(i,:)));
end
winner=winner';
% combine alternating breaks and combobreaks; binarize method: predict winner
altb=find(winner==4); %find altbreaks
winner2=winner;
winner2(altb)=5; %put as combo breaks

[ranks,weights] = relieff(table2array(tableselect),winner2(1:split),5); 
figure,bar(weights(ranks))
ylabel('predictor importance weight')
xlabel('predictor rank')

MdlReduced=fitcecoc(table2array(tableselect(:,{'Var9' 'Var18' 'Var8'})),winner2(1:split));
isLoss = resubLoss(MdlReduced) 
tableselecttest=array2table(testsetmatrices);
selecttest=tableselecttest(:,{'testsetmatrices61' 'testsetmatrices70' 'testsetmatrices60'});
selecttest.Properties.VariableNames={'matrices61' 'matrices70' 'matrices60'};
[Label,Score] = predict(MdlReduced,table2array(selecttest));
index = Label == winner2(1+split:end);
predaccuracy=sum(index)/length(testsetmatrices) 
