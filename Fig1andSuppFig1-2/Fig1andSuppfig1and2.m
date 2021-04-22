%% load data
alldatamat=readmatrix('1drugsimulationdata_paper_rev.csv');
advantagescut=alldatamat(:,11);
%% Unsupervised analysis
matrices=alldatamat(:,1:10);
[coeff,score,latent,tsquared,explained,mu] = pca(matrices);

conditions=ones(size(advantagescut,1),1);
conditionsadd=advantagescut(:,1);
conditions=[conditions,conditionsadd];
logc=log(conditionsadd);

color=[];
pulse=0;
for i=1:size(conditions,1)
    if conditions(i,2)>1 
        color(i,:)=[1,0,1]; %magenta
        pulse=pulse+1;
    else
        color(i,:)=[0,0,0]; %black
    end
end

%plot
figure
scatter(score(:,1),score(:,2),50,color,'filled')

%color on degree of response
colord=color.*(logc/max(logc));
figure
scatter(score(:,1),score(:,2),50,colord,'filled')
colormap(colord)
xlabel('PC1')
ylabel('PC2')

%histogram
figure, histogram(log(conditionsadd))
y=logspace(0,6,20);
xlabel('log(pulse advantage)')
ylabel('%simulated tumors')

%% Make matrix

matrices=alldatamat(:,1:10);

conditions=ones(size(advantagescut,1),1);
conditionsadd=advantagescut(:,1);
conditions=[conditions,conditionsadd];
logc=log(conditionsadd);

matrices(:,11)=matrices(:,10)./matrices(:,9); %n02/n01
matrices(:,12)=matrices(:,6)./matrices(:,2); %u2/u1
matrices(:,13)=matrices(:,1)./matrices(:,4); %g1/g2 off
matrices(:,14)=matrices(:,5)./matrices(:,4); %g1/g2 on
matrices(:,15)=alldatamat(:,13); %sum(Moff*n0)
matrices(:,16)=alldatamat(:,14); %sum(M*n0)
matrices(:,17)=matrices(:,1)-matrices(:,5); %f
%%
testsetmatrices=matrices(40001:end,:);
matrices=matrices(1:40000,:);
testsetlogc=logc(40001:end);
logc=logc(1:40000);

%% plot
C =[];
x=matrices(:,13);
y=logc;
R=.5;
for i = 1:length(x)
    v = (abs(x-x(i)) < R) & (abs(y-y(i)) < R);    
    C(i) = sum(v)-1;                                          
end
figure,scatter(x,y,20,C)
xlabel('g1off/g2off')
ylabel('log(pulse advantage)')

C =[];
x=log(matrices(:,15));
y=logc;
R=.5;
for i = 1:length(x)
    v = (abs(x-x(i)) < R) & (abs(y-y(i)) < R);  
    C(i) = sum(v)-1;                                          
end
figure,scatter(x,y,20,C)
xlabel('log(sum(Moff*n0))')
ylabel('log(pulse advantage)')

C =[];
x=matrices(:,16);
y=logc;
R=.5;
for i = 1:length(x)
    v = (abs(x-x(i)) < R) & (abs(y-y(i)) < R);    
    C(i) = sum(v)-1;                                         
end
figure,scatter(x,y,20,C)
xlabel('sum(M*n0)')
ylabel('log(pulse advantage)')

%% machine learning
table=array2table(matrices);
tableselect=array2table(matrices(:,11:17));
t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature');
rng(1); 
Mdl = fitrensemble(tableselect,logc,'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2 %.91

impOOB = oobPermutedPredictorImportance(Mdl);
figure
bar(sort(impOOB))
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = table.Properties.VariableNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';   

MdlReduced = fitrensemble(table(:,{'matrices13' 'matrices15' 'matrices16'}),logc,'Method','Bag', ...
    'NumLearningCycles',200,'Learners',t);
yHatReduced = oobPredict(MdlReduced);
r2Reduced = corr(Mdl.Y,yHatReduced)^2 %.9, 

MdlReduced2 = fitrensemble(table(:,{'matrices15' 'matrices16'}),logc,'Method','Bag', ...
    'NumLearningCycles',200,'Learners',t);
yHatReduced2 = oobPredict(MdlReduced2);
r2Reduced2 = corr(Mdl.Y,yHatReduced2)^2 % .8 

%% repeat importance predictor with normalized matrices
meanmat=mean(matrices,1);
matricesnorm=matrices(:,11:end)./meanmat(11:end);
tableselect=array2table(matricesnorm);
Mdl = fitrensemble(tableselect,logc,'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2 %.91

impOOB = oobPermutedPredictorImportance(Mdl);
figure
bar(sort(impOOB))
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = table.Properties.VariableNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';   %same as non-normalized
%% testset
tableselecttest=array2table(testsetmatrices);
selecttest=tableselecttest(:,{'testsetmatrices13' 'testsetmatrices15' 'testsetmatrices16'});
selecttest.Properties.VariableNames={'matrices13' 'matrices15' 'matrices16'};
Yfit = predict(MdlReduced,selecttest);
mdltest=fitlm(Yfit,testsetlogc)  

selecttest2=tableselecttest(:,{'testsetmatrices15' 'testsetmatrices16'});
selecttest2.Properties.VariableNames={'matrices15' 'matrices16'};
Yfit2 = predict(MdlReduced2,selecttest2);
mdltest2=fitlm(Yfit2,testsetlogc)
%% figure out differences between tumors that require matrices13 for prediction vs ones that don't. Also which tumors aren't predicted well and why
%using all 3 predictors: which ones have difference greater than root mean squared error of .132?

largediff=find(abs(Yfit-testsetlogc)>.132);
smalldiff=find(abs(Yfit-testsetlogc)<=.132);
selectlargediff=testsetmatrices(largediff,:);
selectsmalldiff=testsetmatrices(smalldiff,:);
[coeff,score,latent,tsquared,explained,mu] = pca(testsetmatrices);
color=repmat([0 0 0],size(testsetlogc,1),1);
color(largediff,:)=repmat([1,0,1],length(largediff),1);
figure,scatter(score(:,1),score(:,2),50,color,'filled')
xlabel('PC1')
ylabel('PC2')
mean(selectlargediff,1);
mean(selectsmalldiff,1);
figure, h1=histogram(selectsmalldiff(:,16));
hold on
h2=histogram(selectlargediff(:,16));
xlabel('sum(M*n0)')
ylabel('#simulated tumors')

color=repmat([0 .45 .74],size(testsetlogc,1),1);
color(largediff,:)=repmat([0,1,1],length(largediff),1);
figure, scatter(Yfit,testsetlogc,50,color)
xlabel('prediction')
ylabel('value')

%do plsr to figure out main differences between big and small diff
X=[selectlargediff;selectsmalldiff];
Y=[repmat(0,1222,1);repmat(1,4399,1)];
[Xloadings,Yloadings,Xscores,Yscores,betaPLS10,PLSPctVar,mse,stats] = plsregress(X,Y,10);

% what are the differences between tumors that are equally predicted with only 15 and 16 vs those that need 13?
%figure, scatter(Yfit2,testsetlogc)
smalldiff2=find(abs(Yfit2-testsetlogc)<=.132); %accurately predicted with same cutoff with only 2 predictors
needs13 = setdiff(smalldiff,smalldiff2); %tumors that need 13 for good prediction
not_needs13=intersect(smalldiff,smalldiff2);
selectneeds13=testsetmatrices(needs13,:);
selectnot=testsetmatrices(not_needs13,:);
color=repmat([0 0 0],size(testsetlogc,1),1);
color(needs13,:)=repmat([1,0,1],length(needs13),1);
figure,scatter(score(:,1),score(:,2),50,color,'filled')
xlabel('PC1')
ylabel('PC2')
mean(selectneeds13,1);
mean(selectnot,1);
figure, h1=histogram(selectnot(:,16));
hold on
h2=histogram(selectneeds13(:,16));
xlabel('sum(M*n0)')
ylabel('#simulated tumors')

color=repmat([0 .45 .74],size(testsetlogc,1),1);
color(needs13,:)=repmat([0,1,1],length(needs13),1);
figure, scatter(Yfit2,testsetlogc,50,color)
xlabel('prediction')
ylabel('value')

%% binarize method
figure, scatter([1:40000],sort(logc)) %inflection point around .7712, exp(.7712)=2.16
xlabel('simulated tumor #')
ylabel('log(pulse advantage)')
binarylogc=ones(length(logc),1);
binarylogc(find(logc<=.7712))=0;
binarylogctest=ones(length(testsetlogc),1);
binarylogctest(find(testsetlogc<=.7712))=0;

Mdlall=fitcecoc(table2array(table),binarylogc);
isLoss = resubLoss(Mdlall) 
[Label,Score] = predict(Mdlall,table2array(tableselecttest));
index = Label == binarylogctest;
predaccuracy=sum(index)/length(testsetmatrices)
[X,Y,T,AUC] = perfcurve(binarylogctest,Score(:,2),1);
figure, plot(X,Y)
hold on

MdlReduced=fitcecoc(table2array(table(:,{'matrices15' 'matrices16'})),binarylogc);
isLoss = resubLoss(MdlReduced) 
[Label,Score] = predict(MdlReduced,table2array(selecttest(:,2:3)));
index = Label == binarylogctest;
predaccuracy=sum(index)/length(testsetmatrices) 
[X,Y,T,AUC] = perfcurve(binarylogctest,Score(:,2),1);
plot(X,Y)
hold off
xlabel('false positive rate')
ylabel('true positive rate')
%% correlate degree of response with percentage of success
figure,hist(sortrows(alldatamat(1:40000,12)))
liklihood_highsuccess=length(find(alldatamat(:,12)>.8))/length(alldatamat); 
xlabel('#pulse strategies better than continuous/#pulse strategies attempted')
ylabel('#simulated tumors')
%% best correlate with pulse advantage for a given set of parameters
M =[-0.0300 ,0; 0.0010, 0.0200];
Moff =[.02 ,0; 0.0010, 0.0200];
n0=[10;0];

results=[];

h = waitbar(0, 'simulation running');
for i=1:100
    waitbar(i/100, h);
[ breakadvantage,x,timeon,timeoff,totaltimeon,totaltimeoff] = test_singledrug_optimalbreakOde_nocap( M,Moff,n0);
    results(i,1)=breakadvantage;
    results(i,2)=x;
    results(i,3)=timeon;
    results(i,4)=timeoff;
    results(i,5)=totaltimeon;
    results(i,6)=totaltimeoff;
end
close(h);
figure, scatter(results(:,4),results(:,1))
xlabel('timeon(x-1)')
ylabel('pulse advantage')
mdl=fitlm(results(:,4),results(:,1))
