%% Load data
alldatamat=readmatrix('Euler1drugsim.csv');

%% Unsupervised analysis
% remove rows where advantage=-100 (continuous went to zero)
advantages=alldatamat(:,11);
advantagescut=advantages;
advantagescut(find(advantages==-100))=[];
alldatamat(find(advantages==-100),:)=[];

matrices=alldatamat(:,1:10);
[coeff,score,latent,tsquared,explained,mu] = pca(matrices);

conditions=ones(size(advantagescut,1),1);
conditionsadd=advantagescut(:,1);
conditions=[conditions,conditionsadd];
logc=log(conditionsadd);

%color on degree of response
color=[];
for i=1:size(conditions,1)
    if conditions(i,2)==200 
        color(i,:)=[.5,0,.5]; %dark magenta
    elseif conditions(i,2)>1 
        color(i,:)=[1,0,1]; %magenta
    else
        color(i,:)=[0,0,0];
    end
end
colord=color.*(logc/max(logc));
colord(find(conditionsadd==200),:)=repmat([1,0,1],length(find(conditionsadd==200)),1);
figure
scatter(score(:,1),score(:,2),50,colord,'filled')

%plot pie chart
curedrate=length(find(conditionsadd==200))/length(conditionsadd);
figure,pie([curedrate,1-curedrate])

%% Machine learning
matrices(:,11)=matrices(:,10)./matrices(:,9); %n02/n01
matrices(:,12)=matrices(:,6)./matrices(:,2); %u2/u1
matrices(:,13)=matrices(:,1)./matrices(:,4); %a1/b
matrices(:,14)=matrices(:,5)./matrices(:,4); %a2/b
matrices(:,15)=alldatamat(:,13); %sum(Moff*n0)
matrices(:,16)=alldatamat(:,14); %sum(M*n0)
matrices(:,17)=matrices(:,1)-matrices(:,5); %f

testsetmatrices=matrices(15001:end,:);
matrices=matrices(1:15000,:);
testsetlogc=logc(15001:end);
logc=logc(1:15000);

%% figure out which tumor properties separate cure from relapse delay
%binarize
table=array2table(matrices);
tableselect=array2table(matrices(:,11:17));
binarylogc=ones(length(logc),1);
binarylogc(find(exp(logc)>199& exp(logc)<201))=0;
binarylogctest=ones(length(testsetlogc),1);
binarylogctest(find(exp(testsetlogc)>199& exp(testsetlogc)<201))=0;

Mdl = fitcecoc(table2array(tableselect),binarylogc);
isLoss = resubLoss(Mdl) 
[Label,Score] = predict(Mdl,testsetmatrices(:,11:17));
index = Label == binarylogctest;
predaccuracy=sum(index)/length(testsetmatrices) %.97

[ranks,weights] = relieff(table2array(tableselect),binarylogc,5);
figure, bar(weights(ranks))

%% Load data 2 drugs simulation
alldatamat=readmatrix('Euler2drugsim.csv');
%% Unsupervised analysis
% remove rows of zeros
advantages=alldatamat(:,53:57);
advantagescut=advantages;
advantagescut( all(~advantages,2), : ) = [];
alldatamat(all(~advantages,2),:)=[];

% remove rows where advantage=-100 in any iteration(continuous went to zero)
[row,col] = find(advantagescut(:,:)==-100);
advantagescut(row,:)=[];
alldatamat(row,:)=[];

matrices=alldatamat(:,1:52);
[coeff,score,latent,tsquared,explained,mu] = pca(matrices);
conditions=ones(size(advantagescut,1),1);
conditionsadd=advantagescut(:,1:5);
conditions=[conditions,conditionsadd];

conditions(:,4)=[];
color=[];
combo=0;
comboc=0; %cured count
combobreaks=0;
combobreaksc=0;
alternatingbreaks=0;
alternatingbreaksc=0;
alternating=0;
alternatingc=0;
soc=0;
winner=[];
for i=1:size(conditions,1)
    if conditions(i,2)==200
        color(i,:)=[0,1,1]; % cyan if combo cured
        comboc=comboc+1;
        winner(i)=2;
    elseif conditions(i,3)==200
        color(i,:)=[.5,0,0]; % dark red if alternating cured
        alternatingc=alternatingc+1;
        winner(i)=3;
    elseif conditions(i,4)==200
        color(i,:)=[.5,0.,5]; % dark magenta if alternatingbreaks cured
        alternatingbreaksc=alternatingbreaksc+1;
        winner(i)=4;
    elseif conditions(i,5)==200
        color(i,:)=[0,0,.5]; % dark blue if combobreaks cured
        combobreaksc=combobreaksc+1;
        winner(i)=5;
    elseif find(conditions(i,:)==max(conditions(i,:)))==2;
        color(i,:)=[0.5843,0.8157,0.9882]; %light blue if combo best
        combo=combo+1;
        winner(i)=2;
    elseif find(conditions(i,:)==max(conditions(i,:)))==5;
        color(i,:)=[0,0,1]; %blue if combobreaksbest
        combobreaks=combobreaks+1;
        winner(i)=5;
    elseif find(conditions(i,:)==max(conditions(i,:)))==4;
        color(i,:)=[1,0,1]; %magenta if alternating breaks best
        alternatingbreaks=alternatingbreaks+1;
        winner(i)=4;
    elseif find(conditions(i,:)==max(conditions(i,:)))==3;
        color(i,:)=[1,0,0]; %grouped colors: red if alternating2 best
        alternating=alternating+1;
        winner(i)=3;
    elseif find(conditions(i,:)==max(conditions(i,:)))==1;
        color(i,:)=[0,0,0];
        soc=soc+1;
        winner(i)=1;
    else
        winner=100;
    end
end


figure
scatter3(score(:,1),score(:,2),score(:,3),50,color,'filled')

%calculate and plot frequencies of strategy success
myColors=[repmat([0,0,0],2,1);repmat([0,1,1],2,1);repmat([.5,0,0],2,1);repmat([.5,0.,5],2,1);...
    repmat([0,0,.5],2,1);repmat([0.5843,0.8157,0.9882],2,1);repmat([1,0,0],2,1);repmat([1,0,1],2,1);repmat([0,0,1],2,1)];
for k=1:9
    b(k).FaceColor = 'flat';
    b(k).CData = myColors(k*2-1:k*2,:);
end
winnercount=tabulate(winner);
figure,pie(winnercount(:,3))
%% plot cure distribution
curerate=[length(find(conditions(:,2)==200)),length(find(conditions(:,3)==200)),length(find(conditions(:,4)==200)),length(find(conditions(:,5)==200))];
cureratepercent=curerate./length(find(conditions==200));
figure,pie(curerate)

%% biplot
logc=log(conditions);
[coeff,score,latent,tsquared,explained,mu] = pca(conditions);
figure
biplot(coeff(:,1:2),'Scores',score(:,1:2));