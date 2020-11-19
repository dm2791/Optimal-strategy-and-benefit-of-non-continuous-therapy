%% load data
tumormarker_data=readtable('tumormarkers_paper.csv');
treatment_data=readtable('treatmentdata_paper.csv');

%% plot
patients=unique(treatment_data.Sample_ID);

for i=1:length(patients)
    index=strcmp(tumormarker_data.Sample_ID,patients(i));
    index2=find(strcmp(treatment_data.Sample_ID,patients(i)));
    figure,gscatter(tumormarker_data.relativedays2(index),tumormarker_data.CA_15_3num(index))
    maxy=ylim;
    maxy=maxy(2);
    hold on
    for j=index2(1):index2(end)
        if rem(j,2)==1
    plot([treatment_data.Start_Time(j),treatment_data.End_Time(j)],[maxy/2,maxy/2])
    text(treatment_data.Start_Time(j),(maxy/2)*1.1,treatment_data.Tx_Regimen(j))
        else
                plot([treatment_data.Start_Time(j),treatment_data.End_Time(j)],[maxy/1.5,maxy/1.5])
    text(treatment_data.Start_Time(j),(maxy/1.5)*1.1,treatment_data.Tx_Regimen(j))
        end
    end
    seq=strcmp(treatment_data.Sample_Tx_Status(index2),'Post-treatment');
    seq_date=index2(find(seq,1,'last'));
    plot([treatment_data.End_Time(seq_date),treatment_data.End_Time(seq_date)],[0,maxy],'k--')
    hold off
    legend([{'CA.15.3'};treatment_data.Tx_Regimen(index2)],'Location','northeastoutside')
    title(patients(i))
    xlabel('days')
    ylabel('U/mL')
end

%% calculate slope at t=0 for each treatment per patient as well as time to resistance (TTR)
% for TTR, estimate time to return to starting size for purposes of model

%patient 468 
Let_468=(tumormarker_data.CA_15_3num(206)-tumormarker_data.CA_15_3num(205))/(tumormarker_data.relativedays2(206)-tumormarker_data.relativedays2(205));
TTR_Let468=(tumormarker_data.relativedays2(212)+tumormarker_data.relativedays2(213))/2-tumormarker_data.relativedays2(205);
NabPac_468=(tumormarker_data.CA_15_3num(217)-tumormarker_data.CA_15_3num(216))/(tumormarker_data.relativedays2(217)-tumormarker_data.relativedays2(216));
TTR_NabPac468=tumormarker_data.relativedays2(222)-tumormarker_data.relativedays2(216);

%patient 288
Let_288=(tumormarker_data.CA_15_3num(52)-tumormarker_data.CA_15_3num(51))/(tumormarker_data.relativedays2(52)-tumormarker_data.relativedays2(51));
TTR_Let288=tumormarker_data.relativedays2(103)-tumormarker_data.relativedays2(51);
Everol_288=(tumormarker_data.CA_15_3num(105)-tumormarker_data.CA_15_3num(104))/(tumormarker_data.relativedays2(105)-tumormarker_data.relativedays2(104));
TTR_Everol288=(tumormarker_data.relativedays2(111)+tumormarker_data.relativedays2(112))/2-tumormarker_data.relativedays2(104);

%patient 15
Gem_15=(tumormarker_data.CA_15_3num(12)-tumormarker_data.CA_15_3num(11))/(tumormarker_data.relativedays2(12)-tumormarker_data.relativedays2(11));
TTR_Gem15=tumormarker_data.relativedays2(28)-tumormarker_data.relativedays2(11);

%patient 422
Cap_422=(tumormarker_data.CA_15_3num(153)-tumormarker_data.CA_15_3num(152))/(tumormarker_data.relativedays2(153)-tumormarker_data.relativedays2(152));
TTR_Cap422=tumormarker_data.relativedays2(167)-tumormarker_data.relativedays2(152);
Erib_422=(tumormarker_data.CA_15_3num(181)-tumormarker_data.CA_15_3num(180))/(tumormarker_data.relativedays2(181)-tumormarker_data.relativedays2(180));
TTR_Erib422=tumormarker_data.relativedays2(190)-tumormarker_data.relativedays2(180);
Pac_422=(tumormarker_data.CA_15_3num(191)-tumormarker_data.CA_15_3num(190))/(tumormarker_data.relativedays2(191)-tumormarker_data.relativedays2(190));
TTR_Pac422=tumormarker_data.relativedays2(199)-tumormarker_data.relativedays2(190);

%% plot on values for all patients
ons=[Let_468;NabPac_468;Let_288;Everol_288;Gem_15;Cap_422;Erib_422;Pac_422];
onsgroup=[{'Let_468'},{'NabPac_468'},{'Let_288'},{'Everol_288'},{'Gem_15'},{'Cap_422'},{'Erib_422'},{'Pac_422'}];
figure,gscatter([1:8],ons,onsgroup')
ylabel('initial slope on treatment')
%% load all clinical data
treatmentTbl=readtable('treatmenttable_paper.csv');
%% pie chart
stop_reason = tabulate(treatmentTbl.Reason_Stop);
stop_reason(6,2)={96}; %combine "death" with "progression/death"
stop_reason(12,:)=[];
stop_reason(11,2)={222}; %combine very small slices: unknown, postmenopause,break,other
stop_reason(8:10,:)=[];
stop_reason(3,2)={1079}; %combine ongoing and Ongoing
stop_reason(9,:)=[];
%redo frequency calc
stopreason_count=cell2mat(stop_reason(:,2));
stopreason_freq=stopreason_count/sum(stopreason_count);
stop_reason(:,4)=num2cell(stopreason_freq*100);
labels = stop_reason(:,1);
figure,pie(cell2mat(stop_reason(:,2)),labels)

%% average time to resistance
resistanceTbl=treatmentTbl(strcmp(treatmentTbl.Reason_Stop,'Disease Progression'),:);
TTR=resistanceTbl.End_Time-resistanceTbl.Start_Time;

% plot as Kaplan meier curve
figure,ecdf(round(TTR,-2),'function','survivor','Bounds','on')
xlabel('days')
ylabel('proportion of patients without resistance')
%% how many patients got resistance
nTreatment = length(unique(treatmentTbl.Patient_ID));
patients=unique(treatmentTbl.Patient_ID);
resistance_patients=unique(resistanceTbl.Patient_ID);
labels_res={'No Resistance','Resistance'};

% how many patients had resistance to more than 1 drug
numres=tabulate(resistanceTbl.Patient_ID);
numres_count=cell2mat(numres(:,2));
numres_one=find(numres_count==1);
numres_two=find(numres_count==2);
numres_multi=find(numres_count>2);
labels_numres={'One','Two','Three or More'};

figure
ax1 = subplot(1,2,1);
pie([nTreatment-length(resistance_patients),length(resistance_patients)],labels_res)
title(ax1,'#Patients with Resistance');

ax2 = subplot(1,2,2);
pie(ax2,[length(numres_one),length(numres_two),length(numres_multi)],labels_numres);
title(ax2,'#Resistant Drugs');