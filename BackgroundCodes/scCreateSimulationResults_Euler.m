rng('shuffle');
scurr = rng;

fileID = fopen(['Eulersim_' num2str(scurr.Seed) '.txt'],'w');


while true
    
    Moff=zeros(4,4);
    Moff(1,1)=.01 + (.9-.01).*rand(1);
    Moff(2,1)=.0001+(.001-.0001).*rand(1);
    Moff(3,1)=.0001+(.001-.0001).*rand(1);
    Moff(2,2)=.01 + (.9-.01).*rand(1);
    Moff(4,2)=.0001+(.001-.0001).*rand(1);
    Moff(3,3)=.01 + (.9-.01).*rand(1);
    Moff(4,3)=.0001+(.001-.0001).*rand(1);
    Moff(4,4)=.01 + (.9-.01).*rand(1);
    MA=Moff;
    MA(1,1)= -.9 + (-Moff(4,4)+.9).*rand(1);
    fA=Moff(1,1)-MA(1,1); %drug effect
    MA(3,3)= MA(3,3)-fA;
    MA(2,1)=Moff(2,1)+(.001-Moff(2,1)).*rand(1);
    MA(4,3)=Moff(4,3)+(.001-Moff(4,3)).*rand(1);
    MB=Moff;
    MB(1,1)= -.9 + (-Moff(4,4)+.9).*rand(1);
    fB=Moff(1,1)-MB(1,1); %drug effect
    MB(2,2)= MB(2,2)-fB;
    MB(3,1)=Moff(3,1)+(.001-Moff(3,1)).*rand(1);
    MB(4,2)=Moff(4,2)+(.001-Moff(4,2)).*rand(1);
    
    n0=[0;0;0;0];
    n0(1)=round(1 + (100-4).*rand(1));
    left=100-n0(1);
    n0(2)=round(0+(left-2).*rand(1));
    left=100-n0(1)-n0(2);
    n0(3)=round(0+(left-1).*rand(1));
    n0(4)=100-n0(1)-n0(2)-n0(3);
    
    if sum(MA*n0)>0 | sum(MB*n0)>0
        continue
    end
    
    
    %make matrix of vectors of matricies
    vector=[];
    for k=1:4
        vector=[vector;Moff(:,k)];
    end
    for k=1:4
        vector=[vector;MA(:,k)];
    end
    for k=1:4
        vector=[vector;MB(:,k)];
    end
    vector=[vector;n0];
    
    N0=n0';
    
    results=[];
    for j=1:50
        [advantagecombo,advantagealternating2,advantageoscillating,advantagealternatingbreak,advantagecombobreak] = test_optimalbreakOde_4celltypes_nocap_Euler( MA,MB,Moff,N0);
        results(j,1)=advantagecombo;
        results(j,2)=advantagealternating2;
        results(j,3)=advantageoscillating;
        results(j,4)=advantagealternatingbreak;
        results(j,5)=advantagecombobreak;
    end
    
    if length(find(results(:,1)==200))>0
        result1=200;
    elseif length(find(results(:,1)==-100))>0
        result1=-100;
    else
        result1=max(results(:,1)); %maxadvantagecombo
    end
    
    if length(find(results(:,2)==200))>0
        result2=200;
    elseif length(find(results(:,2)==-100))>0
        result2=-100;
    else
        result2=max(results(:,2)); %maxadvantagealternating
    end
    
    if length(find(results(:,3)==200))>0
        result3=200;
    elseif length(find(results(:,3)==-100))>0
        result3=-100;
    else
        result3=max(results(:,3)); %maxadvantageoscillating
    end
    
    if length(find(results(:,4)==200))>0
        result4=200;
    elseif length(find(results(:,4)==-100))>0
        result4=-100;
    else
        result4=max(results(:,4)); %maxadvantagealternatingbreaks
    end
    
    if length(find(results(:,5)==200))>0
        result5=200;
    elseif length(find(results(:,5)==-100))>0
        result5=-100;
    else
        result5=max(results(:,5)); %maxadvantagecombobreaks
    end
    
    result6=length(find(results(:,1)>1))/length(results(:,1));
    result7=length(find(results(:,2)>1))/length(results(:,2));
    result8=length(find(results(:,3)>1))/length(results(:,3));
    result9=length(find(results(:,4)>1))/length(results(:,4));
    result10=length(find(results(:,5)>1))/length(results(:,5));
    result11=sum(Moff*n0);
    result12=sum(MA*n0);
    result13=sum(MB*n0);
    
    fprintf(fileID,'%f\t',vector);
    fprintf(fileID,'%f\t', result1);
    fprintf(fileID,'%f\t', result2);
    fprintf(fileID,'%f\t', result3);
    fprintf(fileID,'%f\t', result4);
    fprintf(fileID,'%f\t', result5);
    fprintf(fileID,'%f\t', result6);
    fprintf(fileID,'%f\t', result7);
    fprintf(fileID,'%f\t', result8);
    fprintf(fileID,'%f\t', result9);
    fprintf(fileID,'%f\t', result10);
    fprintf(fileID,'%f\t', result11);
    fprintf(fileID,'%f\t', result12);
    fprintf(fileID,'%f\t', result13);
    fprintf(fileID,'\n');
end
