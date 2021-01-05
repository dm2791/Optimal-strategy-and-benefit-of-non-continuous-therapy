rng('shuffle');
scurr = rng;

fileID = fopen(['Euler1drugsim_' num2str(scurr.Seed) '.txt'],'w');

while true
    
    Moff=zeros(2,2);
    Moff(1,1)=.01 + (.9-.01).*rand(1);
    Moff(2,1)=.0001+(.001-.0001).*rand(1);
    Moff(2,2)=.01 + (.9-.01).*rand(1);
    M=Moff;
    M(1,1)= -.9 + (-Moff(2,2)+.9).*rand(1);
    M(2,1)=Moff(2,1)+(.001-Moff(2,1)).*rand(1);
    N0=[0,0];
    N0(1)=1 + (100-1).*rand(1);
    N0(2)=100-N0(1);
    
    if sum(M*N0')>0
        continue
    end
    
    %make matrix of vectors of matricies
    vector=[];
    for k=1:2
        vector=[vector;Moff(:,k)];
    end
    for k=1:2
        vector=[vector;M(:,k)];
    end
    vector=[vector;N0'];
    
    results=[];
    
    for j=1:50
        [ breakadvantage,~,~,~,~] = test_singledrug_optimalbreak_Euler_nocap( M,Moff,N0);
        results(j,1)=breakadvantage;
    end
    if length(find(results==200))>0
        result1=200;
    elseif length(find(results==-100))>0
        result1=-100;
    elseif length(find(results==100))>0
        result1=100;
    else
        result1=max(results(:,1));
    end
    result2=length(find(results>1))/length(results);
    result3=sum(Moff*N0');
    result4=sum(M*N0');
    
    
    fprintf(fileID,'%f\t',vector);
    fprintf(fileID,'%f\t', result1);
    fprintf(fileID,'%f\t', result2);
    fprintf(fileID,'%f\t', result3);
    fprintf(fileID,'%f\t', result4);
    fprintf(fileID,'\n');
end


