rng('shuffle');
scurr = rng;
alldata=[];

while true
    
    Moff=zeros(2,2);
    Moff(1,1)=.01 + (.9-.01).*rand(1);
    Moff(2,1)=.0001+(.001-.0001).*rand(1); 
    Moff(2,2)=.01 + (.9-.01).*rand(1);
    
    M=Moff;
    M(1,1)= -.9 + (-Moff(2,2)+.9).*rand(1); 
    M(2,1)=Moff(2,1)+(.001-Moff(2,1)).*rand(1); 
    
    n0=[0;0];
    n0(1)=1 + (100-1).*rand(1);
    n0(2)=100-n0(1);
    
    if sum(M*n0)>0
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
    vector=[vector;n0];
    
    results=[];
    for j=1:50
        [ breakadvantage,~,~,~,~,~] = test_singledrug_optimalbreakOde_nocap( M,Moff,n0);
        results(j,1)=breakadvantage;
    end
    result1=max(results(:,1));
    result2=length(find(results>1))/length(results);
    result3=sum(Moff*n0);
    result4=sum(M*n0);

    data=[vector',result1,result2,result3,result4];
    alldata=[alldata;data];
end


