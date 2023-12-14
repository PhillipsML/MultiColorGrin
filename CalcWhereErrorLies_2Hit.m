function [Corr, CNan] = CalcWhereErrorLies_2Hit(Cheat,CCT_c)



for j = 1:10
    clear id nid
    [id(:,1),id(:,2)] = find(Cheat==j);
    nid = id(isnan(Cheat(id(:,1),2)),1);
    Corr(j,1,1) = sum(CCT_c(nid,1)==j)/(2*size(nid,1));
    Corr(j,1,2) = sum(isnan(CCT_c(nid,2)))/(2*size(nid,1));
    Corr(j,1,3) = sum(sum(CCT_c(nid,:)==4))/(2*size(nid,1));
    CNan(j,1,1) = 0.5-Corr(j,1,1);
    CNan(j,1,2) = 1-Corr(j,1,1)-Corr(j,1,2)-Corr(j,1,3) -CNan(j,1,1);
    for k = 1:10
        if j==k
            Corr(j,k+1,1)=NaN; 
            Corr(j,k+1,2)=NaN;
            Corr(j,k+1,3)=NaN;
            CNan(j,k+1,1)=NaN;
            CNan(j,k+1,2)=NaN;
        else
            clear kd
            kd = id(find(sum(Cheat(id(:,1),:)==k,2)),1);
            Corr(j,k+1,1) = sum(sum(CCT_c(kd,:)==j))/(2*size(kd,1));
            Corr(j,k+1,2) = sum(sum(CCT_c(kd,:)==k))/(2*size(kd,1));
            Corr(j,k+1,3) = sum(sum(CCT_c(kd,:)==4))/(2*size(kd,1));
            CNan(j,k+1,1) = sum(sum(isnan(CCT_c(kd,:))))/(size(kd,1)*2);
            CNan(j,k+1,2) = 1-Corr(j,k+1,1)-Corr(j,k+1,2)-Corr(j,k+1,3)-CNan(j,k+1,1);
        end
    end
    figure(j)
    %bar([Corr(j,:,1)',Corr(j,:,2)',[Corr(j,:,3)+CNan(j,:,1)]',CNan(j,:,2)'],'stacked')
    bar([Corr(j,:,1)',[Corr(j,:,3)+CNan(j,:,1)]',CNan(j,:,2)',Corr(j,:,2)'],'stacked')
    ylim([0 1])
end
