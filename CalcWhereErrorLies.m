function Error_c = CalcWhereErrorLies(Cheat,CCT)
%Will find where the errors are for the modeling dataset. This works for
%single hits, not for the dual hit.

%Cheat is the cheat sheat for what the matches should be
%CCT is what the model spit out


marks =[find(diff(Cheat));250];
Corr = CCT==Cheat;
CCT_nan = isnan(CCT);
CCT_gcamp = CCT==4;

for i = 1:10
    if i ==1
        Error_c(i,1) = mean(sum(Corr(1:marks(i),:))/sum(Cheat==i));
        Error_c(i,2) = mean(sum(CCT_nan(1:marks(i),:))/sum(Cheat==i));
        Error_c(i,3) = mean(sum(CCT_gcamp(1:marks(i),:))/sum(Cheat==i));
        Error_c(i,4) = 1-sum(Error_c(i,1:3));
    else
        Error_c(i,1) = mean(sum(Corr(marks(i-1)+1:marks(i),:))/sum(Cheat==i));
        Error_c(i,2) = mean(sum(CCT_nan(marks(i-1)+1:marks(i),:))/sum(Cheat==i));
        Error_c(i,3) = mean(sum(CCT_gcamp(marks(i-1)+1:marks(i),:))/sum(Cheat==i));
        Error_c(i,4) = 1-sum(Error_c(i,1:3));
    end
end