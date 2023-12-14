function CCT_comb = IdentifyFluorophore_2StepCorr(hek, cells_line, cells_n_line,PD)
%%regular
for i = 1:size(cells_line)
    beta(i,:) = linfit(hek, cells_line(i,:));
end

cutoff = mean(beta)+1.5*std(beta);
ispos = beta>cutoff;

for i =1:size(cells_line)
    clear ota
    if sum(ispos(i,:))==1
        ota = (beta(i,:)-mean(beta))./std(beta);
        CellType(i,1) = find(ispos(i,:)==1);
        CellType(i,2) = max(ota);
    elseif sum(ispos(i,:))>1
        ota = (beta(i,:)-mean(beta))./std(beta);
        CellType(i,1) = find(ota==max(ota));
        CellType(i,2) = max(ota);
    else
        CellType(i,1:2) = NaN;
    end
end

for i = 1:size(cells_n_line)
    beta_n(i,:) = linfit(hek, cells_n_line(i,:));
end

cutoff_n = mean(beta_n)+1.5*std(beta_n);
ispos_n = beta_n>cutoff_n;

for i =1:size(cells_n_line)
    clear ota
    if sum(ispos_n(i,:))==1
        ota = (beta_n(i,:)-mean(beta_n))./std(beta_n);
        CellType_n(i,1) = find(ispos_n(i,:)==1);
        CellType_n(i,2) = max(ota);
    elseif sum(ispos_n(i,:))>1
        ota = (beta_n(i,:)-mean(beta_n))./std(beta_n);
        CellType_n(i,1) = find(ota==max(ota));
        CellType_n(i,2) = max(ota);
    else
        CellType_n(i,1:2) = NaN;
    end
end


for i = 1:size(CellType,1)
    if isnan(CellType(i,1))
        if isnan(CellType_n(i,1))
            CCT_og(i,1) = NaN;
        else
            CCT_og(i,1) =CellType_n(i,1);
        end
    else
        if isnan(CellType_n(i,1))
            CCT_og(i,1) = CellType(i,1);
        else
            if CellType(i,2)>CellType_n(i,2)
                 CCT_og(i,1) = CellType(i,1);
            else
                 CCT_og(i,1) = CellType_n(i,1);
            end
        end
    end
end

%% Correction
for i = 1:size(cells_line)
    beta(i,:) = linfit(hek, cells_line(i,:));
end

dTheor = abs((mean(beta)/sum(mean(beta)))./PD(1,:));
nstd = 1.5./dTheor;
cutoff_c = mean(beta)+nstd.*std(beta);
ispos = beta>cutoff_c;

for i =1:size(cells_line)
    clear ota
    if sum(ispos(i,:))==1
        ota = (beta(i,:)-cutoff_c)./cutoff_c;
        CellType(i,1) = find(ispos(i,:)==1);
        CellType(i,2) = max(ota);
    elseif sum(ispos(i,:))>1
        ota = (beta(i,:)-cutoff_c)./cutoff_c;
        CellType(i,1) = find(ota==max(ota));
        CellType(i,2) = max(ota);
    else
        CellType(i,1:2) = NaN;
    end
end

for i = 1:size(cells_n_line)
    beta_n(i,:) = linfit(hek, cells_n_line(i,:));
end

dTheor_n = abs((mean(beta_n)/sum(mean(beta_n)))./PD(2,:));
nstd_n = 1.5./dTheor_n;
cutoff_cn = mean(beta_n)+nstd_n.*std(beta_n);
ispos_n = beta_n>cutoff_cn;

for i =1:size(cells_n_line)
    clear ota
    if sum(ispos_n(i,:))==1
        ota = (beta_n(i,:)-cutoff_cn)./cutoff_cn;
        CellType_n(i,1) = find(ispos_n(i,:)==1);
        CellType_n(i,2) = max(ota);
    elseif sum(ispos_n(i,:))>1
        ota = (beta_n(i,:)-cutoff_cn)./cutoff_cn;
        CellType_n(i,1) = find(ota==max(ota));
        CellType_n(i,2) = max(ota);
    else
        CellType_n(i,1:2) = NaN;
    end
end


for i = 1:size(CellType,1)
    if isnan(CellType(i,1))
        if isnan(CellType_n(i,1))
            CCT_c(i,1) = NaN;
        else
            CCT_c(i,1) =CellType_n(i,1);
        end
    else
        if isnan(CellType_n(i,1))
            CCT_c(i,1) = CellType(i,1);
        else
            if CellType(i,2)>CellType_n(i,2)
                 CCT_c(i,1) = CellType(i,1);
            else
                 CCT_c(i,1) = CellType_n(i,1);
            end
        end
    end
end
%% compare
for i =1:size(CCT_og,1)
    if isnan(CCT_og(i))
           CCT_comb(i,1)=CCT_c(i,1);
    else
        CCT_comb(i,1)=CCT_og(i,1);
    end
end

           
