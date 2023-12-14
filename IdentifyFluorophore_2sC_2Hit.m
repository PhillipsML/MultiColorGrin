function CCT_comb = IdentifyFluorophore_2sC_2Hit(hek, cells_line, cells_n_line,PD)

%% Description
%Will find flourophore matches using the hits over ROI and then also
%corrected to theoretical distributions to id overrepresented
%fluorophores. 

%hek is the hek cell data, 10 x 204
%cells_line are the ROI spectra, n x 204
%cells_n_line are the normalized ROI spectra, n x 204
%PD is the theoretical distribution of beta, 2 x 10. First row is raw,
%second row is normalized

%% housekeeping
%assign variables space
beta = zeros(size(cells_line,1),10);
ota = zeros(size(cells_line,1),10);
ota_n = zeros(size(cells_line,1),10);
CellType = zeros(size(cells_line,1),2);
beta_n = zeros(size(cells_line,1),10);
CellType_n = zeros(size(cells_line,1),2);
CCT_og = zeros(size(cells_line,1),1);
CCT_c = zeros(size(cells_line,1),1);
CCT_comb = zeros(size(cells_line,1),1);

%% regular
for i = 1:size(cells_line)
    beta(i,:) = linfit(hek, cells_line(i,:));
end

cutoff = mean(beta)+1.5*std(beta);
ispos = beta>cutoff;

for i =1:size(cells_line)
    ota(i,:) = abs((beta(i,:)-mean(beta))./std(beta));
    ota(i,:) = ota(i,:).*ispos(i,:);
    if sum(ispos(i,:))<1
        CellType(i,1:4) = 0;
    elseif sum(ispos(i,:))==1
        CellType(i,1) = find(ispos(i,:)==1);
        CellType(i,2) = max(ota(i,:));
        CellType(i,3:4) = 0;
    else
        sota = sort(ota(i,:),'descend');
        CellType(i,1) = find(ota(i,:)==sota(1));
        CellType(i,3) = find(ota(i,:)==sota(2));
        CellType(i,2) = sota(1);
        CellType(i,4) = sota(2);
    end
end

for i = 1:size(cells_n_line)
    beta_n(i,:) = linfit(hek, cells_n_line(i,:));
end

cutoff_n = mean(beta_n)+1.5*std(beta_n);
ispos_n = beta_n>cutoff_n;

for i =1:size(cells_n_line)
    ota_n(i,:) = abs((beta_n(i,:)-mean(beta_n))./std(beta_n));
    ota_n(i,:) = ota_n(i,:).*ispos_n(i,:);
    if sum(ispos_n(i,:))<1
        CellType_n(i,1:4) = 0;
    elseif sum(ispos_n(i,:))==1
        CellType_n(i,1) = find(ispos_n(i,:)==1);
        CellType_n(i,2) = max(ota_n(i,:));
        CellType_n(i,3:4) = 0;
    else
        sota = sort(ota_n(i,:),'descend');
        CellType_n(i,1) = find(ota_n(i,:)==sota(1));
        CellType_n(i,3) = find(ota_n(i,:)==sota(2));
        CellType_n(i,2) = sota(1);
        CellType_n(i,4) = sota(2);
    end
end

CTc = [CellType,CellType_n];
for i = 1:size(cells_line,1)
    os = CTc(i,[2,4,6,8]);fs = CTc(i,[1,3,5,7]);
    if max(os>0)
        sos = sort(os,'descend');
        CCT_og(i,1) = fs(find(os==sos(1)));
        if sos(2)>0
            if fs(find(os==sos(1)))==fs(find(os==sos(2)))
                if sos(3)>0
                    CCT_og(i,2) = fs(find(os==sos(3)));
                else
                    CCT_og(i,2) =NaN;
                end
            else
                CCT_og(i,2) = fs(find(os==sos(2)));
            end
        else
            CCT_og(i,2) = NaN;
        end
    else
        CCT_og(i,1:2) = NaN;
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
    ota(i,:) = abs((beta(i,:)-cutoff_c)./cutoff_c);
    ota(i,:) = ota(i,:).*ispos(i,:);
    if sum(ispos(i,:))<1
        CellType(i,1:4) = 0;
    elseif sum(ispos(i,:))==1
        CellType(i,1) = find(ispos(i,:)==1);
        CellType(i,2) = max(ota(i,:));
        CellType(i,3:4) = 0;
    else
        sota = sort(ota(i,:),'descend');
        CellType(i,1) = find(ota(i,:)==sota(1));
        CellType(i,3) = find(ota(i,:)==sota(2));
        CellType(i,2) = sota(1);
        CellType(i,4) = sota(2);
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
    ota_n(i,:) = abs((beta_n(i,:)-cutoff_cn)./cutoff_cn);
    ota_n(i,:) = ota_n(i,:).*ispos_n(i,:);
    if sum(ispos_n(i,:))<1
        CellType_n(i,1:4) = 0;
    elseif sum(ispos_n(i,:))==1
        CellType_n(i,1) = find(ispos_n(i,:)==1);
        CellType_n(i,2) = max(ota_n(i,:));
        CellType_n(i,3:4) = 0;
    else
        sota = sort(ota_n(i,:),'descend');
        CellType_n(i,1) = find(ota_n(i,:)==sota(1));
        CellType_n(i,3) = find(ota_n(i,:)==sota(2));
        CellType_n(i,2) = sota(1);
        CellType_n(i,4) = sota(2);
    end
end


CTc_c = [CellType,CellType_n];
for i = 1:size(cells_line,1)
    os = CTc_c(i,[2,4,6,8]);fs = CTc_c(i,[1,3,5,7]);
    if max(os>0)
        sos = sort(os,'descend');
        CCT_c(i,1) = fs(find(os==sos(1)));
        if sos(2)>0
            if fs(find(os==sos(1)))==fs(find(os==sos(2)))
                if sos(3)>0
                    CCT_c(i,2) = fs(find(os==sos(3)));
                else
                    CCT_c(i,2) =NaN;
                end
            else
                CCT_c(i,2) = fs(find(os==sos(2)));
            end
        else
            CCT_c(i,2) = NaN;
        end
    else
        CCT_c(i,1:2) = NaN;
    end
end

%% compare

for i =1:size(CCT_og,1)
    if isnan(CCT_og(i,1))
        CCT_comb(i,1:2)=CCT_c(i,1:2);
    else
        CCT_comb(i,1) = CCT_og(i,1);
        if isnan(CCT_og(i,2))
            if CCT_og(i,1)==CCT_c(i,1)
                CCT_comb(i,2)= CCT_c(i,2);
            else
                CCT_comb(i,2)= CCT_c(i,1);
            end
        else
            CCT_comb(i,2)=CCT_og(i,2);
        end

    end
end

           
