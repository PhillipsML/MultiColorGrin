
FIm = zeros(460,204,100);
FIm_n = zeros(460,204,100);
Cheat = zeros(460,2);

m = 1;
    EC_g = EC*1;

    bg = [mean(EC_g)+mean(TBC*m)+mean(TuC*m)+mean(SaC*m)+mean(VeC*m)+...
        mean(OrC*m)+mean(ScC*m)+mean(FRC*m)+mean(CyC*m)+mean(NeC*m)];
%bgm_add = bg_measured*0.3+mean(EC)*0.3;
%bgm_add = mean(EC)*0.3;
bgm_add = 0;
for i=1:100


    c=1;
    for k = 1:9
        %bg_p = bg-eval(strcat('mean(',FLs{k,1},')*m'));
        FIm(c:c+9,1:204,i)= ([eval(strcat(FLs{k,1},'(randsample(FLs{k,2},10),:)'))]+bgm_add);
        Cheat(c:c+9,1) = FLs{k,3};
        Cheat(c:c+9,2) = NaN;
        c=c+10;
        if k<9
            for j=(k+1):9
                %bg_p_s = bg_p - eval(strcat('mean(',FLs{j,1},')*m'));
                FIm(c:c+9,1:204,i)= (eval(strcat(FLs{k,1},'(randsample(FLs{k,2},10),:)')) +...
                    eval(strcat(FLs{j,1},'(randsample(FLs{j,2},10),:)')) + bgm_add);
                Cheat(c:c+9,1) = FLs{k,3};
                Cheat(c:c+9,2) = FLs{j,3};
                c=c+10;
            end
        end
    end
    FIm(c:c+9,:,i) =  EC(randsample(259,10),:)*0.3;
    Cheat(c:c+9,1) = 4;
    Cheat(c:c+9,2) = NaN;

         %FIm_noise = awgn(FIm(:,:,1),6,'measured');
    %     FIm_n_noise = zeros(250,204);
    for m = 1:460
        %FIm_n(m,:,i)=FIm_noise(m,:,i)/max(FIm_noise(m,:,i));
        FIm_n(m,:,i)=FIm(m,:,i)/max(FIm(m,:,i));
    end
end

CCT = IdentifyFluorophore_2sC_2Hit(hek, FIm_noise, FIm_n,PD);
PercCorr_c(i,3) = sum(CCT==Cheat)/250;
CCT_c(:,i)=CCT;

mean(PercCorr_c)