wn = [10 9 8 7 6 5 4 3 2 1 0];
clear Cheat CCT
dist = [25;25;25;25;25;25;25;25;25;25];
Cheat(1:dist(1),1)=1;Cheat(end+1:end+dist(2),1)=2;Cheat(end+1:end+dist(3),1)=3;
Cheat(end+1:end+dist(4),1)=4;Cheat(end+1:end+dist(5),1)=5;
Cheat(end+1:end+dist(6),1)=6;Cheat(end+1:end+dist(7),1)=7;
Cheat(end+1:end+dist(8),1)=8;Cheat(end+1:end+dist(9),1)=9;Cheat(end+1:end+dist(10),1)=10;
for k = 1:11
    parfor i=1:100
        tb = randsample(166,dist(1));
        tu = randsample(258,dist(2));
        sa = randsample(361,dist(3));
        eg = randsample(259,dist(4));
        ve = randsample(272,dist(5));
        or = randsample(424,dist(6));
        sc = randsample(670,dist(7));
        fr = randsample(701,dist(8));
        cy = randsample(442,dist(9));
        ne = randsample(387,dist(10));
        FIm = [TBC(tb,:);TuC(tu,:);SaC(sa,:);EC(eg,:);VeC(ve,:);OrC(or,:);...
            ScC(sc,:);FRC(fr,:);CyC(cy,:);NeC(ne,:)];
        FIm_n = [TBC_n(tb,:);TuC_n(tu,:);SaC_n(sa,:);EC_n(eg,:);VeC_n(ve,:);...
            OrC_n(or,:);ScC_n(sc,:);FRC_n(fr,:);CyC_n(cy,:);NeC_n(ne,:)];

        FIm_noise = awgn(FIm,wn(k),'measured');
        FIm_n_noise = zeros(250,204);
        for j = 1:250
            FIm_n_noise(j,:)=FIm_noise(j,:)/max(FIm_noise(j,:));
        end
        CCT = IdentifyFluorophore_2StepCorr(hek, FIm_noise, FIm_n_noise,PD);
        PercCorr_sn(i,k) = sum(CCT==Cheat)/250;
        CCT_sn(:,i)=CCT;
    end
end
    mean(PercCorr_sn)