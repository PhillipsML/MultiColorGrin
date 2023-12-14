clear Cheat 
    dist = [15;34;35;0;17;49;52;17;7;24];
    Cheat(1:dist(1),1)=1;Cheat(end+1:end+dist(2),1)=2;Cheat(end+1:end+dist(3),1)=3;
    Cheat(end+1:end+dist(4),1)=4;Cheat(end+1:end+dist(5),1)=5;
    Cheat(end+1:end+dist(6),1)=6;Cheat(end+1:end+dist(7),1)=7;
    Cheat(end+1:end+dist(8),1)=8;Cheat(end+1:end+dist(9),1)=9;Cheat(end+1:end+dist(10),1)=10;
mults = 0:0.1:1;
for k = 1:10
    clear CCT
%     m = mults(k);
    EC_g = EC*0.3;
bg = bg_GRIN*0.35+mean(EC_g);
    %bg = [mean(EC_g)+mean(TBC*m)+mean(TuC*m)+mean(SaC*m)+mean(VeC*m)+...
        %mean(OrC*m)+mean(ScC*m)+mean(FRC*m)+mean(CyC*m)+mean(NeC*m)];
    EC_g_b = EC_g + bg;
    % ExG = mean(EC_g);

%     TBC_g_b = TBC+bg-mean(TBC*m);
%     TuC_g_b = TuC+bg-mean(TuC*m);
%     SaC_g_b = SaC+bg-mean(SaC*m);
%     VeC_g_b = VeC+bg-mean(VeC*m);
%     OrC_g_b = OrC+bg-mean(OrC*m);
%     ScC_g_b = ScC+bg-mean(ScC*m);
%     FRC_g_b = FRC+bg-mean(FRC*m);
%     CyC_g_b = CyC+bg-mean(CyC*m);
%     NeC_g_b = NeC+bg-mean(NeC*m);
  TBC_g_b = TBC+bg;
    TuC_g_b = TuC+bg;
    SaC_g_b = SaC+bg;
    VeC_g_b = VeC+bg;
    OrC_g_b = OrC+bg;
    ScC_g_b = ScC+bg;
    FRC_g_b = FRC+bg;
    CyC_g_b = CyC+bg;
    NeC_g_b = NeC+bg;


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
        FIm_g_b = [TBC_g_b(tb,:);TuC_g_b(tu,:);SaC_g_b(sa,:);EC_g_b(eg,:);VeC_g_b(ve,:);...
            OrC_g_b(or,:);ScC_g_b(sc,:);FRC_g_b(fr,:);CyC_g_b(cy,:);NeC_g_b(ne,:)];
        %FIm_n = [TBC_n(tb,:);TuC_n(tu,:);SaC_n(sa,:);EC_n(eg,:);VeC_n(ve,:);...
        %OrC_n(or,:);ScC_n(sc,:);FRC_n(fr,:);CyC_n(cy,:);NeC_n(ne,:)];

        FIm_wn_g_b = awgn(FIm_g_b,6,'measured');
        FIm_n_wn_g_b = zeros(250,204);
        for j = 1:250
            FIm_n_wn_g_b(j,:)=FIm_wn_g_b(j,:)/max(FIm_wn_g_b(j,:));
        end
        CCT = IdentifyFluorophore_2sC_2Hit(hek, FIm_wn_g_b, FIm_n_wn_g_b,PD);
        CCT_noG= CCT(:,1);
        CCT_noG(CCT(:,1)==4,1)=CCT(CCT(:,1)==4,2);
        %CCT = IdentifyFluorophore(hek, FIm_g_b, FIm_n_g_b);
        CCT_ratio_match(:,i,k)=CCT_noG;
        PercCorr_ratio_match(i,k) = sum(CCT_noG==Cheat)/250;

    end

end
mean(PercCorr_ratio_match)