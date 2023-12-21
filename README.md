
System Requirements: Matlab Version 2022b
Install Guide: Download all files and add to single folder in Matlab path. Install time is less than 5 minutes. 
Demo Neurocipher:
  1. Open files "Detection_Constants" and "YAS21272_Cells" into workspace. Detection constantss include PD: expected distribution of betas, hek: the reference fluorophore spectra in order. YAS21272_cells are the spectra 
     from each cell identified from animal YAS21272. Cells_line are the raw values and cells_n_line are the cells normalized to their maximum spectral bin.
  2. Type "CCT_comb = IdentifyFluorophore_2StepCorr(hek, cells_line, cells_n_line,PD);" Run time is 5 seconds.
  3. Output variable is the expected fluorophore in spectral order: 1=mTagBFP2, 2=mTurquoise, 3=t-Sapphire, 4=GCaMP, 5=mVenus, 6=mOrange2, 7=mScarlet, 8=FusionRed, 9=CyRFP, 10=mNeptune2.5.
Demo Modeling:
  1. Open files "Detection_Constants" and "Modeling". Detection constantss include PD: expected distribution of betas, hek: the reference fluorophore spectra in order. Modeling contains spctrum from single fluorophore 
     hek cell samples.
  2. Open "ModelingRatios". Change line 3 to alter distributions. Change lines 25-30 to incorporate white noise. Play the script. Once connected to the Parallel Pool, the script shuld take less thana 10 seconds.
  3.  The final value in Perc_corr tells the percent correct from Neurocipher. Cheat is the map for the example dataset. CCT_c is what Neurocipher predicts.
  4.  To find what types of errors Neurocipher produced type: "Error_c = CalcWhereErrorLies(Cheat,CCT_c)". This will run in under one second.
  5.  Error_c has the percentages of each error for ech fluoruophore 1-10. The first column is correct, the second is a false negative, the third is positive for GCaMP, and the fourth an incorrect match.
To run the software on the data and recopitulate results: Cycle through running the animal cell data into the IdentifyFluorophore_2StepCorr code. Run various conditions through the modeling codes and idenntify errors.  
