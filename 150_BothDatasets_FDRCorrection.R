
#Random thought: Would it be better to be doing FDR correction across the two datasets? Are we treating these datasets as independent or related?
#e.g.,

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval,MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)

str(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)
# 'data.frame':	139 obs. of  9 variables:
#   $ (Intercept)                : num  0.31173 0.24539 0.00119 0.3062 0.34497 ...
# $ Diagnosis                  : num  0.00166 0.51735 0.98042 0.67369 0.16369 ...
# $ Age                        : num  6.20e-01 2.02e-08 2.98e-03 8.55e-04 8.76e-01 ...
# $ Gender                     : num  0.391 0.118 0.746 0.111 0.553 ...
# $ pH                         : num  0.66 0.368 0.141 0.653 0.889 ...
# $ Hours.Final                : num  0.322 0.089 0.7 0.539 0.634 ...
# $ TZP_BioAnalyzer_RIN        : num  7.66e-01 1.62e-02 8.09e-05 1.67e-02 3.95e-01 ...
# $ TZP_Average.RNAConc..ng.uL.: num  3.34e-01 2.11e-01 8.07e-01 4.08e-01 1.09e-05 ...
# $ Card                       : num  3.61e-55 2.79e-40 1.31e-16 1.27e-17 2.64e-60 ...

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR<-MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR.csv")

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.05]
#[1] "GPHN"  "HTR2B"
row.names(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.10]
#[1] "ABAT"  "SST"   "GPHN"  "MAPK1" "DRD4"  "HTR2B"


#Seems like almost exactly the same results - and probably an easier way to compile/summarize them.

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare,MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare)

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare.csv")

#Making a version excluding housekeeping genes in the FDR calculations:

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK<-MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(83:93, 104, 124:139),]

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(83:93, 104, 124:139),i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(83:93, 104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare[-c(83:93,  104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare_NoHK.csv")

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK[,2]<0.05]
#[1] "ABAT"  "SST"   "GPHN"  "DRD4"  "HTR2B"
#Wow - that actually made a big difference in the results (which I guess makes sense, since the HK should hypothetically not be significant...)

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR_NoHK[,2]<0.10]
#[1]  "ABAT"  "SST"   "GPHN"  "MAPK1" "DRD4"  "HTR2B"

#Although in the end it is basically just the same genes. So that's comforting.

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(83:93, 104, 124:139),])[MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(83:93, 104, 124:139),2]<0.05]
# [1] "ABAT"    "SST"     "BDNF"    "GABRB1"  "GAD1"    "GNAQ"    "GPHN"    "GRIK2"   "GRM5"    "HOMER1"  "GFAP"    "MAPK1"  
# [13] "NSF"     "SLC38A1" "SNCA"    "DRD2"    "DRD4"    "HTR2B" 
length(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(83:93, 104, 124:139),2])
#[1] 111
length(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,2])
#[1] 139


#Let's do it for the LM results too:

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas,MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas)

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE,MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE)

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat,MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat)

MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval,MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)


MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR<-MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval.csv")

#As an aside: The 5 genes included in both datasets now have a 1 after their names for the results from the DA5HT dataset.


#FDR calculations excluding housekeeping genes:
MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR_NoHK<-MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[-c(83:93, 104, 124:139),]

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[-c(83:93,104,124:139),i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR_NoHK, "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104,124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE[-c(83:93, 104,124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat[-c(83:93, 104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[-c(83:93, 104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval_NoHK.csv")


pdf("BothDatasets_BipolarVsSchizophreniaBetas.pdf", width=5, height=5)
plot(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2]~MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3], xlab="Log(2) Fold Change: Schizophrenia", ylab="Log(2) Fold Change: Bipolar Disorder", main="GABA-GLU Dataset")
BestFitLine<-lm(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2]~MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3])
abline(BestFitLine)
dev.off()

cor(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2],MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3])
#[1] 0.7904778

pdf("BothDatasets_BipolarVsSchizophreniaBetas_noHK.pdf", width=5, height=5)
plot(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104, 124:139),2]~MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104, 124:139),3], xlab="Log(2) Fold Change: Schizophrenia", ylab="Log(2) Fold Change: Bipolar Disorder", main="GABA-GLU Dataset")
BestFitLine<-lm(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104, 124:139),2]~MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104, 124:139),3])
abline(BestFitLine)
dev.off()

cor(MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104, 124:139),2],MLM_BothDatasets_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93, 104, 124:139),3])
#[1] 0.8222087



###################################

#One more thing:  Doing the FDR correction across datasets without Card (although I included CardBlock2 for GabaGlu) to avoid overfitting:

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval[,-9],MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)

str(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)
# 'data.frame':	139 obs. of  8 variables:
#   $ (Intercept)                : num  0.288852 0.238289 0.000527 0.175374 0.161874 ...
# $ Diagnosis                  : num  0.00706 0.52316 0.75738 0.93991 0.27913 ...
# $ Age                        : num  2.05e-01 1.19e-09 4.47e-03 2.25e-03 7.97e-01 ...
# $ Gender                     : num  0.9714 0.0213 0.8186 0.8933 0.1763 ...
# $ pH                         : num  0.883 0.437 0.134 0.692 0.659 ...
# $ Hours.Final                : num  0.5022 0.0557 0.8089 0.9794 0.17 ...
# $ TZP_BioAnalyzer_RIN        : num  0.319997 0.347157 0.000307 0.113899 0.549273 ...
# $ TZP_Average.RNAConc..ng.uL.: num  0.90364 0.17147 0.61678 0.47566 0.00706 ...

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR<-MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval

for(i in c(1:8)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[,i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR.csv")

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,2]<0.05]
#[1] character(0)
row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,2]<0.10]
#[1] "HTR2B"


#Seems like almost exactly the same results - and probably an easier way to compile/summarize them.

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval.csv")

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare[,-9],MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare)

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare.csv")

#Making a version excluding housekeeping genes in the FDR calculations:

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK<-MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[-c(83:93, 104, 124:139),]

for(i in c(1:8)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[-c(83:93, 104, 124:139),i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[-c(83:93, 104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare[-c(83:93,  104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare_NoHK.csv")

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,2]<0.05]
#character(0)

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,2]<0.10]
#[1]  "HTR2B"

#The same as before.

#############################


#One more thing:  Doing the FDR correction across datasets for the just Diagnosis analysis:

MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval<-rbind.data.frame(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval,MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval)

str(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval)
# 'data.frame':	139 obs. of  2 variables:
#   $ (Intercept): num  3.61e-319 1.34e-273 3.85e-171 0.00 2.54e-100 ...
# $ Diagnosis  : num  0.0153 0.3208 0.7946 0.9561 0.44 ...

MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR<-MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval

for(i in c(1:2)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval[,i], proc="BH")
  MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR, "MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR.csv")

row.names(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.05]
#[1] "HTR2B"
row.names(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.10]
#[1] "HTR2B"


#Seems like almost exactly the same results - and probably an easier way to compile/summarize them.

write.csv(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval, "MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval.csv")

MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_ChiSquare<-rbind.data.frame(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare,MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare)

write.csv(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_ChiSquare, "MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_ChiSquare.csv")

#Making a version excluding housekeeping genes in the FDR calculations:

row.names(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval)

MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK<-MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval[-c(83:93, 104, 124:139),]

for(i in c(1:2)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval[-c(83:93, 104, 124:139),i], proc="BH")
  MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK, "MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval[-c(83:93, 104, 124:139),], "MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_Pval_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_ChiSquare[-c(83:93,  104, 124:139),], "MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_ChiSquare_NoHK.csv")

row.names(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK[,2]<0.05]
#[1] "HTR2B"

row.names(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK[,2]<0.10]
#[1]  "HTR2B"

row.names(MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByDiagnosis_NegDeltaCq_FDR_NoHK[,2]<0.20]
#[1] "HTR2B"

#The same as before.


###################################

#Again, worrying about overfitting:

#A version without Gender or Card (but with CardBlock2)



#I went back and also outputted the results with CardBlock2 (instead of Card) as a technical co-variate because I was worried about overfitting. It showed a similar pattern of results but appeared to simply not control for technical variation as well.

i<-1
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2 + (1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare<-matrix(0, 93, 8)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval<-matrix(0, 93, 8)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2 + (1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval

for(i in c(1:8)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR.csv")

head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR)
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval)

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR[,2]<0.05]
#character(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR[,2]<0.10]
#character(0)

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_FDR[,2]<0.20]
#character(0)"

#What about nominal relationships with diagnosis?
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval[,2]<0.05]
#[1] "ABAT"  "SST"   "GNAQ"  "GFAP"  "MAPK1" "NSF"   "P2RX7" "SNCA" 


#*********************

#DA5HT Without Card or Gender:

i<-1
Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare<-matrix(0, 46, 7)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval<-matrix(0, 46, 7)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)
  MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare, "MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval, "MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR<-MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval

for(i in c(1:7)){
  tempPvalAdj<-mt.rawp2adjp(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[,i], proc="BH")
  MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR, "MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR.csv")

head(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)
head(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)

row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,2]<0.05]
#[1] "HTR2B"
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,2]<0.10]
#[1] "DRD4"  "HTR2B"

#What about nominal relationships with diagnosis?
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)[MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[,2]<0.05]
# [1] "DBH"   "DRD2"  "DRD4"  "HTR1B" "HTR2B" "MAOB" 


##############

#A version without card or Gender:

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval[,-8],MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)

str(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)
# 'data.frame':	139 obs. of  7 variables:
#   $ (Intercept)                : num  0.288974 0.204638 0.000544 0.172397 0.192028 ...
# $ Diagnosis                  : num  0.00539 0.75519 0.7708 0.92533 0.20828 ...
# $ Age                        : num  2.05e-01 6.58e-09 4.37e-03 2.27e-03 8.27e-01 ...
# $ pH                         : num  0.886 0.317 0.138 0.701 0.773 ...
# $ Hours.Final                : num  0.496 0.113 0.786 0.967 0.23 ...
# $ TZP_BioAnalyzer_RIN        : num  0.304939 0.110815 0.000129 0.091686 0.822305 ...
# $ TZP_Average.RNAConc..ng.uL.: num  0.8968 0.3491 0.5821 0.4836 0.0139 ...

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR<-MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval

for(i in c(1:7)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[,i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR.csv")

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,2]<0.05]
#[1] character(0)
row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR[,2]<0.10]
#[1] "HTR2B"


#Seems like almost exactly the same results - and probably an easier way to compile/summarize them.

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval.csv")

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare<-rbind.data.frame(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare[,-8],MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare)

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare.csv")


#Making a version excluding housekeeping genes in the FDR calculations:

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)

MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK<-MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[-c(83:93, 104, 124:139),]

for(i in c(1:7)){
  tempPvalAdj<-mt.rawp2adjp(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[-c(83:93, 104, 124:139),i], proc="BH")
  MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK, "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[-c(83:93, 104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_Pval_NoHK.csv")

write.csv(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare[-c(83:93,  104, 124:139),], "MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare_NoHK.csv")

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,2]<0.05]
#[1] "HTR2B"

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,2]<0.10]
#[1]  "HTR2B"

row.names(MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK)[MLM_BothDatasets_AllGenesByTheUsualSuspects_NegDeltaCq_FDR_NoHK[,2]<0.20]
#[1] "ABAT"  "SST"   "MAPK1" "DRD4"  "HTR2B"

#The same as before, just noisier.

#########################################
