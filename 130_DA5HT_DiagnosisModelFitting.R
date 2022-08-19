#*********************************************


#Models - skipped this code, it is old and inefficient. Delete?



for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByJustDiagnosis_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByJustDiagnosiss_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByJustDiagnosisRINRNAConc_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByJustDiagnosisRINRNAConc_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByJustDiagnosisRINRNAConcCard_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByJustDiagnosisRINRNAConcCard_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByTheUsualSuspects_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByTheUsualSuspects_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL. + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByTheUsualSuspectsAndRNAconc_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByTheUsualSuspectsAndRNAconc_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#RNAconc seems to actually matter for this dataset. I wonder if I should add it back in to the GabaGlu analyses to make them parallel.


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Wow - card super-duper matters for almost everything.

for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_AllGenesByTheUsualSuspectsRNAConcAndCard_NegDeltaCq.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed3)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_AllGenesByTheUsualSuspectsRNAConcAndCard_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Let's output things better so we can add FDR, etc, plot trends, etc.

i<-1
Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
#Example - to get df for Model:
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 36.1613  1  1.816e-09 ***
# Diagnosis                    5.1395  2    0.07655 .  
# Age                          4.1935  1    0.04058 *  
# Gender                       1.8208  1    0.17722    
# pH                           3.2468  1    0.07156 .  
# Hours.Final                  0.1529  1    0.69575    
# TZP_BioAnalyzer_RIN          3.0480  1    0.08084 .  
# TZP_Average.RNAConc..ng.uL.  0.4282  1    0.51287    
# Card                        71.2767 19  5.635e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR<-MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,i], proc="BH")
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR.csv")

head(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)
head(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)

row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.05]
#[1] "DRD4"  "HTR2B"
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.10]
#[1] "DRD4"  "HTR2B"

row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.20]
#[1] "DRD4"  "HTR2B" "TFRC" 

#What about nominal relationships with diagnosis?
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,2]<0.05]
# [1] "DRD2"  "DRD4"  "HTR2B" "TFRC" 


#Detailed results for these genes - using nlme so that I get p-values for the individual coefficients.
#I looped this later, but it is still included here in order to get the df for the model.

Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)=="HTR2B"], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC    logLik
# 378.3289 469.4375 -159.1645
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept) Residual
# StdDev:   0.7742004 0.424006
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +      TZP_Average.RNAConc..ng.uL. + Card 
# Value Std.Error DF   t-value p-value
# (Intercept)                 -22.319932  5.579326 66 -4.000471  0.0002
# DiagnosisBP                  -0.680805  0.300313 60 -2.266984  0.0270
# DiagnosisSchiz               -0.979576  0.280438 60 -3.493016  0.0009
# Age                          -0.000754  0.008870 60 -0.084972  0.9326
# GenderF                      -0.368348  0.408498 60 -0.901713  0.3708
# pH                            1.915952  0.804924 60  2.380290  0.0205
# Hours.Final                   0.005908  0.017461 60  0.338342  0.7363
# TZP_BioAnalyzer_RIN           0.057568  0.228895 60  0.251505  0.8023
# TZP_Average.RNAConc..ng.uL.   0.000561  0.001942 60  0.289110  0.7735
#I didn't include all of the stats for card.

# Number of Observations: 154
# Number of Groups: 69 

Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)=="DRD4"], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC    logLik
# 302.6371 393.7457 -121.3186
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:   0.3070512 0.4556526
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +      TZP_Average.RNAConc..ng.uL. + Card 
#                                  Value Std.Error DF   t-value p-value
# (Intercept)                 -10.400298 2.9588796 66 -3.514944  0.0008
# DiagnosisBP                  -0.465775 0.1602148 60 -2.907192  0.0051
# DiagnosisSchiz               -0.347007 0.1502869 60 -2.308965  0.0244
# Age                          -0.014961 0.0046783 60 -3.197909  0.0022
# GenderF                      -0.255250 0.2217705 60 -1.150964  0.2543
# pH                           -0.446930 0.4294802 60 -1.040629  0.3022
# Hours.Final                   0.006553 0.0093840 60  0.698295  0.4877
# TZP_BioAnalyzer_RIN           0.547508 0.1224543 60  4.471121  0.0000
# TZP_Average.RNAConc..ng.uL.  -0.003361 0.0010354 60 -3.246236  0.0019

# Number of Observations: 154
# Number of Groups: 69 

#Let's just loop it:
MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[i,]<-summary(Model)$tTable[,1][c(1:9)]
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE[i,]<-summary(Model)$tTable[,2][c(1:9)]
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat[i,]<-summary(Model)$tTable[,4][c(1:9)]
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[i,]<-summary(Model)$tTable[,5][c(1:9)]
  rm(Temp)
  rm(Model)
}


write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas.csv")
write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE.csv")
write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat.csv")
write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval.csv")

#False detection correction:
MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR<-MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,i], proc="BH")
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR.csv")

row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,2]<0.05]
#character(0)
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,3]<0.05]
#[1] "HTR2B"
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,2]<0.10]
#character(0)
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,3]<0.10]
#[1] "HTR2B"

#Interesting - when considered independently (vs. in an ANOVA), only HTR2B has an effect that surpasses false detection correction.

#Genes nominally related to BP:
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,2]<0.05]
#[1] "DRD4"  "HTR2B" "TFRC" 

#Genes nominally related to Schiz:
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,3]<0.05]
#[1] "DBH"   "DRD2"  "DRD3"  "DRD4"  "HTR2B" "MAOB"  

pdf("DA5HT_BipolarVsSchizophreniaBetas.pdf", width=5, height=5)
plot(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2]~MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3], xlab="Coefficients (Betas): Schizophrenia", ylab="Coefficients (Betas): Bipolar Disorder", main="DA-5HT Dataset")
BestFitLine<-lm(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2]~MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3])
abline(BestFitLine)
dev.off()

cor(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2],MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3])
#[1] 0.8492607

pdf("DA5HT_BipolarVsSchizophreniaBetas_NoHK.pdf", width=5, height=5)
plot(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(11,31:46),2]~MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(11,31:46),3], xlab="Log(2) Fold Change: Schizophrenia", ylab="Log(2) Fold Change: Bipolar Disorder", main="DA-5HT Dataset")
BestFitLine<-lm(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(11,31:46),2]~MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(11,31:46),3])
abline(BestFitLine)
dev.off()

cor(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(11,31:46),2],MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(11,31:46),3])
#[1] 0.8787403

#*********************

#Without Card:

i<-1
Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare<-matrix(0, 46, 8)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_ChiSquare)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval<-matrix(0, 46, 8)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)
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

for(i in c(1:8)){
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
#[1] "HTR2B"

#What about nominal relationships with diagnosis?
row.names(MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval)[MLM_DA5HT_AllGenesByTheUsualSuspects_NegDeltaCq_Pval[,2]<0.05]
# [1] "DBH"   "DRD1"  "DRD4"  "HTR2B" "MAOB" 

#****************************


#Just Diagnosis - because I know that someone is going to ask for it:
i<-1
Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)

MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare<-matrix(0, 46, 2)
colnames(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval<-matrix(0, 46, 2)
colnames(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)
  MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare, "MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval, "MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR<-MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval

for(i in c(1:2)){
  tempPvalAdj<-mt.rawp2adjp(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval[,i], proc="BH")
  MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR, "MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR.csv")

head(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR)
head(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval)

row.names(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.05]
#[1] "HTR2B"
row.names(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.10]
#[1] "HTR2B"

row.names(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.20]
#[1] "HTR2B" 

#What about nominal relationships with diagnosis?
row.names(MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval)[MLM_DA5HT_AllGenesByDiagnosis_NegDeltaCq_Pval[,2]<0.05]
#[1] "HTR1B" "HTR1F" "HTR2A" "HTR2B" "HTR5A" "MAOB" 

#*********************
#A version without correcting for HK:

i<-1
Temp<-data.frame(y=(-1*DA5HT_Cq_AllSubjects_QCed3[i,]), SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare)<-row.names(DA5HT_Cq_AllSubjects_QCed3)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval<-matrix(0, 46, 9)
colnames(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)<-row.names(DA5HT_Cq_AllSubjects_QCed3)


for(i in c(1:length(DA5HT_Cq_AllSubjects_QCed3[,1]))){
  Temp<-data.frame(y=(-1*DA5HT_Cq_AllSubjects_QCed3[i,]), SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare.csv")

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR<-MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval[,i], proc="BH")
  MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR, "MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR.csv")

head(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)
head(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)

row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,2]<0.05]
#"HTR2B"
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,2]<0.10]
#"DRD4"  "HTR2B"

#What about nominal relationships with diagnosis?
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval[,2]<0.05]
#[1] "DRD2"  "DRD4"  "HTR2B"

#These are almost exactly the same results identified with -DeltaCq.


#*******************************

#I would like to formally set up a sensitivity analysis (vs. my bootleg output from earlier)

MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels<-matrix(0, 46, 16)
colnames(MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRINConc_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels<-matrix(0, 46, 16)
colnames(MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRINConc_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels<-matrix(0, 46, 16)
colnames(MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRINConc_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)

MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels<-matrix(0, 46, 16)
colnames(MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRINConc_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(DA5HT_NegDeltaCq_AllSubjects_QCed2[1,]))){
  print(i)
  Temp<-data.frame(y=DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
  
  Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,1]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,1]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,1]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,1]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,2]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,2]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,2]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,2]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,3]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,3]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,3]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,3]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,4]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,4]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,4]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,4]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,5]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,5]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,5]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,5]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL., random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,6]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,6]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,6]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,6]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,7]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,7]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,7]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,7]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,8]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,8]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,8]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,8]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  rm(Temp)
  
  print(i)
  Temp<-data.frame(y=-(DA5HT_Cq_AllSubjects_QCed3[i,]), SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, HK=DA5HT_MeanHousekeeping_QCed2)
  
  Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,9]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,9]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,9]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,9]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,10]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,10]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,10]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,10]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,11]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,11]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,11]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,11]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+Card+HK, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,12]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,12]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,12]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,12]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL., random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,13]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,13]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,13]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,13]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,14]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,14]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,14]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,14]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,15]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,15]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,15]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,15]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+HK, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,16]<-summary(Model)$tTable[,1][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,16]<-summary(Model)$tTable[,5][2]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,16]<-summary(Model)$tTable[,1][3]
  MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,16]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  rm(Temp)
  
}

write.csv(cbind(MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels, MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels), "MLM_DA5HT_AllGenes_NegDeltaCq_LM_Bipolar_AcrossModels.csv")

write.csv(cbind(MLM_DA5HT_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels, MLM_DA5HT_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels), "MLM_DA5HT_AllGenes_NegDeltaCq_LM_Schiz_AcrossModels.csv")

#The top results are remarkably stable across models (HTR2B in BP/Schiz, DRD4 in BP)
#There are a few results that are less stable: DBH, TFRC

#************************************************
