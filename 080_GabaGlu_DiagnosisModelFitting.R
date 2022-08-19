
#At this point I skipped ahead to run QC for the DA5HT dataset before choosing models, so that I could maintain consistency across the datasets.


#Model time!

#Note: this is some earlier code that I did not re-run. I should probably come back and delete it.

#The most basic model:

for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByJustDiagnosis_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByJustDiagnosiss_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}


#The usual suspects and RIN:
for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final +BioAnalyzer_RIN+ (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspects_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspects_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Basic observation: These co-variates matter a lot for many of these genes.


#Adding in Card:
for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Wow - card super-duper matters for almost everything. But it uses up so much df... ouch.


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+CardBlock2 + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCardBlock_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCardBlock_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Man... supposedly this model fits PC1-3 well, but it is much less powerful than using Card explicitly. Almost all of the biological/technical effects (age, PMI, pH, gender, RIN) in addition to diagnosis effects are better detected when Card is included instead of CardBlock. Interesting.

#Adding in RNA Conc (I'm waffling on this one):
for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+CardBlock2 + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsRNAConcAndCardBlock_NegDeltaCq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsRNAConcAndCardBlock_NegDeltaCq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Average RNA Conc seems to rarely matter anymore. I'm inclined to toss it as a co-variate.


#################

#Checking to see if normalization using housekeeping genes is artificially producing diagnosis effects:

str(GabaGlu_Cq_AllSubjects_QCed4) 

for(i in c(1:length(GabaGlu_Cq_AllSubjects_QCed4[,1]))){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed4[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, HK=GabaGlu_MeanHousekeeping_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+HK+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspects_AndHK_Cq.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed4)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspects_AndHK_Cq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#The only results that this changed was Grik5 & Homer1 (just barely) dropped below the threshold for nominal significance. The relationship with Mapk1 became notably stronger.

for(i in c(1:length(GabaGlu_Cq_AllSubjects_QCed4[,1]))){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed4[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, HK=GabaGlu_MeanHousekeeping_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+CardBlock2+HK+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCardBlock_AndHK_Cq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCardBlock_AndHK_Cq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#The results using that form of normalization are not notably different from the results for -deltaCq.
#SNCA and PGK1 drop (barely) below the threshold for nominal significance


for(i in c(1:length(GabaGlu_Cq_AllSubjects_QCed4[,1]))){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed4[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, HK=GabaGlu_MeanHousekeeping_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Card+HK+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCard_AndHK_Cq.txt")
  stats_output <- c(
    print(colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCard_AndHK_Cq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#The results using that form of normalization are not notably different from the results for -deltaCq.
#GABRA4 (barely) slips up above the threshold for nominal significance. Most of the other effects are slightly stronger.
#Altogether, this method appears to be slightly more powerful (as would be expected), but does not fundamentally change the pattern of results. Given the amount of pushback that we are likely to get from reviewers for using a non-traditional method, it is probably best to just present this as a sensitivity analysis (vs. a main result)

for(i in c(1:length(GabaGlu_Cq_AllSubjects_QCed4[,1]))){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed4[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Card+BioAnalyzer_RIN+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCard_Cq.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed4)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsAndCard_Cq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Hmm... hard to say. Without controlling for HK, it is pretty hard to see effects for any of the variables because the data are so noisy. SST is still significant. PVALB & BDNF gain nominal significance.
#...stopped reviewing the results, because I just realized that I should probably include RNAconc if we aren't including HK:

for(i in c(1:length(GabaGlu_Cq_AllSubjects_QCed4[,1]))){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed4[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_AllGenesByTheUsualSuspectsRNAConcAndCard_Cq.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed4)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_AllGenesByTheUsualSuspectsRNAConcAndCard_Cq.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#SST & TFRC significant, GRIK1, HOMER1 nominally significant, PVALB & GAD1 very nominally significant.
#ABAT & GPHYN are no longer significant
#Interesting. Many of those top results are well-validated (SST, PVALB, GAD1, GRIK1).



#Let's output things better so we can add FDR, etc, plot trends, etc.
#Note: After discovering that RNA Conc really mattered for the DA5HT dataset, I ended up re-running with RNAConc added to try to keep the analyses in the two datasets parallel

i<-1
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR.csv")

head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.05]
#[1] "GPHN"
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,2]<0.10]
#[1] "ABAT"  "SST"   "GPHN"  "MAPK1" 

#What about nominal relationships with diagnosis?
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,2]<0.05]
# [1] "ABAT"    "SST"     "BDNF"    "GABRB1"  "GAD1"    "GNAQ"    "GPHN"    "GRIK2"   "GRM5"    "HOMER1"  "GFAP"    "MAPK1"   "NSF"    
# [14] "SLC38A1" "SNCA"    "PGK1"    "TBP"     "TFRC"  

#I skipped this for now because I looped the code below. Delete? I included one below to get the df for the model:
#Detailed results for these genes - using nlme so that I get p-values for the individual coefficients.
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)=="ABAT"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC   logLik
# -46.02123 86.93483 69.01062
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:  0.09766105 0.1149919
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +      TZP_Average.RNAConc..ng.uL. + Card 
# Value Std.Error DF   t-value p-value
# (Intercept)                  0.9227425 1.1150765 60  0.827515  0.4112
# DiagnosisBP                  0.1568544 0.0544520 60  2.880600  0.0055
# DiagnosisSchiz               0.0837282 0.0500128 60  1.674137  0.0993
# Age                         -0.0007066 0.0017438 60 -0.405217  0.6868
# GenderF                     -0.0531804 0.0758618 60 -0.701016  0.4860
# pH                           0.0555198 0.1541991 60  0.360052  0.7201
# Hours.Final                  0.0027189 0.0033544 60  0.810544  0.4208
# TZP_BioAnalyzer_RIN          0.0105308 0.0433238 60  0.243073  0.8088
# TZP_Average.RNAConc..ng.uL.  0.0002864 0.0003622 60  0.790927  0.4321
# CardCard 10.eds              0.0068211 0.1448911 29  0.047078  0.9628
#...
# Number of Observations: 133
# Number of Groups: 69 

Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)=="SST"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC    logLik
# 152.4616 286.1042 -30.23078
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:   0.4225919 0.1446237
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +     TZP_Average.RNAConc..ng.uL. + Card 
# Value Std.Error DF   t-value p-value
# (Intercept)             -4.376455  3.736823 61 -1.171170  0.2461
# DiagnosisBP             -0.285232  0.182031 61 -1.566936  0.1223
# DiagnosisSchiz          -0.438492  0.167683 61 -2.615001  0.0112
# Age                     -0.019678  0.005793 61 -3.396832  0.0012
# GenderF                 -0.513915  0.256730 61 -2.001771  0.0498
# pH                       0.510316  0.516563 61  0.987906  0.3271
# Hours.Final             -0.003784  0.011355 61 -0.333228  0.7401
# TZP_BioAnalyzer_RIN          0.285781  0.143949 61  1.985296  0.0516
#TZP_Average.RNAConc..ng.uL.  0.000186  0.001077 61  0.173182  0.8631

Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)=="GPHN"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC   logLik
# -37.33917 96.30347 64.66959
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:   0.0639299 0.1368579
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +     TZP_Average.RNAConc..ng.uL. + Card 
# Value Std.Error DF   t-value p-value
# (Intercept)             -3.685792 1.0125739 61 -3.640023  0.0006
# DiagnosisBP              0.161451 0.0496722 61  3.250322  0.0019
# DiagnosisSchiz           0.108084 0.0453908 61  2.381178  0.0204
# Age                     -0.002655 0.0015928 61 -1.666780  0.1007
# GenderF                 -0.101145 0.0686118 61 -1.474160  0.1456
# pH                       0.032367 0.1400517 61  0.231105  0.8180
# Hours.Final              0.004733 0.0030570 61  1.548270  0.1267
# TZP_BioAnalyzer_RIN          0.101317 0.0387125 61  2.617154  0.0112
#TZP_Average.RNAConc..ng.uL.  0.000352 0.0002953 61  1.192566  0.2377


Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)=="MAPK1"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC  logLik
# -82.5862 51.05644 87.2931
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)   Residual
# StdDev:   0.1004412 0.09379397
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +     TZP_Average.RNAConc..ng.uL. + Card 
# Value Std.Error DF    t-value p-value
# (Intercept)             -0.9126354 1.0424181 61 -0.8754984  0.3847
# DiagnosisBP              0.1354882 0.0510107 61  2.6560750  0.0101
# DiagnosisSchiz           0.0870200 0.0467559 61  1.8611543  0.0675
# Age                     -0.0051427 0.0016274 61 -3.1600849  0.0025
# GenderF                  0.0137454 0.0711362 61  0.1932263  0.8474
# pH                       0.2097016 0.1442019 61  1.4542228  0.1510
# Hours.Final              0.0030273 0.0031575 61  0.9587679  0.3415
# TZP_BioAnalyzer_RIN          0.0490035 0.0400250 61  1.2243241  0.2255
#TZP_Average.RNAConc..ng.uL.  0.0003194 0.0003022 61  1.0570040  0.2947

#Honestly, I should probably just loop this:
#I'm going to skip outputting all of the results for card though.

i<-1
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")

str(summary(Model))
summary(Model)$coefficients
summary(Model)$coefficients$fixed[c(1:9)]
names(summary(Model)$coefficients$fixed[c(1:9)])

print(Model)
summary(Model)$tTable
#there we go

summary(Model)$tTable[,1][c(1:9)]#Coefficients
summary(Model)$tTable[,2][c(1:9)]#SE
summary(Model)$tTable[,4][c(1:9)]#Tstat
summary(Model)$tTable[,5][c(1:9)]#Pval
#DF is 60 for all fixed effects that aren't card.

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)<-names(summary(Model)$coefficients$fixed[c(1:9)])
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[i,]<-summary(Model)$tTable[,1][c(1:9)]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE[i,]<-summary(Model)$tTable[,2][c(1:9)]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat[i,]<-summary(Model)$tTable[,4][c(1:9)]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[i,]<-summary(Model)$tTable[,5][c(1:9)]
  rm(Temp)
  rm(Model)
}


write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas.csv")
write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMSE.csv")
write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMTstat.csv")
write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval.csv")

#False detection correction:
MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR.csv")


#Interesting - when considered independently (vs. in an ANOVA), neither diagnosis has any effects that surpass FDR correction:

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,2]<0.10]
#character(0)

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMFDR[,3]<0.10]
#character(0)

#Genes nominally related to BP:
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,2]<0.05]
#[1]  "ABAT"  "GPHN"  "GRM5"  "MAPK1" "NSF"   "PGK1"  "TFRC" 

#Genes nominally related to Schiz:
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMPval[,3]<0.05]
#[1] "SST"    "BDNF"   "GABRB1" "GAD1"   "GNAQ"   "GPHN"   "HOMER1" "GFAP"   "NSF"    "SNCA"   


pdf("GabaGlu_BipolarVsSchizophreniaBetas.pdf", width=5, height=5)
plot(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2]~MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3], xlab="Coefficients (Betas): Schizophrenia", ylab="Coefficients (Betas): Bipolar Disorder", main="GABA-GLU Dataset")
BestFitLine<-lm(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2]~MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3])
abline(BestFitLine)
dev.off()

cor(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,2],MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[,3])
#[1] 0.7029425

pdf("GabaGlu_BipolarVsSchizophreniaBetas_noHK.pdf", width=5, height=5)
plot(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93),2]~MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93),3], xlab="Log(2) Fold Change: Schizophrenia", ylab="Log(2) Fold Change: Bipolar Disorder", main="GABA-GLU Dataset")
BestFitLine<-lm(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93),2]~MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93),3])
abline(BestFitLine)
dev.off()

cor(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93),2],MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_LMBetas[-c(83:93),3])
#[1] 0.7281593


############

#I went back and also outputted the results with CardBlock2 (instead of Card) as a technical co-variate because I was worried about overfitting. It showed a similar pattern of results but appeared to simply not control for technical variation as well.

i<-1
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2 + (1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCardBlock2_NegDeltaCq_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2 + (1 | ID), data = Temp, REML=F)
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

for(i in c(1:9)){
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
#[1] "ABAT"  "SST"   "BDNF"  "GNAQ"  "GFAP"  "MAPK1" "SNCA"  "TFRC" 

##############

#Exploratory analysis: Adding in RNA metrics that seem to matter in the PCA analyses:

i<-1
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_MBNI_Nanodrop_260.230+CardBlock2 + (1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_ChiSquare<-matrix(0, 93, 11)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval<-matrix(0, 93, 11)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_MBNI_Nanodrop_260.230+CardBlock2 + (1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_ChiSquare, "MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval, "MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval

for(i in c(1:11)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR, "MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR.csv")

head(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR)
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval)

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR[,2]<0.05]
#character(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR[,2]<0.10]
#character(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_FDR[,2]<0.20]
#character(0)"

#What about nominal relationships with diagnosis?
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval)[MLM_GabaGlu_AllGenesByTheUsualSuspectsRNAMetricsAndCardBlock2_NegDeltaCq_Pval[,2]<0.05]
#[1] "ABAT"   "SST"    "BDNF"   "GABRB1" "GNAQ"   "GFAP"   "MAPK1"  "P2RX7"  "SNCA"   "TFRC"  
#Basically the same top results as before, the results are just noisier.

#Looking at the output, for the most part the exploratory RNA metrics are not as significantly related to gene expression as the other variables already in the model (with the exception of gender, which never matters...)
#Probably worth not including.

###########################

#Just Diagnosis - because I'm going to assume that someone is going to ask us for it:

i<-1
Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare<-matrix(0, 93, 2)
colnames(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval<-matrix(0, 93, 2)
colnames(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)


for(i in c(1:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis + (1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare, "MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_ChiSquare.csv")

write.csv(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval, "MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR<-MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval

for(i in c(1:2)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR, "MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR.csv")

head(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR)
head(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval)

row.names(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.05]
#character(0)
row.names(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.10]
#character(0)

row.names(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_FDR[,2]<0.20]
#character(0)"

#What about nominal relationships with diagnosis?
row.names(MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval)[MLM_GabaGlu_AllGenesByDiagnosis_NegDeltaCq_Pval[,2]<0.05]
#[1] "ABAT"  "SST"   "GRIK5" "MAPK1" "NSF"



##############

#And then finally also outputted the results for the version using just Cq (intead of delta Cq) after discovering that HK related normalization and co-variates were dramatically changing the results


i<-1
Temp<-data.frame(y=(-1*GabaGlu_Cq_AllSubjects_QCed4[i,]), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare)<-row.names(GabaGlu_Cq_AllSubjects_QCed4)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval<-matrix(0, 93, 9)
colnames(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)<-row.names(car::Anova(Model, type="III"))
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)<-row.names(GabaGlu_Cq_AllSubjects_QCed4)


for(i in c(1:length(GabaGlu_Cq_AllSubjects_QCed4[,1]))){
  Temp<-data.frame(y=(-1*GabaGlu_Cq_AllSubjects_QCed4[i,]), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare[i,]<-car::Anova(Model, type="III")[[1]]
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_ChiSquare.csv")

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval.csv")


#Calculating FDR:
library(multtest)

MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR<-MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval

for(i in c(1:9)){
  tempPvalAdj<-mt.rawp2adjp(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval[,i], proc="BH")
  MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,i]<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
  rm(tempPvalAdj)
}

write.csv(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR, "MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR.csv")

head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)
head(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,2]<0.05]
#character(0)
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,2]<0.10]
#[1] "SST"

row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_FDR[,2]<0.20]
#[1] "SST"  "TFRC"

#What about nominal relationships with diagnosis?
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegCq_Pval[,2]<0.05]
#[1] "SST"    "BDNF"   "GAD1"   "GRIK1"  "HOMER1" "TFRC" 


#************************************************

#I would like to formally set up a sensitivity analysis (vs. my bootleg output from earlier)

MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels<-matrix(0, 93, 20)
colnames(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "DiagnosisRNAMetrics_CardBlock2","TheUsualSuspects_CardBlock2", "TheUsualSuspectsRINConc_CardBlock2", "TheUsualSuspectsRNAMetrics_CardBlock2", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRNAMetrics_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels<-matrix(0, 93, 20)
colnames(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "DiagnosisRNAMetrics_CardBlock2","TheUsualSuspects_CardBlock2", "TheUsualSuspectsRINConc_CardBlock2", "TheUsualSuspectsRNAMetrics_CardBlock2", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRNAMetrics_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels<-matrix(0, 93, 20)
colnames(MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "DiagnosisRNAMetrics_CardBlock2","TheUsualSuspects_CardBlock2", "TheUsualSuspectsRINConc_CardBlock2", "TheUsualSuspectsRNAMetrics_CardBlock2", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRNAMetrics_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels<-matrix(0, 93, 20)
colnames(MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels)<-c("JustDiagnosis", "DiagnosisRNAMetrics", "DiagnosisRNAMetrics_Card", "TheUsualSuspects", "TheUsualSuspects_Card", "TheUsualSuspectsRINConc", "TheUsualSuspectsRNAMetrics", "TheUsualSuspectsRINConc_Card", "DiagnosisRNAMetrics_CardBlock2","TheUsualSuspects_CardBlock2", "TheUsualSuspectsRINConc_CardBlock2", "TheUsualSuspectsRNAMetrics_CardBlock2", "Cq_JustDiagnosis", "Cq_DiagnosisRNAMetrics", "Cq_DiagnosisRNAMetrics_Card", "Cq_DiagnosisRNAMetrics_Card_HK","Cq_TheUsualSuspectsRINConc", "Cq_TheUsualSuspectsRNAMetrics", "Cq_TheUsualSuspectsRINConc_Card", "Cq_TheUsualSuspectsRINConc_Card_HK")
row.names(MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2)

#Note: I used -Cq to make the direction of effect comparable

for(i in c(2:length(GabaGlu_NegDeltaCq_AllSubjects_QCed2[1,]))){
  print(i)
  
  Temp<-data.frame(y=GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
  
  #Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,1]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,1]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,1]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,1]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,2]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,2]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,2]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,2]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+Card, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,3]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,3]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,3]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,3]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,4]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,4]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,4]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,4]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+Card, random=~1|ID, data=Temp, na.action = na.omit,  method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,5]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,5]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,5]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,5]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL., random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,6]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,6]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,6]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,6]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,7]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,7]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,7]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,7]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,8]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,8]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,8]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,8]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+CardBlock2, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,9]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,9]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,9]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,9]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+CardBlock2, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,10]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,10]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,10]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,10]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,11]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,11]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,11]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,11]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+CardBlock2, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,12]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,12]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,12]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,12]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  rm(Temp)
  
  print(i)
  
  Temp<-data.frame(y=-(GabaGlu_Cq_AllSubjects_QCed4[i,]), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, HK=GabaGlu_MeanHousekeeping_QCed2)
  
  Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,13]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,13]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,13]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,13]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,14]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,14]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,14]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,14]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+Card, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,15]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,15]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,15]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,15]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230+Card+HK, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,16]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,16]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,16]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,16]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL., random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,17]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,17]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,17]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,17]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.230, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,18]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,18]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,18]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,18]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,19]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,19]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,19]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,19]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+HK, random=~1|ID, data=Temp, na.action = na.omit, method="ML")
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[i,20]<-summary(Model)$tTable[,1][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels[i,20]<-summary(Model)$tTable[,5][2]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[i,20]<-summary(Model)$tTable[,1][3]
  MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels[i,20]<-summary(Model)$tTable[,5][3]
  rm(Model)
  
  rm(Temp)
  print(i)
}

#Multiple convergence errors - I had to keep restarting the loop:

#1 - first model, fine for others

#Ah - control = lmeControl(opt = 'optim') makes a difference. It was causing problems the first time I did this, I removed it and everything went fine except gene 1 model 1.. so I added it back in for that one gene/model.

write.csv(cbind(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels, MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Bipolar_AcrossModels), "MLM_GabaGlu_AllGenes_NegDeltaCq_LM_Bipolar_AcrossModels.csv")

write.csv(cbind(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels, MLM_GabaGlu_AllGenes_NegDeltaCq_LMPval_Schiz_AcrossModels), "MLM_GabaGlu_AllGenes_NegDeltaCq_LM_Schiz_AcrossModels.csv")

#Schiz:
#Wow - there is a strong negative bias in the results from just Cq, including *all HK genes* showing a negative coefficient. Correcting for RNA conc and RIN helps a little bit.
#...so does this bias reflect something real?  Or a technical artifact?

#Following normalization using -DeltaCq, that bias is gone... but is that real?
#Notably, known effects in SST, PVALB, and GAD1 weaken.

#Top results that do not show the same direction of effect in the normalized vs. unnormalized results: ABAT, GNAQ, GPHN, MAPK1, NSF.
#Top results that maintain direction of effect: SST

#BP:
#This sort of bias is not present in the BP results, despite the fact that the BP coefficients are also overwhelmingly negative for the HK genes.
#Although notably the results showing upregulation are definitely more significant after controlling for HK expression (ABAT, GPHN, MAPK1, NSF)

#Notably, it seems like the correlation between the diagnoses (which we have observed in many datasets) is weaker in the non-normalized results:

colnames(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels)

#"Cq_JustDiagnosis" 
plot(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,13]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,13])
cor(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,13], MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,13])
#[1] 0.7967488
summary.lm(lm(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,13]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,13]))
#Intercept: 0.105271 - note shifted intercept.
#Slope: 0.665812

"Cq_TheUsualSuspectsRINConc_Card"  
plot(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,19]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,19])
cor(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,19], MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,19])
#[1] 0.6977904
summary.lm(lm(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,19]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,19]))
#Intercept: 0.039395 - note shifted intercept is smaller.
#Slope: 0.650218

#Diagnosis
plot(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,1]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,1])
cor(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,1], MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,1])
#[1] 0.7973004

summary.lm(lm(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,1]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,1]))
#Intercept: 0.02763 - note shifted intercept is mostly gone.
#Slope: 0.66477

#"TheUsualSuspectsRINConc_CardBlock2"
plot(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,11]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,11])
cor(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,11], MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,11])
#[1] 0.8044492

summary.lm(lm(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,11]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,11]))
#Intercept: 0.008121 - shift in intercept is now gone.
#Slope: 0.739603

#"TheUsualSuspectsRINConc_Card" 
plot(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,8]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,8])
cor(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,8], MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,8])
#[1] 0.7029422

summary.lm(lm(MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Bipolar_AcrossModels[,8]~MLM_GabaGlu_AllGenes_NegDeltaCq_LMBetas_Schiz_AcrossModels[,8]))
#Intercept: 0.031692 - shift in intercept is smaller
#Slope: 0.653607


#... or maybe not. The results for diagnosis (no co-variates) are just shifted. Co-variates for RNA quality/concentration seem to help.



#******************************************