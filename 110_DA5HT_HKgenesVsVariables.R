#**************************************************************

#Before normalizing by housekeeping gene:
#Double-checking whether there is a relationship between housekeeping genes and diagnosis: 


#Before doing this, I did not remove any samples, even though some of them had higher Cq, since they didn't seem super extreme, and the replicates were in agreement.
dim(SubjectInfo_OrderedForDA5HTCqMatrix_QCed)
#[1] 154  47
SubjectInfo_OrderedForDA5HTCqMatrix_QCed2<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed

dim(DA5HT_Cq_AllSubjects_QCed2)
#[1]  46 154
DA5HT_Cq_AllSubjects_QCed3<-DA5HT_Cq_AllSubjects_QCed2


#Note - I skipped this for now, because I want to expand other sections first.

#Running analyses requires an MLM because there are 1-2 samples per ID.  

library(lme4)

#The usual suspects:

for(i in c(33:48)){
  Temp<-data.frame(y=DA5HT_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL. + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndRNAConc.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndRNAConc.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Many of the housekeeping genes show a realtionship with age, RIN, and RNAconc.
#Some of the housekeeping genes show relationships with pH and PMI
#None of the housekeeping genes are significantly related to diagnosis.
#Notably, this dataset also contains TFRC, GUSB, and IPO8, but they do not show trends towards a relationship with diagnosis in this dataset like they did in the GabaGlu dataset.
#I went back and re-ran this after making sure ID was defined as character instead of numeric and the output looks exactly the same as before. 

for(i in c(33:48)){
  Temp<-data.frame(y=DA5HT_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Dissecton.Group + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndDissection.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndDissection.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Dissection group seems to matter for a lot of the housekeeping genes. HPRT1 and TFRC and RAB7A show trends towards a relationship with diagnosis. 

for(i in c(33:48)){
  Temp<-data.frame(y=DA5HT_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndCard.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndCard.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Similar to the GABA-GLU dataset, card shows a strong relationship for every gene. 
#After controlling for card, there are no relationships with diagnosis.

for(i in c(33:48)){
  Temp<-data.frame(y=DA5HT_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Dissecton.Group+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndDissectionCard.txt")
  stats_output <- c(
    print(row.names(DA5HT_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="DA5HT_MLM_HousekeepingGenesByTheUsualSuspectsAndDissectionCard.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Dissection still matters even after controlling for card - although they are so collinear I'm not sure I completely believe it.
#No warnings, but this is a huge amount of df to include in a model for this little sample size (20 cards + 12 dissection groups + 6 other variables = 38!)


#Before controlling for housekeeping gene expression let's see if there is any reason to impute gene expression values:

sum(is.na(DA5HT_Cq_AllSubjects_QCed2[c(33:48),]))
#[1] 0
#Nope, looks good.



###Note - I rearranged the code here so this is going to change!!!

# DA5HT_Cq_AllSubjects_QCed3<-DA5HT_Cq_AllSubjects_QCed2[DA5HT_NumberNA_PerGene<20,]
# dim(DA5HT_Cq_AllSubjects_QCed3)

#So I'll skip imputation for this card before calculating average housekeeping gene Cq...
DA5HT_Cq_AllSubjects_QCed3_Imputed<-DA5HT_Cq_AllSubjects_QCed3

DA5HT_TrimmedMeanHousekeeping_QCed<-apply(DA5HT_Cq_AllSubjects_QCed3_Imputed[c(11,31:46),], 2, function(y) mean(y, trim=0.20, na.rm=T))
hist(DA5HT_TrimmedMeanHousekeeping_QCed, breaks=20)

DA5HT_TrimmedMeanNonHousekeeping_QCed<-apply(DA5HT_Cq_AllSubjects_QCed3_Imputed[-c(11,31:46),], 2, function(y) mean(y, trim=0.20, na.rm=T))
hist(DA5HT_TrimmedMeanNonHousekeeping_QCed, breaks=20)

DA5HT_MeanHousekeeping_QCed<-apply(DA5HT_Cq_AllSubjects_QCed3_Imputed[c(11,31:46),], 2, function(y) mean(y, na.rm=T))
hist(DA5HT_MeanHousekeeping_QCed, breaks=20)

DA5HT_MeanNonHousekeeping_QCed<-apply(DA5HT_Cq_AllSubjects_QCed3_Imputed[-c(11,31:46),], 2, function(y) mean(y, na.rm=T))
hist(DA5HT_MeanNonHousekeeping_QCed, breaks=20)


pdf("Plot_DA5HT_MeanNonHK_vs_MeanHK.pdf", width=5, height=5)
plot(DA5HT_MeanNonHousekeeping_QCed~DA5HT_MeanHousekeeping_QCed)
abline(a=0,b=1, col="red")
dev.off()
#weird - it is a little zig-zaggy...
#Also, housekeeping gene expression is much higher (lower Cq) than DA5HT gene expression.

cor(DA5HT_MeanNonHousekeeping_QCed,DA5HT_MeanHousekeeping_QCed)
#[1] 0.9185234
plot(DA5HT_TrimmedMeanNonHousekeeping_QCed~DA5HT_TrimmedMeanHousekeeping_QCed)
cor(DA5HT_TrimmedMeanNonHousekeeping_QCed, DA5HT_TrimmedMeanHousekeeping_QCed)
#[1] 0.9043258
#Very pretty

#Sanity check:
hist(DA5HT_MeanHousekeeping_QCed)
#Centered around 22.25
hist(DA5HT_MeanNonHousekeeping_QCed)
#Centered around 27.75
#So typically housekeeping genes are more highly expressed than non-housekeeping genes.


#Alright, back to looking at the relationship with diagnosis and the usual suspects:

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)


#Let's take a peek at the correlation between housekeeping gene expression and other variables:

Temp<-cor(cbind(DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[, c(7,9,11,16,24:28,38:42)]))
#Similar to PC1, looks like the strongest correlation is with RNAconc (R=-0.68) followed by RIN (-0.49) and 260/280 (0.33)

write.csv(Temp, "DA5HT_CorMatrix_HK_vs_NumericVar_Cq_QCed2.csv")


#To make plots more intuitive, I inverted the Cq measurements (so that increase = greater expression)

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_RIN_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_RIN, xlab="TZP RNA Integrity Number (RIN)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_RIN)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#HK expression increases with RIN

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_RNAConc_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL., xlab="TZP RNA Concentration (ng/uL)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#HK expression increases with RNA conc - this is interesting because the same *quantity* of RNA was included for each sample in the reaction.

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_260280_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.280, xlab="TZP RNA Purity (260/280)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.280)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#negative relationship, maybe useful for correcting RNA conc?

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_260230_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.230, xlab="TZP RNA Purity (260/230)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.230)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#slight positive relationship, but maybe useful for correcting RNA conc?

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_28s18s_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s., xlab="TZP RNA Integrity (28s/18s rRNA)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#Very slight positive relationship.


pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_PMI_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Hours.Final, xlab="Post-Mortem Interval (PMI)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Hours.Final)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#A slight decrease with higher PMI - probably mediated by RIN

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_pH_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$pH, xlab="Brain pH", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$pH)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#No relationship with pH

pdf("DA5HT_Scatterplot_MeanHousekeepingCq_vs_Age_QCed2.pdf", width=6, height=6)
plot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Age, xlab="Age (yrs)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Age)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#A slight increase in HK with age.

pdf("DA5HT_Boxplot_MeanHousekeepingCq_vs_DissectionGroup_QCed2.pdf", width=12, height=6)
boxplot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Dissecton.Group, xlab="Dissection/Extraction Group", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Dissecton.Group, xlab="Dissection/Extraction Group", ylab="Mean Housekeeping Gene Expression (-Cq)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


pdf("Boxplot_MeanHousekeepingCq_vs_Card_QCed2.pdf", width=25, height=6)
par(mar=c(10,10,5,5))
boxplot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey", las=2)
stripchart(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black', las=2)
dev.off()
#Card is less of a source of variation in this dataset.

pdf("DA5HT_Boxplot_MeanHousekeepingCq_vs_Diagnosis_QCed2.pdf", width=4, height=6)
boxplot(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Diagnosis, xlab="Diagnosis", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(-DA5HT_MeanHousekeeping_QCed~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

#So Bipolar and Schizophrenia are associated with slightly lower housekeeping gene expression (higher Cq)


#Alright, back to looking at the relationship with diagnosis and the usual suspects:

#Just diagnosis:

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 1.0968e+05  1     <2e-16 ***
# Diagnosis   2.9453e+00  2     0.2293    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC       BIC   logLik
# -100.1367 -84.95195 55.06836
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)   Residual
# StdDev:   0.3385596 0.06984557
# 
# Fixed effects: y ~ Diagnosis 
# Value  Std.Error DF  t-value p-value
# (Intercept)    22.204577 0.06770822 85 327.9451  0.0000
# DiagnosisBP    -0.034293 0.10131056 66  -0.3385  0.7361
# DiagnosisSchiz  0.133246 0.10002960 66   1.3321  0.1874


#Diagnosis and the usual suspects
Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 153.8883  1  < 2.2e-16 ***
# Diagnosis     4.6035  2   0.100086    
# Age           1.9952  1   0.157802    
# Gender        0.0011  1   0.973833    
# pH            0.2998  1   0.584003    
# Hours.Final   6.7966  1   0.009133 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Diagnosis and the usual suspects+dissection/extraction
Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     182.0287  1    < 2e-16 ***
# Diagnosis         4.6026  2    0.10013    
# Age               3.2988  1    0.06933 .  
# Gender            0.0581  1    0.80952    
# pH                0.0981  1    0.75416    
# Hours.Final       6.0093  1    0.01423 *  
# Dissecton.Group  25.9613 12    0.01087 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Diagnosis and RNA technical variables:

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                          88.2778  1  < 2.2e-16 ***
# Diagnosis                             0.4982  2   0.779515    
# Block.Weight..g.                      0.0694  1   0.792261    
# TZP_Average.RNAConc..ng.uL.          63.5977  1  1.526e-15 ***
# TZP_BioAnalyzer_RIN                  55.7298  1  8.315e-14 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  9.7459  1   0.001797 ** 
# TZP_Average.260.280                   0.0658  1   0.797524    
# TZP_Average.260.230                   0.0702  1   0.791095    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Mostly RNAConc & RIN (maybe 28s/18s)

#How much of the variation is just RNAconc and RIN?
#Note that this quick analysis doesn't take into account the replicate structure of the dataset
Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
summary.lm(lm(y~TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN, data = Temp))
# Residual standard error: 0.2127 on 151 degrees of freedom
# Multiple R-squared:  0.6315,	Adjusted R-squared:  0.6266 
# F-statistic: 129.4 on 2 and 151 DF,  p-value: < 2.2e-16

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
summary.lm(lm(y~TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data = Temp))
# Residual standard error: 0.1974 on 150 degrees of freedom
# Multiple R-squared:  0.6847,	Adjusted R-squared:  0.6783 
# F-statistic: 108.6 on 3 and 150 DF,  p-value: < 2.2e-16

#28s/18s does actually add a good 5% to the R-squared for the model. Interesting. I had assumed they would be redundant.

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN, data = Temp))
# Residual standard error: 0.2099 on 145 degrees of freedom
# Multiple R-squared:  0.6553,	Adjusted R-squared:  0.6363 
# F-statistic: 34.46 on 8 and 145 DF,  p-value: < 2.2e-16

#Biological variables don't add much (2%) - a little bit though.

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data = Temp))
# Residual standard error: 0.1969 on 144 degrees of freedom
# Multiple R-squared:  0.6987,	Adjusted R-squared:  0.6799 
# F-statistic: 37.11 on 9 and 144 DF,  p-value: < 2.2e-16 

#The biological variables add even less if 28s/18s is in the model.

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+Card, data = Temp))
# Residual standard error: 0.2098 on 126 degrees of freedom
# Multiple R-squared:  0.7007,	Adjusted R-squared:  0.6366 
# F-statistic: 10.92 on 27 and 126 DF,  p-value: < 2.2e-16

#Notably, the jump in R-Squared after adding Card is only 5%, despite using up a large number of df
#And there is no increase in Adjusted R-squared.

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Card, data = Temp))
# Residual standard error: 0.1944 on 125 degrees of freedom
# Multiple R-squared:  0.7452,	Adjusted R-squared:  0.6882 
# F-statistic: 13.06 on 28 and 125 DF,  p-value: < 2.2e-16

#Same verdict for the model with 28s/18s. Small gain in adj.R-squared

#Technical variables and card:
Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           69.6880  1  < 2.2e-16 ***
# Diagnosis                              0.5453  2     0.7614    
# Block.Weight..g.                       0.0589  1     0.8083    
# TZP_Average.RNAConc..ng.uL.           50.7337  1  1.058e-12 ***
# TZP_BioAnalyzer_RIN                   34.6910  1  3.864e-09 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   2.0373  1     0.1535    
# TZP_Average.260.280                    0.1453  1     0.7031    
# TZP_Average.260.230                    0.2120  1     0.6452    
# Card                                 230.6905 19  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Most important technical variables + the usual suspects
Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 485.1953  1  < 2.2e-16 ***
# Diagnosis                     0.8321  2    0.65966    
# Age                           4.1011  1    0.04286 *  
# Gender                        0.0459  1    0.83029    
# pH                            0.0228  1    0.88002    
# Hours.Final                   0.9431  1    0.33149    
# TZP_Average.RNAConc..ng.uL.  70.0591  1  < 2.2e-16 ***
# TZP_BioAnalyzer_RIN          29.9698  1  4.388e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                          562.9830  1  < 2.2e-16 ***
# Diagnosis                              1.1703  2  0.5570158    
# Age                                    2.1648  1  0.1412044    
# Gender                                 0.3411  1  0.5592149    
# pH                                     0.9775  1  0.3228282    
# Hours.Final                            0.5586  1  0.4548244    
# TZP_Average.RNAConc..ng.uL.           89.9396  1  < 2.2e-16 ***
# TZP_BioAnalyzer_RIN                   45.7645  1  1.334e-11 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  11.0158  1  0.0009034 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Most important technical variables + the usual suspects+Card
Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 397.1573  1  < 2.2e-16 ***
# Diagnosis                     1.0600  2    0.58859    
# Age                           4.4269  1    0.03538 *  
# Gender                        0.0493  1    0.82435    
# pH                            0.0025  1    0.95992    
# Hours.Final                   0.2698  1    0.60348    
# TZP_Average.RNAConc..ng.uL.  61.8920  1  3.628e-15 ***
# TZP_BioAnalyzer_RIN          27.0860  1  1.946e-07 ***
# Card                        240.3334 19  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                          400.7031  1  < 2.2e-16 ***
# Diagnosis                              1.1664  2    0.55811    
# Age                                    3.1308  1    0.07683 .  
# Gender                                 0.1359  1    0.71234    
# pH                                     0.1704  1    0.67976    
# Hours.Final                            0.1583  1    0.69070    
# TZP_Average.RNAConc..ng.uL.           66.6731  1  3.205e-16 ***
# TZP_BioAnalyzer_RIN                   28.8931  1  7.649e-08 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   2.1675  1    0.14096    
# Card                                 230.7033 19  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC       BIC   logLik
# -231.9141 -140.8056 145.9571
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)   Residual
# StdDev:   0.2115941 0.03521094
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_Average.RNAConc..ng.uL. +      TZP_BioAnalyzer_RIN + Card 
# Value Std.Error DF   t-value p-value
# (Intercept)                 25.733699 1.4275652 66 18.026286  0.0000
# DiagnosisBP                 -0.036777 0.0762935 60 -0.482044  0.6315
# DiagnosisSchiz               0.038432 0.0711543 60  0.540128  0.5911
# Age                         -0.004316 0.0022676 60 -1.903159  0.0618
# GenderF                      0.020568 0.1024492 60  0.200764  0.8416
# pH                          -0.009314 0.2049097 60 -0.045453  0.9639
# Hours.Final                  0.002041 0.0043439 60  0.469817  0.6402
# TZP_Average.RNAConc..ng.uL. -0.003446 0.0004842 60 -7.116103  0.0000
# TZP_BioAnalyzer_RIN         -0.272424 0.0578693 60 -4.707574  0.0000

#I didn't copy and paste all of the card values..
#Also, note that the betas are still inverted (Cq).


Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#Again, the other technical variables don't add much and just use up df
#260/280 and 260/230  should hypothetically correct for bias in RNAconc, but don't seem to be doing much.
#I guess that makes sense, since 260/280 and 260/230 are strongly multicollinear


Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
#TZP_Average.RNAConc..ng.uL._Uncontaminated2<-lm(Temp$TZP_Average.RNAConc..ng.uL.~Temp$TZP_Average.260.280+Temp$TZP_Average.260.230)$residuals+mean(Temp$TZP_Average.RNAConc..ng.uL.)
#TZP_Average.RNAConc..ng.uL._Uncontaminated2<-lm(Temp$TZP_Average.RNAConc..ng.uL.~Temp$TZP_Average.260.280)$residuals+mean(Temp$TZP_Average.RNAConc..ng.uL.)
TZP_Average.RNAConc..ng.uL._Uncontaminated2<-Temp$TZP_Average.RNAConc..ng.uL./Temp$TZP_Average.260.280
Temp2<-data.frame(Temp,TZP_Average.RNAConc..ng.uL._Uncontaminated2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL._Uncontaminated2+TZP_BioAnalyzer_RIN+Card + (1 | ID), data = Temp2, REML=F)
car::Anova(Model, type="III")
#Not as good of a predictor.

#Hmmm... another possibility: 
#I believe the protocol is actually set up to use exactly 1 ug of RNA
#Is it possible that the issue isn't the concentration, but the miscalculation of concentration (based on contamination) which was then used to calculate the 1 ug for the protocol?

Temp<-data.frame(y=DA5HT_MeanHousekeeping_QCed, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
#TZP_Average.RNAConc..ng.uL._Uncontaminated2<-lm(Temp$TZP_Average.RNAConc..ng.uL.~Temp$TZP_Average.260.280+Temp$TZP_Average.260.230)$residuals+mean(Temp$TZP_Average.RNAConc..ng.uL.)
TZP_Average.RNAConc..ng.uL._Uncontaminated2<-Temp$TZP_Average.RNAConc..ng.uL./Temp$TZP_Average.260.280
Temp2<-data.frame(Temp,TZP_Average.RNAConc..ng.uL._Uncontaminated2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_Average.RNAConc..ng.uL._Uncontaminated2+TZP_BioAnalyzer_RIN+Card + (1 | ID), data = Temp2, REML=F)
car::Anova(Model, type="III")

#Come back to this - I'm not even sure that is interpretable.

#Take away:
#1. Diagnosis is never related to HK expression in this dataset.
#2. The most important co-variates seem to be Age, RNAConc, and RIN
#3. Card is highly significant in all models, but does not necessarily improve adj.R-squared much (due to using up so much df)
#4. TZP_BioAnalyzer_rRNA.Ratio..28s.18s. may be worth including models moving forward - it does seem to add additional information (not just redundant with RIN)



#******************************

#Out of curiousity, five of the housekeeping genes are found in both datasets (HPRT1, GUSB, B2M, IPO8, TFRC) - how well do their values correlate?
#Note: my first round of analysis for this missed gene IPO8. Later when I came back to re-do it, I didn't re-do everything (because quite a few of the analyses turned out to be less than insightful...) so some of the analyses below still don't include it.

row.names(GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),])
#[1] "HPRT1" "GUSB"  "B2M"   "IPO8"  "TFRC" 

dim(GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),])
#[1]   5 133
dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
#[1] 133  50

row.names(DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),])
#[1] "HPRT1" "GUSB"  "B2M"   "IPO8"  "TFRC" 
dim(DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),])
#[1] 5 154
dim(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
#[1] 154  47

#Let's try a simple version of correlation vs. a more complex model controlling for sample-level technical co-variates.

#Calculating average housekeeping gene expression by subject for shared genes
GabaGlu_SharedHousekeepingGenes_AverageBySubject<-matrix(0,5, 69)
row.names(GabaGlu_SharedHousekeepingGenes_AverageBySubject)<-row.names(GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),])
colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID))

for(i in c(1:5)){
  Temp<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),]
  GabaGlu_SharedHousekeepingGenes_AverageBySubject[i,]<-tapply(Temp[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$ID, function(y) mean(y, na.rm=TRUE))
  rm(Temp)
}
head(GabaGlu_SharedHousekeepingGenes_AverageBySubject)
#         1      10      13      14      15     16      17      18     19       2      20     21      22      23      24      25
# HPRT1 21.643 21.4275 21.4685 21.1390 20.9800 21.723 21.2975 22.8940 21.396 21.4505 21.5400 20.939 21.4710 21.4010 20.9385 21.5385
# GUSB  25.015 24.7520 24.6240 24.9145 24.2265 25.154 25.1315 24.9825 24.638 24.5710 24.5830 24.532 25.2375 25.0955 24.4225 23.8840
# B2M   21.165 21.5055 20.6730 21.2545 20.8705 22.023 21.2005 21.9200 21.393 21.3100 21.5265 21.019 21.4600 22.2555 21.1320 21.6025
# IPO8  24.359 24.0825 24.1540 24.0245 23.8485 24.469 24.1160 24.5280 24.353 24.0645 24.2120 23.804 24.2785 24.3805 23.6265 23.6430
# TFRC  23.955 23.4500 23.6490 23.1540 22.8540 23.415 23.7835 24.0900 23.597 23.2715 23.4705 22.956 23.8035 24.2325 23.3425 23.2205

DA5HT_SharedHousekeepingGenes_AverageBySubject<-matrix(0,5, 69)
row.names(DA5HT_SharedHousekeepingGenes_AverageBySubject)<-row.names(DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),])
colnames(DA5HT_SharedHousekeepingGenes_AverageBySubject)<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$ID))
#Interesting - that isn't the same order as created by tapply in the other dataset.

for(i in c(1:5)){
  Temp<-DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),]
  DA5HT_SharedHousekeepingGenes_AverageBySubject[i,]<-tapply(Temp[i,], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$ID, function(y) mean(y, na.rm=TRUE))
  rm(Temp)
}
head(DA5HT_SharedHousekeepingGenes_AverageBySubject)
#           1      10      13      14      15      16      17       18      19       2      20      21       22      23      24
# HPRT1 21.25100 21.4035 21.3125 20.8055 21.0575 21.3560 21.2545 22.63000 21.1615 21.3550 21.5555 20.5265 21.47575 21.3035 20.9115
# GUSB  24.06275 24.4280 24.0175 24.3125 24.1655 24.3690 24.6005 24.51275 24.2000 24.4695 24.6620 24.0000 25.11325 25.0170 23.8315
# B2M   20.89225 21.3570 20.3530 21.0860 20.7795 21.5470 21.0730 21.50925 21.0620 21.2385 21.4575 20.6905 21.60775 22.0765 20.9735
# IPO8  23.45525 23.7735 23.7525 23.6925 23.8225 23.9245 23.8890 23.72400 23.5780 23.6560 23.7430 23.1390 23.96225 23.8010 23.0290
# TFRC  23.01475 22.6375 22.8470 22.4215 22.5355 22.4680 23.1400 23.34975 22.7650 23.0635 22.9560 22.5350 23.41850 23.8430 22.8975

DA5HT_SharedHousekeepingGenes_AverageBySubject<-DA5HT_SharedHousekeepingGenes_AverageBySubject[,order(colnames(DA5HT_SharedHousekeepingGenes_AverageBySubject))]

cbind(colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject), colnames(DA5HT_SharedHousekeepingGenes_AverageBySubject))
#They're both in the same order. :)

SubjectInfo_forHKplots<-SubjectInfo[SubjectInfo$ID%in%colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject),]
dim(SubjectInfo_forHKplots)
#[1] 69 44

cbind(SubjectInfo_forHKplots$ID[order(as.character(SubjectInfo_forHKplots$ID))], colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject))
#They are both the same order.

SubjectInfo_forHKplots<-SubjectInfo_forHKplots[order(as.character(SubjectInfo_forHKplots$ID)),]
SubjectInfo_forHKplots$ID

#Let's see how they correlate:

#HPRT1:
pdf("HPRT1_AverageCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject[1,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[1,], xlab="DA-5HT Dataset: Average Cq by Subject", ylab="GABA-GLU Dataset: Average Cq by Subject", main="HPRT1", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject[1,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[1,])
abline(BestFitLine)
dev.off()
#Some correlation, mostly driven by a handful of extreme samples
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject[1,],DA5HT_SharedHousekeepingGenes_AverageBySubject[1,])
#[1] 0.6911992

#GUSB:
pdf("GUSB_AverageCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject[2,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[2,], xlab="DA-5HT Dataset: Average Cq by Subject", ylab="GABA-GLU Dataset: Average Cq by Subject", main="GUSB", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject[2,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[2,])
abline(BestFitLine)
dev.off()
#less driven by just a handful of extreme samples:
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject[2,],DA5HT_SharedHousekeepingGenes_AverageBySubject[2,])
#[1] 0.618492

#B2M
pdf("B2M_AverageCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject[3,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[3,], xlab="DA-5HT Dataset: Average Cq by Subject", ylab="GABA-GLU Dataset: Average Cq by Subject", main="B2M", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject[3,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[3,])
abline(BestFitLine)
dev.off()
#Very pretty
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject[3,],DA5HT_SharedHousekeepingGenes_AverageBySubject[3,])
#[1] 0.8873383


pdf("IPO8_AverageCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject[4,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[4,], xlab="DA-5HT Dataset: Average Cq by Subject", ylab="GABA-GLU Dataset: Average Cq by Subject", main="IPO8", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject[4,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[4,])
abline(BestFitLine)
dev.off()
#Very pretty. There is an outlier data points in the GabaGlu dataset
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject[4,],DA5HT_SharedHousekeepingGenes_AverageBySubject[4,])
#[1] 0.6579012


#TFRC
pdf("TFRC_AverageCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject[5,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[5,], xlab="DA-5HT Dataset: Average Cq by Subject", ylab="GABA-GLU Dataset: Average Cq by Subject", main="TFRC", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject[5,]~DA5HT_SharedHousekeepingGenes_AverageBySubject[5,])
abline(BestFitLine)
dev.off()
#Very pretty
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject[5,],DA5HT_SharedHousekeepingGenes_AverageBySubject[5,])
#[1] 0.8366621


#I did not re-run all of the following code because it seemed like a pretty ill-thought out group of analyses. Delete?
#I didn't go back and re-do this when I re-ran the analysis after obtaining RIN values - it seemed like a particularly ill-thought out set of analyses:
#Alright, now let's try two manners of controlling for sample-level technical variation and see if it improves the correlation:
# 
# GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject<-matrix(0, 4, 72)
# row.names(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject)<-row.names(GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 92),])
# colnames(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID))
# 
# for(i in c(1:4)){
#   Temp<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 92),]
#   Temp2<-data.frame(y=Temp[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#   Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Card +(1 | ID), data = Temp2, REML=F)
#   
#   Temp_Coefficients<-data.frame(CardRenamed=names(fixef(Model)), fixef(Model))
#   
#   SubjectInfo_OrderedForGabaGluCqMatrix_QCed2_CardCoefficients<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[,1:2] 
#   
#   SubjectInfo_OrderedForGabaGluCqMatrix_QCed2_CardCoefficients$CardRenamed<-paste("Card", SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, sep="")
#   
#   SubjectInfo_OrderedForGabaGluCqMatrix_QCed2_CardCoefficients<-join(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2_CardCoefficients, Temp_Coefficients, by="CardRenamed", type="left")
#   
#   Temp_CorrectedForCard<-Temp[i,]-SubjectInfo_OrderedForGabaGluCqMatrix_QCed2_CardCoefficients$fixef.Model.
#   
#   GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[i,]<-tapply(Temp_CorrectedForCard, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID, function(y) mean(y, na.rm=TRUE))
#   rm(Temp)
#   rm(Temp2)
#   rm(Model)
#   rm(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2_CardCoefficients)
#   rm(Temp_CorrectedForCard)
# }
# head(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject)
# #       1       10       12       13       14       15       16       17       18       19
# # HPRT1 21.76475 21.37877 22.40408 21.94302 21.61352 21.45452 21.59302 21.11401 22.71051 21.21251
# # GUSB  23.92200 23.78768 24.24299 23.82306 24.11356 23.42556 23.76206 24.06507 23.91607 23.57157
# # B2M   20.91600 20.96687 22.38815 20.69226 21.27376 20.88976 21.36476 20.70289 21.42239 20.89539
# # TFRC  23.48900 23.05283 23.87463 23.77885 23.28385 22.98385 22.92485 23.27355 23.58005 23.08705
# 
# 
# #Previous Sanity check:
# # summary.lm(lm(Temp[i,]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))
# # summary.lm(lm(Temp_CorrectedForCard~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))
# 
# 
# 
# DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject<-matrix(0, 4, 72)
# row.names(DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject)<-row.names(DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:33,35),])
# colnames(DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject)<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$ID))
# 
# for(i in c(1:4)){
#   Temp<-DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:33,35),]
#   Temp2<-data.frame(y=Temp[i,], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
#   Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Card +(1 | ID), data = Temp2, REML=F)
#   
#   Temp_Coefficients<-data.frame(CardRenamed=names(fixef(Model)), fixef(Model))
#   
#   SubjectInfo_OrderedForDA5HTCqMatrix_QCed2_CardCoefficients<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[,1:2] 
#   
#   SubjectInfo_OrderedForDA5HTCqMatrix_QCed2_CardCoefficients$CardRenamed<-paste("Card", SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, sep="")
#   
#   SubjectInfo_OrderedForDA5HTCqMatrix_QCed2_CardCoefficients<-join(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2_CardCoefficients, Temp_Coefficients, by="CardRenamed", type="left")
#   
#   Temp_CorrectedForCard<-Temp[i,]-SubjectInfo_OrderedForDA5HTCqMatrix_QCed2_CardCoefficients$fixef.Model.
#   
#   DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[i,]<-tapply(Temp_CorrectedForCard, SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$ID, function(y) mean(y, na.rm=TRUE))
#   rm(Temp)
#   rm(Temp2)
#   rm(Model)
#   rm(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2_CardCoefficients)
#   rm(Temp_CorrectedForCard)
# }
# head(DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject)
# # 1        2        3        4        5        6        7        8        9       10
# # HPRT1 21.11097 21.35400 21.35100 20.98100 20.93100 21.15100 21.76700 21.59100 20.64139 20.91939
# # GUSB  23.95638 24.36062 24.67862 24.22762 24.06462 24.42262 24.64662 24.55262 23.80342 24.23192
# # B2M   20.81375 21.16475 21.23375 20.89975 20.85475 21.32675 21.31975 21.29875 20.26169 21.01719
# # TFRC  22.95162 23.00537 22.87037 23.26737 22.89337 23.17337 23.98237 23.07637 22.86259 22.59509
# # 
# 
# DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject<-DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[, order(colnames(DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject))]
# 
# #HPRT1:
# plot(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[1,]~DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[1,])
# #Some correlation, still mostly driven by a handful of extreme genes
# cor(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[1,],DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[1,])
# #[1] 0.8048282
# #Better than before (it was 0.71 before)
# 
# #GUSB:
# plot(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[2,]~DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[2,])
# #perhaps a little uglier?
# cor(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[2,],DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[2,])
# #[1] 0.5456544
# #Worse than before (it was 0.6146429)
# 
# #B2M
# plot(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[3,]~DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[3,])
# #pretty
# cor(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[3,],DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[3,])
# #[1] 0.7721487
# #Not as good as before: 0.8977061
# 
# #TFRC
# plot(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[4,]~DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[4,])
# #o.k.
# cor(GabaGlu_SharedHousekeepingGenes_AverageMinusCardBySubject[4,],DA5HT_SharedHousekeepingGenes_AverageMinusCardBySubject[4,])
# #[1] 0.6336549
# #Also not as good as before: 0.8427888

#Upon further thought, I think that I am subtracting out much of the variation by subject that isn't due to the variables included in the model.  
#...but hypothetically the estimates for the coefficients for the variables included in the model should be better (more replicable) if controlling for the technical variation is improving things. Does that hold true?

colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

TempForScaledOutput<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3[, c(1:2,6:9,11,24,26,32)]
colnames(TempForScaledOutput)
TempForScaledOutput[,4]<-scale(TempForScaledOutput[,4])
TempForScaledOutput[,6]<-scale(TempForScaledOutput[,6])
TempForScaledOutput[,7]<-scale(TempForScaledOutput[,7])
TempForScaledOutput[,8]<-scale(TempForScaledOutput[,8])
TempForScaledOutput[,9]<-scale(TempForScaledOutput[,9])
head(TempForScaledOutput)
# ID        Card Diagnosis        Age Gender         pH Hours.Final TZP_BioAnalyzer_RINTZP_Average.RNAConc..ng.uL. CardBlock2
# 97   2  Card 1.eds     Schiz -0.4491722      M -0.2996711  0.25663894       0.4718285             -0.68134932   OddCards
# 193  3  Card 1.eds   Control -2.2615231      M  1.0205538  0.02132677       0.1246256              0.35332425   OddCards
# 289  4  Card 1.eds        BP  0.2304594      F -0.0356261 -0.69845165       0.2982271              0.50356298   OddCards
# 385 17 Card 10.eds        BP  1.6652373      M  0.9545426  0.78263086       0.9926329              0.03093258  EvenCards
# 481 18 Card 10.eds     Schiz -0.3736575      M -1.4218622 -0.98913139      -2.1321933              0.44020845  EvenCards
# 577 19 Card 10.eds        BP  0.2304594      F -1.2238285 -1.75043549       0.8190315             -0.05812900  EvenCards

colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
TempForScaledOutput2<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[, c(1:2,6:9,11,24,26)]
colnames(TempForScaledOutput2)
TempForScaledOutput2[,4]<-scale(TempForScaledOutput2[,4])
TempForScaledOutput2[,6]<-scale(TempForScaledOutput2[,6])
TempForScaledOutput2[,7]<-scale(TempForScaledOutput2[,7])
TempForScaledOutput2[,8]<-scale(TempForScaledOutput2[,8])
TempForScaledOutput2[,9]<-scale(TempForScaledOutput2[,9])
head(TempForScaledOutput2)
# ID        Card Diagnosis        Age Gender            pH Hours.Final TZP_BioAnalyzer_RINTZP_Average.RNAConc..ng.uL.
# 1    1 Card 39.eds   Control  2.1194462      M -1.342291e-01 -0.73657200      -0.1811594               1.3875207
# 49   2 Card 39.eds     Schiz -0.4833825      M -2.685189e-01  0.27247504       0.5206004              -0.6851077
# 97   3 Card 39.eds   Control -2.2681793      M  1.074379e+00  0.03422782       0.1697205               0.3303203
# 145  4 Card 39.eds        BP  0.1859163      F  6.062195e-05 -0.69452838       0.3451604               0.4777645
# 193  5 Card 39.eds   Control -0.1115498      M -1.409982e+00 -0.10591760       0.6960404               0.2627676
# 241  6 Card 39.eds     Schiz -0.2602829      M  1.611538e+00 -0.24606302       0.5206004              -0.1393091


for(i in c(1:5)){
  Temp<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_GabaGlu<-fixef(Model)
  
  Temp<-DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_DA5HT<-fixef(Model)
  
  print(row.names(Temp)[i])
  Temp_Results<-data.frame(Temp_Coefficients_GabaGlu[c(2:9)], Temp_Coefficients_DA5HT[c(2:9)])
  colnames(Temp_Results)<-c("GabaGlu", "DA5HT")
  Temp_Results$ABS_Diff<-abs(Temp_Results[,1]-Temp_Results[,2])
  print(Temp_Results)
  print(sum(Temp_Results$ABS_Diff))
}

# [1] "HPRT1"
# GabaGlu        DA5HT   ABS_Diff
# DiagnosisBP              0.01615612  0.043586370 0.02743025
# DiagnosisSchiz           0.16440957  0.238046235 0.07363667
# Age                     -0.08359619 -0.069679921 0.01391627
# GenderF                  0.07459426  0.158251935 0.08365768
# pH                      -0.11575942 -0.146108541 0.03034913
# Hours.Final             -0.05789164  0.009952882 0.06784452
# TZP_BioAnalyzer_RIN         -0.26297701 -0.236802137 0.02617488
#TZP_Average.RNAConc..ng.uL. -0.15651077 -0.237603952 0.08109318
# [1] 0.4041026
# [1] "GUSB"
# GabaGlu       DA5HT    ABS_Diff
# DiagnosisBP              0.070224084 -0.01831737 0.088541449
# DiagnosisSchiz           0.102176144  0.04752418 0.054651961
# Age                     -0.036647859 -0.09268905 0.056041192
# GenderF                  0.030081140  0.03586568 0.005784538
# pH                       0.019743835  0.01820641 0.001537421
# Hours.Final              0.003821189  0.02864815 0.024826963
# TZP_BioAnalyzer_RIN         -0.066787414 -0.07525374 0.008466323
#TZP_Average.RNAConc..ng.uL. -0.162659515 -0.34533790 0.182678380
# [1] 0.4225282
# [1] "B2M"
# GabaGlu        DA5HT     ABS_Diff
# DiagnosisBP              0.033910830  0.007823377 2.608745e-02
# DiagnosisSchiz           0.172023065  0.128245179 4.377789e-02
# Age                     -0.114286396 -0.121782181 7.495785e-03
# GenderF                  0.144884635  0.042289680 1.025950e-01
# pH                       0.001085318  0.025907534 2.482222e-02
# Hours.Final              0.110953478  0.155474768 4.452129e-02
# TZP_BioAnalyzer_RIN         -0.072349076 -0.072379022 2.994525e-05
#TZP_Average.RNAConc..ng.uL. -0.164347015 -0.280998037 1.166510e-01
# [1] 0.3659806
# [1] "IPO8"
# GabaGlu       DA5HT    ABS_Diff
# DiagnosisBP              0.074489869  0.02426468 0.050225186
# DiagnosisSchiz           0.131177967  0.04455420 0.086623771
# Age                     -0.028305265 -0.02555793 0.002747331
# GenderF                  0.023224124 -0.03878835 0.062012469
# pH                      -0.008357002 -0.02792195 0.019564945
# Hours.Final             -0.033451019  0.01225136 0.045702380
# TZP_BioAnalyzer_RIN         -0.133963849 -0.11295351 0.021010340
#TZP_Average.RNAConc..ng.uL. -0.161558479 -0.27564758 0.114089098
# [1] 0.4019755
# [1] "TFRC"
# GabaGlu       DA5HT   ABS_Diff
# DiagnosisBP              0.217541568  0.13206005 0.08548152
# DiagnosisSchiz           0.155623478  0.05820778 0.09741570
# Age                     -0.001162099 -0.01845856 0.01729646
# GenderF                  0.010182762 -0.07622468 0.08640744
# pH                       0.038549225  0.07013366 0.03158444
# Hours.Final             -0.045678604 -0.00172413 0.04395447
# TZP_BioAnalyzer_RIN         -0.168044624 -0.13475094 0.03329369
#TZP_Average.RNAConc..ng.uL. -0.164880961 -0.27559242 0.11071146
# [1] 0.5061452

for(i in c(1:5)){
  Temp<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_GabaGlu<-fixef(Model)
  
  Temp<-DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_DA5HT<-fixef(Model)
  
  print(row.names(Temp)[i])
  Temp_Results<-data.frame(Temp_Coefficients_GabaGlu[c(2:9)], Temp_Coefficients_DA5HT[c(2:9)])
  colnames(Temp_Results)<-c("GabaGlu", "DA5HT")
  Temp_Results$ABS_Diff<-abs(Temp_Results[,1]-Temp_Results[,2])
  print(Temp_Results)
  print(sum(Temp_Results$ABS_Diff))
}

# [1] "HPRT1"
# GabaGlu       DA5HT    ABS_Diff
# DiagnosisBP             -0.01105059 -0.01750684 0.006456253
# DiagnosisSchiz           0.15348769  0.16451546 0.011027770
# Age                     -0.12244505 -0.08100377 0.041441288
# GenderF                 -0.11698517 -0.01576043 0.101224733
# pH                      -0.13922220 -0.12354851 0.015673690
# Hours.Final             -0.06673975 -0.04687048 0.019869269
# TZP_BioAnalyzer_RIN         -0.25765053 -0.24585660 0.011793934
#TZP_Average.RNAConc..ng.uL. -0.24192950 -0.40647915 0.164549647
# [1] 0.3720366
# [1] "GUSB"
# GabaGlu        DA5HT    ABS_Diff
# DiagnosisBP              0.129857519 -0.004541197 0.134398716
# DiagnosisSchiz           0.137462456  0.083080944 0.054381512
# Age                     -0.058692425 -0.088197760 0.029505336
# GenderF                 -0.008451706  0.055256806 0.063708512
# pH                       0.028172625  0.019354024 0.008818601
# Hours.Final              0.014890095  0.048447297 0.033557202
# TZP_BioAnalyzer_RIN         -0.111025227 -0.077466328 0.033558899
#TZP_Average.RNAConc..ng.uL. -0.120364326 -0.320995280 0.200630954
# [1] 0.5585597
# [1] "B2M"
# GabaGlu        DA5HT    ABS_Diff
# DiagnosisBP              0.08536097 -0.023121144 0.108482116
# DiagnosisSchiz           0.15255053  0.121442819 0.031107713
# Age                     -0.09484797 -0.119810417 0.024962447
# GenderF                  0.06718083  0.074833596 0.007652766
# pH                       0.01205885  0.009172672 0.002886180
# Hours.Final              0.06921853  0.158182975 0.088964441
# TZP_BioAnalyzer_RIN         -0.10560812 -0.047872906 0.057735215
#TZP_Average.RNAConc..ng.uL. -0.11671202 -0.306141190 0.189429169
# [1] 0.51122
# [1] "IPO8"
# GabaGlu         DA5HT    ABS_Diff
# DiagnosisBP              0.13160524  7.151814e-05 0.131533727
# DiagnosisSchiz           0.12089803  5.500766e-02 0.065890370
# Age                     -0.05519571 -3.494635e-02 0.020249365
# GenderF                  0.02837613 -1.047457e-02 0.038850706
# pH                      -0.01418843 -2.267234e-02 0.008483911
# Hours.Final             -0.03369970  4.159662e-03 0.037859362
# TZP_BioAnalyzer_RIN         -0.15113829 -1.261692e-01 0.024969127
#TZP_Average.RNAConc..ng.uL. -0.12032516 -2.733625e-01 0.153037287
# [1] 0.4808739
# [1] "TFRC"
# GabaGlu        DA5HT    ABS_Diff
# DiagnosisBP              0.23907393  0.170184375 0.068889554
# DiagnosisSchiz           0.16583181  0.047287009 0.118544796
# Age                     -0.01138582 -0.001650185 0.009735637
# GenderF                 -0.05703319 -0.165754800 0.108721610
# pH                       0.02232081  0.064973343 0.042652532
# Hours.Final             -0.05094705 -0.005941599 0.045005453
# TZP_BioAnalyzer_RIN         -0.24421186 -0.134462664 0.109749201
#TZP_Average.RNAConc..ng.uL. -0.11805533 -0.307456764 0.189401433
# [1] 0.6927002



for(i in c(1:5)){
  Temp<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+CardBlock2+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_GabaGlu<-fixef(Model)
  
  Temp<-DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_DA5HT<-fixef(Model)
  
  print(row.names(Temp)[i])
  Temp_Results<-data.frame(Temp_Coefficients_GabaGlu[c(2:9)], Temp_Coefficients_DA5HT[c(2:9)])
  colnames(Temp_Results)<-c("GabaGlu", "DA5HT")
  Temp_Results$ABS_Diff<-abs(Temp_Results[,1]-Temp_Results[,2])
  print(Temp_Results)
  print(sum(Temp_Results$ABS_Diff))
}

# [1] "HPRT1"
# GabaGlu       DA5HT    ABS_Diff
# DiagnosisBP              0.02839429 -0.01750684 0.045901132
# DiagnosisSchiz           0.17734836  0.16451546 0.012832894
# Age                     -0.08811318 -0.08100377 0.007109412
# GenderF                  0.07641472 -0.01576043 0.092175151
# pH                      -0.11347295 -0.12354851 0.010075558
# Hours.Final             -0.05418564 -0.04687048 0.007315166
# TZP_BioAnalyzer_RIN         -0.26468366 -0.24585660 0.018827065
#TZP_Average.RNAConc..ng.uL. -0.15529581 -0.40647915 0.251183337
# [1] 0.4454197
# [1] "GUSB"
# GabaGlu        DA5HT    ABS_Diff
# DiagnosisBP              0.09140312 -0.004541197 0.095944316
# DiagnosisSchiz           0.12587347  0.083080944 0.042792528
# Age                     -0.04553685 -0.088197760 0.042660909
# GenderF                  0.03457636  0.055256806 0.020680449
# pH                       0.02375489  0.019354024 0.004400868
# Hours.Final              0.01062787  0.048447297 0.037819425
# TZP_BioAnalyzer_RIN         -0.06969829 -0.077466328 0.007768040
#TZP_Average.RNAConc..ng.uL. -0.16107689 -0.320995280 0.159918386
# [1] 0.4119849
# [1] "B2M"
# GabaGlu        DA5HT    ABS_Diff
# DiagnosisBP              0.049096848 -0.023121144 0.072217992
# DiagnosisSchiz           0.190714532  0.121442819 0.069271713
# Age                     -0.122375001 -0.119810417 0.002564584
# GenderF                  0.145290613  0.074833596 0.070457017
# pH                       0.002313521  0.009172672 0.006859151
# Hours.Final              0.115741814  0.158182975 0.042441161
# TZP_BioAnalyzer_RIN         -0.074474716 -0.047872906 0.026601810
#TZP_Average.RNAConc..ng.uL. -0.162868294 -0.306141190 0.143272896
# [1] 0.4336863
# [1] "IPO8"
# GabaGlu         DA5HT    ABS_Diff
# DiagnosisBP              0.081137532  7.151814e-05 0.081066014
# DiagnosisSchiz           0.138725425  5.500766e-02 0.083717767
# Age                     -0.031058704 -3.494635e-02 0.003887643
# GenderF                  0.024133979 -1.047457e-02 0.034608551
# pH                      -0.007282331 -2.267234e-02 0.015390013
# Hours.Final             -0.031359075  4.159662e-03 0.035518737
# TZP_BioAnalyzer_RIN         -0.134925029 -1.261692e-01 0.008755869
#TZP_Average.RNAConc..ng.uL. -0.160748109 -2.733625e-01 0.112614341
# [1] 0.3755589
# [1] "TFRC"
# GabaGlu        DA5HT    ABS_Diff
# DiagnosisBP              0.228246492  0.170184375 0.058062117
# DiagnosisSchiz           0.166657574  0.047287009 0.119370565
# Age                     -0.004123214 -0.001650185 0.002473029
# GenderF                  0.011268457 -0.165754800 0.177023257
# pH                       0.040646559  0.064973343 0.024326784
# Hours.Final             -0.042554106 -0.005941599 0.036612507
# TZP_BioAnalyzer_RIN         -0.169751646 -0.134462664 0.035288982
#TZP_Average.RNAConc..ng.uL. -0.162662113 -0.307456764 0.144794651
# [1] 0.5979519


for(i in c(1:5)){
  Temp<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:84, 86, 88, 92),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+CardBlock2+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_GabaGlu<-fixef(Model)
  
  Temp<-DA5HT_Cq_AllSubjects_QCed3_Imputed[c(31:35),]
  Temp2<-data.frame(y=Temp[i,], TempForScaledOutput2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+BioAnalyzer_RIN+Average.RNAConc..ng.uL.+(1 | ID), data = Temp2, REML=F)
  Temp_Coefficients_DA5HT<-fixef(Model)
  
  print(row.names(Temp)[i])
  Temp_Results<-data.frame(Temp_Coefficients_GabaGlu[c(2:9)], Temp_Coefficients_DA5HT[c(2:9)])
  colnames(Temp_Results)<-c("GabaGlu", "DA5HT")
  Temp_Results$ABS_Diff<-abs(Temp_Results[,1]-Temp_Results[,2])
  print(Temp_Results)
  print(sum(Temp_Results$ABS_Diff))
}

# [1] "HPRT1"
# GabaGlu        DA5HT   ABS_Diff
# DiagnosisBP              0.02839429  0.043586370 0.01519208
# DiagnosisSchiz           0.17734836  0.238046235 0.06069788
# Age                     -0.08811318 -0.069679921 0.01843326
# GenderF                  0.07641472  0.158251935 0.08183722
# pH                      -0.11347295 -0.146108541 0.03263559
# Hours.Final             -0.05418564  0.009952882 0.06413853
# TZP_BioAnalyzer_RIN         -0.26468366 -0.236802137 0.02788153
#TZP_Average.RNAConc..ng.uL. -0.15529581 -0.237603952 0.08230814
# [1] 0.3831242
# [1] "GUSB"
# GabaGlu       DA5HT    ABS_Diff
# DiagnosisBP              0.09140312 -0.01831737 0.109720484
# DiagnosisSchiz           0.12587347  0.04752418 0.078349289
# Age                     -0.04553685 -0.09268905 0.047152199
# GenderF                  0.03457636  0.03586568 0.001289321
# pH                       0.02375489  0.01820641 0.005548478
# Hours.Final              0.01062787  0.02864815 0.018020280
# TZP_BioAnalyzer_RIN         -0.06969829 -0.07525374 0.005555449
#TZP_Average.RNAConc..ng.uL. -0.16107689 -0.34533790 0.184261002
# [1] 0.4498965
# [1] "B2M"
# GabaGlu        DA5HT     ABS_Diff
# DiagnosisBP              0.049096848  0.007823377 0.0412734708
# DiagnosisSchiz           0.190714532  0.128245179 0.0624693522
# Age                     -0.122375001 -0.121782181 0.0005928199
# GenderF                  0.145290613  0.042289680 0.1030009328
# pH                       0.002313521  0.025907534 0.0235940136
# Hours.Final              0.115741814  0.155474768 0.0397329547
# TZP_BioAnalyzer_RIN         -0.074474716 -0.072379022 0.0020956942
#TZP_Average.RNAConc..ng.uL. -0.162868294 -0.280998037 0.1181297425
# [1] 0.390889
# [1] "IPO8"
# GabaGlu       DA5HT   ABS_Diff
# DiagnosisBP              0.081137532  0.02426468 0.05687285
# DiagnosisSchiz           0.138725425  0.04455420 0.09417123
# Age                     -0.031058704 -0.02555793 0.00550077
# GenderF                  0.024133979 -0.03878835 0.06292232
# pH                      -0.007282331 -0.02792195 0.02063962
# Hours.Final             -0.031359075  0.01225136 0.04361044
# TZP_BioAnalyzer_RIN         -0.134925029 -0.11295351 0.02197152
#TZP_Average.RNAConc..ng.uL. -0.160748109 -0.27564758 0.11489947
# [1] 0.4205882
# [1] "TFRC"
# GabaGlu       DA5HT   ABS_Diff
# DiagnosisBP              0.228246492  0.13206005 0.09618644
# DiagnosisSchiz           0.166657574  0.05820778 0.10844979
# Age                     -0.004123214 -0.01845856 0.01433534
# GenderF                  0.011268457 -0.07622468 0.08749314
# pH                       0.040646559  0.07013366 0.02948710
# Hours.Final             -0.042554106 -0.00172413 0.04082998
# TZP_BioAnalyzer_RIN         -0.169751646 -0.13475094 0.03500071
#TZP_Average.RNAConc..ng.uL. -0.162662113 -0.27559242 0.11293031
# [1] 0.5247128

#Comparison summary:

# # The usual suspects (+RIN & RNAconc)
# [1] 0.4041026
# [1] 0.4225282
# [1] 0.3659806
# [1] 0.4019755
# [1] 0.5061452
# # 
# # The usual suspects (+RIN & RNAconc) + Card
# [1] 0.3720366
# [1] 0.5585597
# [1] 0.51122
# [1] 0.4808739
# [1] 0.6927002
# # 
# # The usual suspects (+RIN & RNAconc) + CardBlock2/Card
# [1] 0.4454197
# [1] 0.4119849
# [1] 0.4336863
# [1] 0.3755589
# [1] 0.5979519
# # 
# # The usual suspects (+RIN & RNAconc) + CardBlock2/NoCard
# [1] 0.3831242
# [1] 0.4498965
# [1] 0.390889
# [1] 0.4205882
# [1] 0.5247128


#Well, that is clear as mud. For some of the housekeeping genes, controlling for Card or CardBlock2 improves the similarity between the effects estimated in the two datasets, and for some genes it makes it worse.
#To make things more confusing, if there is actually *no real effect* introducing more variables should just increase noise.
#So let's look at two subject-level variables that I'm certain are likely to have real effects: RIN and RNAconc

# The usual suspects (+RIN & RNAconc)
# TZP_BioAnalyzer_RIN         -0.26297701 -0.236802137 0.02617488
#TZP_Average.RNAConc..ng.uL. -0.15651077 -0.237603952 0.08109318
# TZP_BioAnalyzer_RIN         -0.066787414 -0.07525374 0.008466323
#TZP_Average.RNAConc..ng.uL. -0.162659515 -0.34533790 0.182678380
# TZP_BioAnalyzer_RIN         -0.072349076 -0.072379022 2.994525e-05
#TZP_Average.RNAConc..ng.uL. -0.164347015 -0.280998037 1.166510e-01
# TZP_BioAnalyzer_RIN         -0.133963849 -0.11295351 0.021010340
#TZP_Average.RNAConc..ng.uL. -0.161558479 -0.27564758 0.114089098
# TZP_BioAnalyzer_RIN         -0.168044624 -0.13475094 0.03329369
#TZP_Average.RNAConc..ng.uL. -0.164880961 -0.27559242 0.11071146
# 
# The usual suspects (+RIN & RNAconc) + Card
# TZP_BioAnalyzer_RIN         -0.25765053 -0.24585660 0.011793934
#TZP_Average.RNAConc..ng.uL. -0.24192950 -0.40647915 0.164549647
# TZP_BioAnalyzer_RIN         -0.111025227 -0.077466328 0.033558899
#TZP_Average.RNAConc..ng.uL. -0.120364326 -0.320995280 0.200630954
# TZP_BioAnalyzer_RIN         -0.10560812 -0.047872906 0.057735215
#TZP_Average.RNAConc..ng.uL. -0.11671202 -0.306141190 0.189429169
# TZP_BioAnalyzer_RIN         -0.15113829 -1.261692e-01 0.024969127
#TZP_Average.RNAConc..ng.uL. -0.12032516 -2.733625e-01 0.153037287
# TZP_BioAnalyzer_RIN         -0.24421186 -0.134462664 0.109749201
#TZP_Average.RNAConc..ng.uL. -0.11805533 -0.307456764 0.189401433
# 
# The usual suspects (+RIN & RNAconc) + CardBlock2/Card
# TZP_BioAnalyzer_RIN         -0.26468366 -0.24585660 0.018827065
#TZP_Average.RNAConc..ng.uL. -0.15529581 -0.40647915 0.251183337
# TZP_BioAnalyzer_RIN         -0.06969829 -0.077466328 0.007768040
#TZP_Average.RNAConc..ng.uL. -0.16107689 -0.320995280 0.159918386
# TZP_BioAnalyzer_RIN         -0.074474716 -0.047872906 0.026601810
#TZP_Average.RNAConc..ng.uL. -0.162868294 -0.306141190 0.143272896
# TZP_BioAnalyzer_RIN         -0.134925029 -1.261692e-01 0.008755869
#TZP_Average.RNAConc..ng.uL. -0.160748109 -2.733625e-01 0.112614341
# TZP_BioAnalyzer_RIN         -0.169751646 -0.134462664 0.035288982
#TZP_Average.RNAConc..ng.uL. -0.162662113 -0.307456764 0.144794651

# 
# The usual suspects (+RIN & RNAconc) + CardBlock2/NoCard
# TZP_BioAnalyzer_RIN         -0.26468366 -0.236802137 0.02788153
#TZP_Average.RNAConc..ng.uL. -0.15529581 -0.237603952 0.08230814
# TZP_BioAnalyzer_RIN         -0.06969829 -0.07525374 0.005555449
#TZP_Average.RNAConc..ng.uL. -0.16107689 -0.34533790 0.184261002
# TZP_BioAnalyzer_RIN         -0.074474716 -0.072379022 0.0020956942
#TZP_Average.RNAConc..ng.uL. -0.162868294 -0.280998037 0.1181297425
# TZP_BioAnalyzer_RIN         -0.134925029 -0.11295351 0.02197152
#TZP_Average.RNAConc..ng.uL. -0.160748109 -0.27564758 0.11489947
# TZP_BioAnalyzer_RIN         -0.169751646 -0.13475094 0.03500071
#TZP_Average.RNAConc..ng.uL. -0.162662113 -0.27559242 0.11293031

#And that is, again, completely clear as mud. Ah well. It was a nice thought, right?



