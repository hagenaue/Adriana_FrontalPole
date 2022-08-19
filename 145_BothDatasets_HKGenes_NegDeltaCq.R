
#Out of curiousity, five of the housekeeping genes are found in both datasets (HPRT1, GUSB, B2M, IPO8, TFRC) - how well do their values correlate?
#Note: my first round of analysis for this missed gene IPO8. Later when I came back to re-do it, I didn't re-do everything (because quite a few of the analyses turned out to be less than insightful...) so some of the analyses below still don't include it.

colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,c(83:84, 86, 88, 92)])
#[1] "HPRT1" "GUSB"  "B2M"   "IPO8"  "TFRC" 

dim(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,c(83:84, 86, 88, 92)])
#[1] 133   5
dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
#[1] 133  50

colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2[,c(31:35)])
#[1] "HPRT1" "GUSB"  "B2M"   "IPO8"  "TFRC" 
dim(DA5HT_NegDeltaCq_AllSubjects_QCed2[,c(31:35)])
#[1] 154   5
dim(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
#[1] 154  47

#Let's try a simple version of correlation vs. a more complex model controlling for sample-level technical co-variates.

#Calculating average housekeeping gene expression by subject for shared genes
GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq<-matrix(0,5, 69)
row.names(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,c(83:84, 86, 88, 92)])
colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID))

for(i in c(1:5)){
  Temp<-GabaGlu_NegDeltaCq_AllSubjects_QCed2[,c(83:84, 86, 88, 92)]
  GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[i,]<-tapply(Temp[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$ID, function(y) mean(y, na.rm=TRUE))
  rm(Temp)
}
head(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq)
# 1         10         13         14         15         16        17
# HPRT1  1.141575  1.1739091  0.9835771  1.2897273  1.2081818  1.3091818  1.324273
# GUSB  -2.230425 -2.1505909 -2.1719229 -2.4857727 -2.0383182 -2.1218182 -2.509727
# B2M    1.619575  1.0959091  1.7790771  1.1742273  1.3176818  1.0091818  1.421273
# IPO8  -1.574425 -1.4810909 -1.7019229 -1.5957727 -1.6603182 -1.4368182 -1.494227
# TFRC  -1.170425 -0.8485909 -1.1969229 -0.7252727 -0.6658182 -0.3828182 -1.161727

DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq<-matrix(0,5, 69)
row.names(DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq)<-colnames(DA5HT_NegDeltaCq_AllSubjects_QCed2[,c(31:35)])
colnames(DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq)<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$ID))
#Interesting - that isn't the same order as created by tapply in the other dataset.

for(i in c(1:5)){
  Temp<-DA5HT_NegDeltaCq_AllSubjects_QCed2[,c(31:35)]
  DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[i,]<-tapply(Temp[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$ID, function(y) mean(y, na.rm=TRUE))
  rm(Temp)
}
head(DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq)
# 1      10         13         14         15          16         17
# HPRT1  0.8348824  0.9710  0.8519412  1.2839118  1.1188824  1.03688235  1.0678824
# GUSB  -1.9768676 -2.0535 -1.8530588 -2.2230882 -1.9891176 -1.97611765 -2.2781176
# B2M    1.1936324  1.0175  1.8114412  1.0034118  1.3968824  0.84588235  1.2493824
# IPO8  -1.3693676 -1.3990 -1.5880588 -1.6030882 -1.6461176 -1.53161765 -1.5666176
# TFRC  -0.9288676 -0.2630 -0.6825588 -0.3320882 -0.3591176 -0.07511765 -0.8176176

DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq<-DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[,order(colnames(DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq))]

cbind(colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq), colnames(DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq))
#They're both in the same order. :)

SubjectInfo_forHKplots<-SubjectInfo[SubjectInfo$ID%in%colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq),]
dim(SubjectInfo_forHKplots)
#[1] 69 44

cbind(SubjectInfo_forHKplots$ID[order(as.character(SubjectInfo_forHKplots$ID))], colnames(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq))
#They are both the same order.

SubjectInfo_forHKplots<-SubjectInfo_forHKplots[order(as.character(SubjectInfo_forHKplots$ID)),]
SubjectInfo_forHKplots$ID

#Let's see how they correlate:

setwd("~/Documents/Microarray Gen/FrontalPole/Output/FinalOutput_AcrossDatasets/HKgene_CorrelationsBetweenDatasets_NegDeltaCq")

#HPRT1:
pdf("HPRT1_AverageNegDeltaCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[1,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[1,], xlab="DA-5HT Dataset: Average NegDeltaCq by Subject", ylab="GABA-GLU Dataset: Average NegDeltaCq by Subject", main="HPRT1", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[1,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[1,])
abline(BestFitLine)
dev.off()
#Some correlation, mostly driven by a handful of extreme samples
#There are some subjects with extreme values that really don't correlate with their values in the other dataset.
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[1,],DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[1,])
#[1] 0.4486994

#GUSB:
pdf("GUSB_AverageNegDeltaCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[2,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[2,], xlab="DA-5HT Dataset: Average NegDeltaCq by Subject", ylab="GABA-GLU Dataset: Average NegDeltaCq by Subject", main="GUSB", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[2,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[2,])
abline(BestFitLine)
dev.off()
#less driven by just a handful of extreme samples:
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[2,],DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[2,], use="pairwise.complete")
#[1] 0.3865809

#B2M
pdf("B2M_AverageNegDeltaCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[3,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[3,], xlab="DA-5HT Dataset: Average NegDeltaCq by Subject", ylab="GABA-GLU Dataset: Average NegDeltaCq by Subject", main="B2M", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[3,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[3,])
abline(BestFitLine)
dev.off()
#Still very pretty, oddly enough.
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[3,],DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[3,])
#[1] 0.8985143

pdf("IPO8_AverageNegDeltaCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[4,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[4,], xlab="DA-5HT Dataset: Average NegDeltaCq by Subject", ylab="GABA-GLU Dataset: Average NegDeltaCq by Subject", main="IPO8", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[4,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[4,])
abline(BestFitLine)
dev.off()
#Definitely noisier than for Cq
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[4,],DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[4,])
#[1] 0.4001958


#TFRC
pdf("TFRC_AverageNegDeltaCq_GabaGluVsDA5HT.pdf", width=5, height=5)
plot(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[5,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[5,], xlab="DA-5HT Dataset: Average NegDeltaCq by Subject", ylab="GABA-GLU Dataset: Average NegDeltaCq by Subject", main="TFRC", col=as.numeric(as.factor(SubjectInfo_forHKplots$Diagnosis))+1)
BestFitLine<-lm(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[5,]~DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[5,])
abline(BestFitLine)
dev.off()
#Still Very pretty
cor(GabaGlu_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[5,],DA5HT_SharedHousekeepingGenes_AverageBySubject_NegDeltaCq[5,])
#[1] 0.7453367


#So in general, subtracting out the sample-related variability in HK genes decreases the correlation between subject-level average HK expression in the two datasets. Not so unexpected, since presumably real subject-level variability in HK expression is also being subtracted out as well as variability introduced by technical factors at the sample level. 

