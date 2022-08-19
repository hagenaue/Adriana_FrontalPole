
#Making some plots for the top diagnosis genes:

#Most people are used to looking at Negative Delta Delta Cq in plots - so let's plot things that way:

tapply(DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))
# Control        BP     Schiz 
# -2.681232 -2.769004 -2.678466 

tapply(DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))[[1]]
#[1] -2.681232

DA5HT_NegDeltaCq_AverageForCTRLs<-matrix(0, 46, 1)

for(i in c(1:46)){
  DA5HT_NegDeltaCq_AverageForCTRLs[i]<-tapply(DA5HT_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))[[1]]
}

DA5HT_NegDeltaCq_AverageForCTRLs

DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2<-DA5HT_NegDeltaCq_AllSubjects_QCed2

for(i in c(1:154)){
  DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[i,]<-DA5HT_NegDeltaCq_AllSubjects_QCed2[i,]-DA5HT_NegDeltaCq_AverageForCTRLs}

#Sanity Check:
DA5HT_NegDeltaDeltaCq_AverageForCTRLs<-DA5HT_NegDeltaCq_AverageForCTRLs

for(i in c(1:46)){
  DA5HT_NegDeltaDeltaCq_AverageForCTRLs[i]<-tapply(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))[[1]]
}
DA5HT_NegDeltaDeltaCq_AverageForCTRLs
#Percent - basically 0 for everything.

pdf("Boxplot_DA5HT_HTR2B_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="HTR2B", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("HTR2B Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

#Double-checking (after finding out it is likely to be meningeal): relationship with cohort
pdf("Boxplot_DA5HT_HTR2B_ByCohort.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Cohort, las=2, cex.axis=0.75, cex.lab=0.75, pch=20, cex=1.7, main="HTR2B", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Cohort, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1)
mtext(expression(paste("HTR2B Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_HTR2B_ByDiagnosisCohort.pdf", width=10, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis*SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Cohort, las=2, cex.axis=0.5, cex.lab=0.5, pch=20, cex=1.7, main="HTR2B", col=rep(c(2,3,4),9), outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis*SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Cohort, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1)
mtext(expression(paste("HTR2B Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()


pdf("Boxplot_DA5HT_DRD4_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="DRD4", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("DRD4 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_DRD2_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD2"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="DRD2", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD2"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("DRD2 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_MAOB_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOB"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="MAOB", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOB"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("MAOB Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_TFRC_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="TFRC"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="TFRC", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="TFRC"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("TFRC Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

#Validation genes:

pdf("Boxplot_DA5HT_GAPDH_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="GAPDH"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GAPDH", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="GAPDH"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GAPDH Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_HTR6_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR6"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="HTR6", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR6"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("HTR6 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_MAOA_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOA"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="MAOA", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOA"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("MAOA Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_DA5HT_REEP5_ByDiagnosis.pdf", width=4, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="REEP5"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="REEP5", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="REEP5"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("REEP5 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()


#***********************************


#Out of curiousity, what are the top genes for some of the other variables?

#Age - much fewer genes relate to age in the DA5HT dataset vs. the GabaGlu dataset:
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,3]<0.005]
#[1] "DRD2"  "DRD4"  "HTR2A" "TH"    "GPI"   "REEP5"

#Gender - basically nothing even vaguely relates to gender:
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,4]<0.20]
#character(0)

#pH: fewer genes relate to pH than in the GabaGlu dataset:
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,5]<0.10]
# character(0)

#PMI: more than in the GabaGlu dataset, but still not much.
row.names(MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_DA5HT_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,6]<0.05]
# [1]  "DDC"   "HTR2A" "B2M"   "GPI" 


#**************************************************