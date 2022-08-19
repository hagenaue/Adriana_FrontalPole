
#Making some plots for the top diagnosis genes:

#Most people are used to looking at Negative Delta Delta Cq in plots - so let's plot things that way:

tapply(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))
# Control       BP    Schiz 
# 1.381810 1.475975 1.443336 

tapply(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))[[1]]
#[1] 1.38181

GabaGlu_NegDeltaCq_AverageForCTRLs<-matrix(0, 93, 1)

for(i in c(1:93)){
  GabaGlu_NegDeltaCq_AverageForCTRLs[i]<-tapply(GabaGlu_NegDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))[[1]]
}

GabaGlu_NegDeltaCq_AverageForCTRLs

GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2<-GabaGlu_NegDeltaCq_AllSubjects_QCed2

for(i in c(1:133)){
  GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[i,]<-GabaGlu_NegDeltaCq_AllSubjects_QCed2[i,]-GabaGlu_NegDeltaCq_AverageForCTRLs}

#Sanity Check:
GabaGlu_NegDeltaDeltaCq_AverageForCTRLs<-GabaGlu_NegDeltaCq_AverageForCTRLs

for(i in c(1:93)){
  GabaGlu_NegDeltaDeltaCq_AverageForCTRLs[i]<-tapply(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,i], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, function(y) mean(y, na.rm=T))[[1]]
}
GabaGlu_NegDeltaDeltaCq_AverageForCTRLs
#Perfect - basically 0 for everything.


pdf("Boxplot_GabaGlu_ABAT_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="ABAT", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("ABAT Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()


pdf("Boxplot_GabaGlu_SST_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SST", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SST Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_GPHN_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GPHN", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GPHN Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_MAPK1_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="MAPK1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("MAPK1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_GNAQ_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GNAQ"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GNAQ", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GNAQ"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GNAQ Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_NSF_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="NSF"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="NSF", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="NSF"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("NSF Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()



#How about some of the other nominally significant (or validation-important) genes:
#TFRC significant, GRIK1, HOMER1 nominally significant, PVALB & GAD1
# "SLC38A1" "PGK1"  

pdf("Boxplot_GabaGlu_TFRC_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="TFRC"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="TFRC", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="TFRC"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("TFRC Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_PVALB_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="PVALB"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="PVALB", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="PVALB"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("PVALB Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_GAD1_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GAD1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GAD1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GAD1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GAD1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_SLC38A1_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC38A1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SLC38A1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC38A1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SLC38A1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_GFAP_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GFAP"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GFAP", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GFAP"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GFAP Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_CALB1_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="CALB1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="CALB1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="CALB1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("CALB1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_GABBR2_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GABBR2"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GABBR2", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GABBR2"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GABBR2 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()


pdf("Boxplot_GabaGlu_GABRD_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GABRD"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GABRD", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GABRD"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GABRD Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

#GABRR1 failed QC, so it wasn't plotted as a validation gene

pdf("Boxplot_GabaGlu_AQP4_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="AQP4"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="AQP4", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="AQP4"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("AQP4 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_GJA1_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GJA1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="GJA1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GJA1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GJA1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_S100B_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="S100B"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="S100B", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="S100B"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("S100B Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_SLC1A2_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC1A2"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SLC1A2", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC1A2"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SLC1A2 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_SLC1A3_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC1A3"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SLC1A3", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC1A3"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SLC1A3 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

#This gene supposedly goes opposite directions in BP vs. Schiz in other datasets:
pdf("Boxplot_GabaGlu_SLC6A1_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC6A1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SLC6A1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC6A1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SLC6A1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_SLC6A11_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC6A11"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SLC6A11", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC6A11"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SLC6A11 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_SLC6A13_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC6A13"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="SLC6A13", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SLC6A13"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SLC6A13 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf("Boxplot_GabaGlu_TBP_ByDiagnosis.pdf", width=4, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="TBP"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), main="TBP", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="TBP"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("TBP Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()


#***********************************

#Out of curiousity, what are the top genes for some of the other variables?

#Age - there are a lot of genes related to age. Here are some of the top ones:
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,3]<0.000005]
#[1] "ADCY7"   "CALB1"   "DLG4"    "GABRA5"  "GRIK1"   "GRIK4"   "GRIN2B"  "HOMER1"  "HOMER2"  "SHANK2"  "SLC17A8"

#Gender
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,4]<0.10]
# "SLC17A6"

#pH:not much - I had to increase the cut-off for FDR
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,5]<0.05]
#[1] "GAD1"    "GNAI1"   "GRIN2C"  "GRM2"    "HOMER1"  "AQP4"    "GJA1"    "SLC1A3"  "SLC6A12" "SLC7A11" "HPRT1"  
#many astrocyte genes

#PMI: 
row.names(MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR)[MLM_GabaGlu_AllGenesByTheUsualSuspectsAndCard_NegDeltaCq_FDR[,6]<0.05]
# [1] "GRIA3"  "GRIA4"  "GRIK2"  "GRM8"   "AQP4"   "SHANK2" "SNCA" 



#############################################
