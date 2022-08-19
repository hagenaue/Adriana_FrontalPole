#Reviewing the metadata for correlated variables:

##Note: I went down a serious rabbit hole here. Much of this code/analysis is pretty irrelevant for the actual Frontal Pole qPCR paper.

#############################

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/FullSample")


#Outputting some basic correlations:
colnames(SubjectInfo)

CorMatrix_RNAMetrics_Pearson<-cor(SubjectInfo[,NumericColumnsInSubjectInfo], use="pairwise.complete.obs")

write.csv(CorMatrix_RNAMetrics_Pearson, "CorMatrix_RNAMetrics_Pearson.csv")

#Double-checking how much the relationships are driven by outliers/distribution
CorMatrix_RNAMetrics_Spearman<-cor(SubjectInfo[,NumericColumnsInSubjectInfo], use="pairwise.complete.obs", method="spearman")

write.csv(CorMatrix_RNAMetrics_Spearman, "CorMatrix_RNAMetrics_Spearman.csv")

rm(CorMatrix_RNAMetrics_Pearson, CorMatrix_RNAMetrics_Spearman)

#With low quality RNA or imputed values removed:

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC")

CorMatrix_RNAMetrics_Pearson<-cor(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData[,NumericColumnsInSubjectInfo], use="pairwise.complete.obs")

write.csv(CorMatrix_RNAMetrics_Pearson, "CorMatrix_RNAMetrics_Pearson_NoBad.csv")

CorMatrix_RNAMetrics_Spearman<-cor(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData[,NumericColumnsInSubjectInfo], use="pairwise.complete.obs", method="spearman")

write.csv(CorMatrix_RNAMetrics_Spearman, "CorMatrix_RNAMetrics_Spearman_NoBad.csv")

pdf("Heatmap_CorMatrix_RNAMetrics_Pearson_NoBad.pdf", width=5, height=5)
#par(mar=c(5.1, 4.1, 8.1, 8.1))
heatmap(CorMatrix_RNAMetrics_Pearson, cexRow=0.3, cexCol=0.3)
dev.off()

pdf("Heatmap_CorMatrix_RNAMetrics_Spearman_NoBad.pdf", width=5, height=5)
#par(mar=c(5.1, 4.1, 8.1, 8.1))
heatmap(CorMatrix_RNAMetrics_Spearman, cexRow=0.3, cexCol=0.3)
dev.off()

rm(CorMatrix_RNAMetrics_Pearson, CorMatrix_RNAMetrics_Spearman)


setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/FullSample")


####################


#Getting into the weeds: Which variables can be consolidated?  Which variables are important for interpreting later gene expression measurements? (or including as co-variates?)


#The basics: Consolidating RNA extraction related-variables:

#Plotting the similarity between measurements made at MNI vs. the DSC:

colnames(SubjectInfo)

pdf("Scatterplot_MBNITTPConcbyDSCTTPConc.pdf", width=6, height=6)
plot(TTP_MBNI_Nanodrop_RNAConc..ng.uL.~TTP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo, xlab="TTP DSC RNA Concentration (ng/uL)", ylab="TTP MBNI RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_MBNI_Nanodrop_RNAConc..ng.uL.~TTP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITZPConcbyDSCTZPConc.pdf", width=6, height=6)
plot(TZP_MBNI_Nanodrop_RNAConc..ng.uL.~TZP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo, xlab="TZP DSC RNA Concentration (ng/uL)", ylab="TZP MBNI RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_MBNI_Nanodrop_RNAConc..ng.uL.~TZP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITTP260230byDSCTTP260230.pdf", width=6, height=6)
plot(TTP_MBNI_Nanodrop_260.230~TTP_DSC_Nanodrop_260.230, data=SubjectInfo, xlab="TTP RNA Purity: DSC 260/230", ylab="TTP RNA Purity: MBNI 260/230", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_MBNI_Nanodrop_260.230~TTP_DSC_Nanodrop_260.230, data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITZP260230byDSCTZP260230.pdf", width=6, height=6)
plot(TZP_MBNI_Nanodrop_260.230~TZP_DSC_Nanodrop_260.230, data=SubjectInfo, xlab="TZP RNA Purity: DSC 260/230", ylab="TZP RNA Purity: MBNI 260/230", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_MBNI_Nanodrop_260.230~TZP_DSC_Nanodrop_260.230, data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITTP260280byDSCTTP260280.pdf", width=6, height=6)
plot(TTP_MBNI_Nanodrop_260.280~TTP_DSC_Nanodrop_260.280, data=SubjectInfo, xlab="TTP RNA Purity: DSC 260/280", ylab="TTP RNA Purity: MBNI 260/280", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_MBNI_Nanodrop_260.280~TTP_DSC_Nanodrop_260.280, data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITZP260280byDSCTZP260280.pdf", width=6, height=6)
plot(TZP_MBNI_Nanodrop_260.280~TZP_DSC_Nanodrop_260.280, data=SubjectInfo, xlab="TZP RNA Purity: DSC 260/280", ylab="TZP RNA Purity: MBNI 260/280", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_MBNI_Nanodrop_260.280~TZP_DSC_Nanodrop_260.280, data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC")

#After removing the samples that failed QC or that were imputed (missing RNA extraction information):

pdf("Scatterplot_MBNITTPConcbyDSCTTPConc_noBad.pdf", width=6, height=6)
plot(TTP_MBNI_Nanodrop_RNAConc..ng.uL.~TTP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP DSC RNA Concentration (ng/uL)", ylab="TTP MBNI RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_MBNI_Nanodrop_RNAConc..ng.uL.~TTP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITZPConcbyDSCTZPConc_noBad.pdf", width=6, height=6)
plot(TZP_MBNI_Nanodrop_RNAConc..ng.uL.~TZP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TZP DSC RNA Concentration (ng/uL)", ylab="TZP MBNI RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_MBNI_Nanodrop_RNAConc..ng.uL.~TZP_DSC_Nanodrop_RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITTP260230byDSCTTP260230_noBad.pdf", width=6, height=6)
plot(TTP_MBNI_Nanodrop_260.230~TTP_DSC_Nanodrop_260.230, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Purity: DSC 260/230", ylab="TTP RNA Purity: MBNI 260/230", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_MBNI_Nanodrop_260.230~TTP_DSC_Nanodrop_260.230, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITZP260230byDSCTZP260230_noBad.pdf", width=6, height=6)
plot(TZP_MBNI_Nanodrop_260.230~TZP_DSC_Nanodrop_260.230, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TZP RNA Purity: DSC 260/230", ylab="TZP RNA Purity: MBNI 260/230", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_MBNI_Nanodrop_260.230~TZP_DSC_Nanodrop_260.230, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITTP260280byDSCTTP260280_noBad.pdf", width=6, height=6)
plot(TTP_MBNI_Nanodrop_260.280~TTP_DSC_Nanodrop_260.280, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Purity: DSC 260/280", ylab="TTP RNA Purity: MBNI 260/280", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_MBNI_Nanodrop_260.280~TTP_DSC_Nanodrop_260.280, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_MBNITZP260280byDSCTZP260280_noBad.pdf", width=6, height=6)
plot(TZP_MBNI_Nanodrop_260.280~TZP_DSC_Nanodrop_260.280, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TZP RNA Purity: DSC 260/280", ylab="TZP RNA Purity: MBNI 260/280", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_MBNI_Nanodrop_260.280~TZP_DSC_Nanodrop_260.280, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()


#MBNI and DSC measurements of RNA concentration and 260/230 can be easily averaged
#MBNI and DSC measurements of 260/280 are less similar, but talking to Adriana it isn't clear which of them is more "correct", so we'll need to just average them too.


#The BioAnalyzer RNA Integrity measures  (RIN and rRNA.Ratio..28s.18s.) are strongly related too and supposedly targeting the same thing - although instead of consolidating these measures, RIN is supposedly a better output, so we may just want to focus on that.


setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/FullSample")

pdf("Scatterplot_TTP_RIN_vs_TTP_28s18s.pdf", width=6, height=6)
plot(TTP_BioAnalyzer_RIN~TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo, xlab="TTP RNA Integrity (28s/18s rRNA)", ylab="TTP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_BioAnalyzer_RIN~TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TZP_RIN_vs_TZP_28s18s.pdf", width=6, height=6)
plot(TZP_BioAnalyzer_RIN~TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo, xlab="TZP RNA Integrity (28s/18s rRNA)", ylab="TZP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_BioAnalyzer_RIN~TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC")

pdf("Scatterplot_TTP_RIN_vs_TTP_28s18s_noBad.pdf", width=6, height=6)
plot(TTP_BioAnalyzer_RIN~TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Integrity (28s/18s rRNA)", ylab="TTP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_BioAnalyzer_RIN~TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TZP_RIN_vs_TZP_28s18s_noBad.pdf", width=6, height=6)
plot(TZP_BioAnalyzer_RIN~TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TZP RNA Integrity (28s/18s rRNA)", ylab="TZP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_BioAnalyzer_RIN~TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()


####

#Showing changes in concentration, purity, and concentration with RNA purification:


pdf("Scatterplot_TZPConcbyTTPConc.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL.~TTP_Average.RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP Average RNA Concentration (ng/uL)", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL.~TTP_Average.RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

pdf("Scatterplot_TZPConcbyTTPConc_ColbyDissectionGroup.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL.~TTP_Average.RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP Average RNA Concentration (ng/uL)", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7, col=as.factor(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$Dissecton.Group%in%c("9","10")))
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL.~TTP_Average.RNAConc..ng.uL., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()


pdf("Scatterplot_TZP260230byTTP260230.pdf", width=6, height=6)
plot(TZP_Average.260.230~TTP_Average.260.230, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Purity: Average 260/230", ylab="TZP RNA Purity: Average 260/230", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.260.230~TTP_Average.260.230, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red")
dev.off()

pdf("Scatterplot_TZP260280byTTP260280.pdf", width=6, height=6)
plot(TZP_Average.260.280~TTP_Average.260.280, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Purity: Average 260/280", ylab="TZP RNA Purity: Average 260/280", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.260.280~TTP_Average.260.280, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red")
dev.off()

pdf("Scatterplot_TZP_RINbyTTP_RIN.pdf", width=6, height=6)
plot(TZP_BioAnalyzer_RIN~TTP_BioAnalyzer_RIN, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Integrity Number (RIN)", ylab="TZP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_BioAnalyzer_RIN~TTP_BioAnalyzer_RIN, data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red")
dev.off()

pdf("Scatterplot_TZP_28s18s_byTTP_28s18s.pdf", width=6, height=6)
plot(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData, xlab="TTP RNA Integrity (28s/18s rRNA)", ylab="TZP RNA Integrity (28s/18s rRNA)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo_SurvivedBasicQC_NoMissingExtractionData)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red")
dev.off()


#There's no relationship between the TTP and TZP purity measurements - so let's plot overlying histograms instead:
#Using stolen code...

MinForHistograms<-min(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.260.230, SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.260.230), na.rm=T)
MaxForHistograms<-max(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.260.230, SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.260.230), na.rm=T)

ax <- pretty(c(MinForHistograms,MaxForHistograms), n=20)

hgA <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.260.230, breaks=ax, plot = FALSE) # Save first histogram data
hgB <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.260.230, breaks=ax, plot = FALSE) # Save 2nd histogram data

c1 <- rgb(173,216,230, max = 255, alpha = 150, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 150, names = "lt.pink")

pdf("OverlappingHistogram_TTPvsTZP_260230_noBad.pdf", width=6, height=6)
plot(hgA, col = c1, ylim=c(0,30), xlab="RNA Purity: 260/230", main="") # Plot 1st histogram using a transparent color
plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
dev.off()


MinForHistograms<-min(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.260.280, SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.260.280), na.rm=T)
MaxForHistograms<-max(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.260.280, SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.260.280), na.rm=T)

ax <- pretty(c(MinForHistograms,MaxForHistograms), n=20)

hgA <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.260.280, breaks=ax, plot = FALSE) # Save first histogram data
hgB <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.260.280, breaks=ax, plot = FALSE) # Save 2nd histogram data

c1 <- rgb(173,216,230, max = 255, alpha = 150, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 150, names = "lt.pink")

pdf("OverlappingHistogram_TTPvsTZP_260280_noBad.pdf", width=6, height=6)
plot(hgA, col = c1, ylim=c(0,30), xlab="RNA Purity: 260/280", main="") # Plot 1st histogram using a transparent color
plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
dev.off()


MinForHistograms<-min(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_BioAnalyzer_RIN, SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_BioAnalyzer_RIN), na.rm=T)
MaxForHistograms<-max(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_BioAnalyzer_RIN, SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_BioAnalyzer_RIN), na.rm=T)

ax <- pretty(c(MinForHistograms,MaxForHistograms), n=20)

hgA <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_BioAnalyzer_RIN, breaks=ax, plot = FALSE) # Save first histogram data
hgB <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_BioAnalyzer_RIN, breaks=ax, plot = FALSE) # Save 2nd histogram data

c1 <- rgb(173,216,230, max = 255, alpha = 150, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 150, names = "lt.pink")

pdf("OverlappingHistogram_TTPvsTZP_BioAnalyzer_RIN_noBad.pdf", width=6, height=6)
plot(hgA, col = c1, ylim=c(0,20), xlab="RNA Integrity (RIN)", main="") # Plot 1st histogram using a transparent color
plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
dev.off()


MinForHistograms<-min(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.RNAConc..ng.uL., SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.RNAConc..ng.uL.), na.rm=T)
MaxForHistograms<-max(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.RNAConc..ng.uL., SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.RNAConc..ng.uL.), na.rm=T)

ax <- pretty(c(MinForHistograms,MaxForHistograms), n=20)

hgA <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_Average.RNAConc..ng.uL., breaks=ax, plot = FALSE) # Save first histogram data
hgB <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_Average.RNAConc..ng.uL., breaks=ax, plot = FALSE) # Save 2nd histogram data

c1 <- rgb(173,216,230, max = 255, alpha = 150, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 150, names = "lt.pink")

pdf("OverlappingHistogram_TTPvsTZP_Average.RNAConc..ng.uL._noBad.pdf", width=6, height=6)
plot(hgA, col = c1, ylim=c(0,30), xlab="RNA Concentration (ng/uL)", main="") # Plot 1st histogram using a transparent color
plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
dev.off()


MinForHistograms<-min(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_BioAnalyzer_rRNA.Ratio..28s.18s., SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.), na.rm=T)
MaxForHistograms<-max(c(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_BioAnalyzer_rRNA.Ratio..28s.18s., SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.), na.rm=T)

ax <- pretty(c(MinForHistograms,MaxForHistograms), n=20)

hgA <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TTP_BioAnalyzer_rRNA.Ratio..28s.18s., breaks=ax, plot = FALSE) # Save first histogram data
hgB <- hist(SubjectInfo_SurvivedBasicQC_NoMissingExtractionData$TZP_BioAnalyzer_rRNA.Ratio..28s.18s., breaks=ax, plot = FALSE) # Save 2nd histogram data

c1 <- rgb(173,216,230, max = 255, alpha = 150, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 150, names = "lt.pink")

pdf("OverlappingHistogram_TTPvsTZP_BioAnalyzer_rRNA.Ratio..28s.18s._noBad.pdf", width=6, height=6)
plot(hgA, col = c1, ylim=c(0,20), xlab="RNA Integrity (28s/18s rRNA)", main="") # Plot 1st histogram using a transparent color
plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
dev.off()


###################################

#Quickly checking for diagnosis-related imbalance in the design for any of the numeric variables:

#Note: The formal evaluation of imbalance should probably be performed both with the NA substitutes and with the NA's removed.


setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/FullSample/Variables_vsDiagnosis")

colnames(SubjectInfo)
Temp<-SubjectInfo[,NumericColumnsInSubjectInfo]
colnames(Temp)
Temp$Diagnosis<-SubjectInfo$Diagnosis

DescriptiveStatistics<-matrix(0,28, 30)
colnames(DescriptiveStatistics)<-c("mean", "sd", "min", "max", "meanCTRL", "meanBP", "meanSCHIZ", "sdCTRL", "sdBP", "sdSCHIZ", "minCTRL", "minBP", "minSCHIZ", "maxCTRL", "maxBP", "maxSCHIZ", "BetaBP", "BetaSCHIZ", "SE_BP", "SE_SCHIZ", "TstatBP", "TstatSCHIZ", "PvalBP", "PvalSCHIZ", "ANOVA_SumOfSquares1", "ANOVA_SumOfSquares2", "ANOVA_df1", "ANOVA_df2", "ANOVA_Fstat","ANOVA_Pval")
row.names(DescriptiveStatistics)<-colnames(Temp)[c(1:28)]

for(i in c(1:28)){
  DescriptiveStatistics[i,1]<-mean(Temp[,i], na.rm=T)
  DescriptiveStatistics[i,2]<-sd(Temp[,i], na.rm=T)
  DescriptiveStatistics[i,3]<-min(Temp[,i], na.rm=T)
  DescriptiveStatistics[i,4]<-max(Temp[,i], na.rm=T)
  DescriptiveStatistics[i,c(5:7)]<-tapply(Temp[,i], SubjectInfo$Diagnosis, function(y) mean(y, na.rm=T))
  DescriptiveStatistics[i,c(8:10)]<-tapply(Temp[,i], SubjectInfo$Diagnosis, function(y) sd(y, na.rm=T))
  DescriptiveStatistics[i,c(11:13)]<-tapply(Temp[,i], SubjectInfo$Diagnosis, function(y) min(y, na.rm=T))
  DescriptiveStatistics[i,c(14:16)]<-tapply(Temp[,i], SubjectInfo$Diagnosis, function(y) max(y, na.rm=T))
  DemographicsLM<-summary.lm(lm(Temp[,i]~SubjectInfo$Diagnosis))
  DescriptiveStatistics[i,c(17:18)]<-DemographicsLM$coefficients[c(2:3),1]
  DescriptiveStatistics[i,c(19:20)]<-DemographicsLM$coefficients[c(2:3),2]
  DescriptiveStatistics[i,c(21:22)]<-DemographicsLM$coefficients[c(2:3),3]
  DescriptiveStatistics[i,c(23:24)]<-DemographicsLM$coefficients[c(2:3),4]
  ANOVAresults<-car::Anova(lm(Temp[,i]~Temp$Diagnosis), type=3)
  DescriptiveStatistics[i,c(25:26)]<-ANOVAresults$`Sum Sq`[c(2,3)]
  DescriptiveStatistics[i,c(27:28)]<-ANOVAresults$Df[c(2,3)]
  DescriptiveStatistics[i,29]<-ANOVAresults$`F value`[2]
  DescriptiveStatistics[i,30]<-ANOVAresults$`Pr(>F)`[2]
}

str(DescriptiveStatistics)
head(DescriptiveStatistics)

write.csv(DescriptiveStatistics, "DescriptiveStatistics.csv")


DescriptiveStatisticsNoNA<-DescriptiveStatistics

TempForPH<-Temp[(SubjectInfo$Subject.Number=="4463")==F,]

for(i in 2){
  DescriptiveStatisticsNoNA[i,1]<-mean(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,2]<-sd(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,3]<-min(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,4]<-max(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,c(5:7)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) mean(y, na.rm=T))
  DescriptiveStatisticsNoNA[i,c(8:10)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) sd(y, na.rm=T))
  DescriptiveStatisticsNoNA[i,c(11:13)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) min(y, na.rm=T))
  DescriptiveStatisticsNoNA[i,c(14:16)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) max(y, na.rm=T))
  DemographicsLM<-summary.lm(lm(TempForPH[,i]~TempForPH$Diagnosis))
  DescriptiveStatisticsNoNA[i,c(17:18)]<-DemographicsLM$coefficients[c(2:3),1]
  DescriptiveStatisticsNoNA[i,c(19:20)]<-DemographicsLM$coefficients[c(2:3),2]
  DescriptiveStatisticsNoNA[i,c(21:22)]<-DemographicsLM$coefficients[c(2:3),3]
  DescriptiveStatisticsNoNA[i,c(23:24)]<-DemographicsLM$coefficients[c(2:3),4]
  rm(DemographicsLM)
  ANOVAresults<-car::Anova(lm(TempForPH[,i]~TempForPH$Diagnosis), type=3)
  DescriptiveStatisticsNoNA[i,c(25:26)]<-ANOVAresults$`Sum Sq`[c(2,3)]
  DescriptiveStatisticsNoNA[i,c(27:28)]<-ANOVAresults$Df[c(2,3)]
  DescriptiveStatisticsNoNA[i,29]<-ANOVAresults$`F value`[2]
  DescriptiveStatisticsNoNA[i,30]<-ANOVAresults$`Pr(>F)`[2]
  rm(ANOVAresults)
}

TempForRNAQual<-Temp[(SubjectInfo$Subject.Number=="4463"|SubjectInfo$Subject.Number=="5162")==F,]

for(i in c(4:28)){
  DescriptiveStatisticsNoNA[i,1]<-mean(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,2]<-sd(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,3]<-min(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,4]<-max(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA[i,c(5:7)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) mean(y, na.rm=T))
  DescriptiveStatisticsNoNA[i,c(8:10)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) sd(y, na.rm=T))
  DescriptiveStatisticsNoNA[i,c(11:13)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) min(y, na.rm=T))
  DescriptiveStatisticsNoNA[i,c(14:16)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) max(y, na.rm=T))
  DemographicsLM<-summary.lm(lm(TempForRNAQual[,i]~TempForRNAQual$Diagnosis))
  DescriptiveStatisticsNoNA[i,c(17:18)]<-DemographicsLM$coefficients[c(2:3),1]
  DescriptiveStatisticsNoNA[i,c(19:20)]<-DemographicsLM$coefficients[c(2:3),2]
  DescriptiveStatisticsNoNA[i,c(21:22)]<-DemographicsLM$coefficients[c(2:3),3]
  DescriptiveStatisticsNoNA[i,c(23:24)]<-DemographicsLM$coefficients[c(2:3),4]
  rm(DemographicsLM)
  ANOVAresults<-car::Anova(lm(TempForRNAQual[,i]~TempForRNAQual$Diagnosis), type=3)
  DescriptiveStatisticsNoNA[i,c(25:26)]<-ANOVAresults$`Sum Sq`[c(2,3)]
  DescriptiveStatisticsNoNA[i,c(27:28)]<-ANOVAresults$Df[c(2,3)]
  DescriptiveStatisticsNoNA[i,29]<-ANOVAresults$`F value`[2]
  DescriptiveStatisticsNoNA[i,30]<-ANOVAresults$`Pr(>F)`[2]
  rm(ANOVAresults)
}

write.csv(DescriptiveStatisticsNoNA, "DescriptiveStatisticsNoNA.csv")

#None of the numeric variables differ by diagnosis, either with NAs replaced with grand means or without.


########################################

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC/Variables_vsDiagnosis")

colnames(SubjectInfo)
Temp<-SubjectInfo[SubjectInfo$LowQualityRNA==F,NumericColumnsInSubjectInfo]
colnames(Temp)
Temp$Diagnosis<-SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F]

DescriptiveStatistics_NoBad<-matrix(0,28, 30)
colnames(DescriptiveStatistics_NoBad)<-c("mean", "sd", "min", "max", "meanCTRL", "meanBP", "meanSCHIZ", "sdCTRL", "sdBP", "sdSCHIZ", "minCTRL", "minBP", "minSCHIZ", "maxCTRL", "maxBP", "maxSCHIZ", "BetaBP", "BetaSCHIZ", "SE_BP", "SE_SCHIZ", "TstatBP", "TstatSCHIZ", "PvalBP", "PvalSCHIZ", "ANOVA_SumOfSquares1", "ANOVA_SumOfSquares2", "ANOVA_df1", "ANOVA_df2", "ANOVA_Fstat","ANOVA_Pval")
row.names(DescriptiveStatistics_NoBad)<-colnames(Temp)[c(1:28)]

for(i in c(1:28)){
  DescriptiveStatistics_NoBad[i,1]<-mean(Temp[,i], na.rm=T)
  DescriptiveStatistics_NoBad[i,2]<-sd(Temp[,i], na.rm=T)
  DescriptiveStatistics_NoBad[i,3]<-min(Temp[,i], na.rm=T)
  DescriptiveStatistics_NoBad[i,4]<-max(Temp[,i], na.rm=T)
  DescriptiveStatistics_NoBad[i,c(5:7)]<-tapply(Temp[,i], Temp$Diagnosis, function(y) mean(y, na.rm=T))
  DescriptiveStatistics_NoBad[i,c(8:10)]<-tapply(Temp[,i], Temp$Diagnosis, function(y) sd(y, na.rm=T))
  DescriptiveStatistics_NoBad[i,c(11:13)]<-tapply(Temp[,i], Temp$Diagnosis, function(y) min(y, na.rm=T))
  DescriptiveStatistics_NoBad[i,c(14:16)]<-tapply(Temp[,i], Temp$Diagnosis, function(y) max(y, na.rm=T))
  DemographicsLM<-summary.lm(lm(Temp[,i]~Temp$Diagnosis))
  DescriptiveStatistics_NoBad[i,c(17:18)]<-DemographicsLM$coefficients[c(2:3),1]
  DescriptiveStatistics_NoBad[i,c(19:20)]<-DemographicsLM$coefficients[c(2:3),2]
  DescriptiveStatistics_NoBad[i,c(21:22)]<-DemographicsLM$coefficients[c(2:3),3]
  DescriptiveStatistics_NoBad[i,c(23:24)]<-DemographicsLM$coefficients[c(2:3),4]
  rm(DemographicsLM)
  ANOVAresults<-car::Anova(lm(Temp[,i]~Temp$Diagnosis), type=3)
  DescriptiveStatistics_NoBad[i,c(25:26)]<-ANOVAresults$`Sum Sq`[c(2,3)]
  DescriptiveStatistics_NoBad[i,c(27:28)]<-ANOVAresults$Df[c(2,3)]
  DescriptiveStatistics_NoBad[i,29]<-ANOVAresults$`F value`[2]
  DescriptiveStatistics_NoBad[i,30]<-ANOVAresults$`Pr(>F)`[2]
  rm(ANOVAresults)
}

str(DescriptiveStatistics_NoBad)
head(DescriptiveStatistics_NoBad)

write.csv(DescriptiveStatistics_NoBad, "DescriptiveStatistics_NoBad.csv")



DescriptiveStatisticsNoNA_NoBad<-DescriptiveStatistics_NoBad

TempForPH<-Temp[(SubjectInfo$Subject.Number[SubjectInfo$LowQualityRNA==F]=="4463")==F,]

for(i in 2){
  DescriptiveStatisticsNoNA_NoBad[i,1]<-mean(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,2]<-sd(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,3]<-min(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,4]<-max(TempForPH[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,c(5:7)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) mean(y, na.rm=T))
  DescriptiveStatisticsNoNA_NoBad[i,c(8:10)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) sd(y, na.rm=T))
  DescriptiveStatisticsNoNA_NoBad[i,c(11:13)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) min(y, na.rm=T))
  DescriptiveStatisticsNoNA_NoBad[i,c(14:16)]<-tapply(TempForPH[,i], TempForPH$Diagnosis, function(y) max(y, na.rm=T))
  DemographicsLM<-summary.lm(lm(TempForPH[,i]~TempForPH$Diagnosis))
  DescriptiveStatisticsNoNA_NoBad[i,c(17:18)]<-DemographicsLM$coefficients[c(2:3),1]
  DescriptiveStatisticsNoNA_NoBad[i,c(19:20)]<-DemographicsLM$coefficients[c(2:3),2]
  DescriptiveStatisticsNoNA_NoBad[i,c(21:22)]<-DemographicsLM$coefficients[c(2:3),3]
  DescriptiveStatisticsNoNA_NoBad[i,c(23:24)]<-DemographicsLM$coefficients[c(2:3),4]
  rm(DemographicsLM)
  ANOVAresults<-car::Anova(lm(TempForPH[,i]~TempForPH$Diagnosis), type=3)
  DescriptiveStatisticsNoNA_NoBad[i,c(25:26)]<-ANOVAresults$`Sum Sq`[c(2,3)]
  DescriptiveStatisticsNoNA_NoBad[i,c(27:28)]<-ANOVAresults$Df[c(2,3)]
  DescriptiveStatisticsNoNA_NoBad[i,29]<-ANOVAresults$`F value`[2]
  DescriptiveStatisticsNoNA_NoBad[i,30]<-ANOVAresults$`Pr(>F)`[2]
  rm(ANOVAresults)
}

TempForRNAQual<-Temp[(SubjectInfo$Subject.Number[SubjectInfo$LowQualityRNA==F]=="4463"|SubjectInfo$Subject.Number[SubjectInfo$LowQualityRNA==F]=="5162")==F,]

for(i in c(4:28)){
  DescriptiveStatisticsNoNA_NoBad[i,1]<-mean(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,2]<-sd(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,3]<-min(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,4]<-max(TempForRNAQual[,i], na.rm=T)
  DescriptiveStatisticsNoNA_NoBad[i,c(5:7)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) mean(y, na.rm=T))
  DescriptiveStatisticsNoNA_NoBad[i,c(8:10)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) sd(y, na.rm=T))
  DescriptiveStatisticsNoNA_NoBad[i,c(11:13)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) min(y, na.rm=T))
  DescriptiveStatisticsNoNA_NoBad[i,c(14:16)]<-tapply(TempForRNAQual[,i], TempForRNAQual$Diagnosis, function(y) max(y, na.rm=T))
  DemographicsLM<-summary.lm(lm(TempForRNAQual[,i]~TempForRNAQual$Diagnosis))
  DescriptiveStatisticsNoNA_NoBad[i,c(17:18)]<-DemographicsLM$coefficients[c(2:3),1]
  DescriptiveStatisticsNoNA_NoBad[i,c(19:20)]<-DemographicsLM$coefficients[c(2:3),2]
  DescriptiveStatisticsNoNA_NoBad[i,c(21:22)]<-DemographicsLM$coefficients[c(2:3),3]
  DescriptiveStatisticsNoNA_NoBad[i,c(23:24)]<-DemographicsLM$coefficients[c(2:3),4]
  rm(DemographicsLM)
  ANOVAresults<-car::Anova(lm(TempForRNAQual[,i]~TempForRNAQual$Diagnosis), type=3)
  DescriptiveStatisticsNoNA_NoBad[i,c(25:26)]<-ANOVAresults$`Sum Sq`[c(2,3)]
  DescriptiveStatisticsNoNA_NoBad[i,c(27:28)]<-ANOVAresults$Df[c(2,3)]
  DescriptiveStatisticsNoNA_NoBad[i,29]<-ANOVAresults$`F value`[2]
  DescriptiveStatisticsNoNA_NoBad[i,30]<-ANOVAresults$`Pr(>F)`[2]
  rm(ANOVAresults)
}

write.csv(DescriptiveStatisticsNoNA_NoBad, "DescriptiveStatisticsNoNA_NoBad.csv")

#None of the numeric variables differ by diagnosis, either with NAs replaced with grand means or without.


#Making a prettier version for publication (or, at least, the basis for a prettier version, to be further prettified by hand in Excel):
colnames(DescriptiveStatisticsNoNA_NoBad)

Demographics_WrittenFormat<-paste(
row.names(DescriptiveStatisticsNoNA_NoBad),
paste("mean (+/-sd): ", signif(DescriptiveStatisticsNoNA_NoBad[,1], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,2], 2), sep=""),
paste("range: ", signif(DescriptiveStatisticsNoNA_NoBad[,3], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,4], 3), sep=""),
paste("CTRL mean (+/-sd): ", signif(DescriptiveStatisticsNoNA_NoBad[,5], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,8], 2), sep=""),
paste("CTRL range: ", signif(DescriptiveStatisticsNoNA_NoBad[,11], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,14], 3), sep=""),
paste("BP mean (+/-sd): ", signif(DescriptiveStatisticsNoNA_NoBad[,6], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,9], 2), sep=""),
paste("BP range: ", signif(DescriptiveStatisticsNoNA_NoBad[,12], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,15], 3), sep=""),
paste("SCHIZ mean (+/-sd): ", signif(DescriptiveStatisticsNoNA_NoBad[,7], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,10], 2), sep=""),
paste("SCHIZ range: ", signif(DescriptiveStatisticsNoNA_NoBad[,13], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,16], 3), sep=""),
paste("Effect of Diagnosis: F(", DescriptiveStatisticsNoNA_NoBad[,27], ",", DescriptiveStatisticsNoNA_NoBad[,28], ")=",  signif(DescriptiveStatisticsNoNA_NoBad[,29], 3), ", p=", signif(DescriptiveStatisticsNoNA_NoBad[,30], 3), sep=""),
paste("Effect of BP: Beta(+/-SE)=", signif(DescriptiveStatisticsNoNA_NoBad[,17], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,19], 2), ", T(67)=", signif(DescriptiveStatisticsNoNA_NoBad[,21], 3), ", p=", signif(DescriptiveStatisticsNoNA_NoBad[,23], 3), sep=""),
paste("Effect of SCHIZ: Beta(+/-SE)=", signif(DescriptiveStatisticsNoNA_NoBad[,18], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,20], 2), ", T(67)=", signif(DescriptiveStatisticsNoNA_NoBad[,22], 3), ", p=", signif(DescriptiveStatisticsNoNA_NoBad[,24], 3), sep=""), sep=", ")

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables")

write.csv(Demographics_WrittenFormat, "Demographics_WrittenFormat.csv")

Demographics_TableFormat<-data.frame(
  Variable=row.names(DescriptiveStatisticsNoNA_NoBad),
  meansd=paste(signif(DescriptiveStatisticsNoNA_NoBad[,1], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,2], 2), sep=""),
  range=paste(signif(DescriptiveStatisticsNoNA_NoBad[,3], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,4], 3), sep=""),
  CTRLmeansd=paste(signif(DescriptiveStatisticsNoNA_NoBad[,5], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,8], 2), sep=""),
  CTRLrange=paste(signif(DescriptiveStatisticsNoNA_NoBad[,11], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,14], 3), sep=""),
  BPmeansd=paste(signif(DescriptiveStatisticsNoNA_NoBad[,6], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,9], 2), sep=""),
  BPrange=paste(signif(DescriptiveStatisticsNoNA_NoBad[,12], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,15], 3), sep=""),
  SCHIZmeansd=paste(signif(DescriptiveStatisticsNoNA_NoBad[,7], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,10], 2), sep=""),
  SCHIZrange=paste(signif(DescriptiveStatisticsNoNA_NoBad[,13], 3), "-", signif(DescriptiveStatisticsNoNA_NoBad[,16], 3), sep=""),
  EffectOfDiagnosis=paste("F(", DescriptiveStatisticsNoNA_NoBad[,27], ",", DescriptiveStatisticsNoNA_NoBad[,28], ")=",  signif(DescriptiveStatisticsNoNA_NoBad[,29], 3), ", p=", signif(DescriptiveStatisticsNoNA_NoBad[,30], 3), sep=""),
  EffectOfBP=paste("Beta(+/-SE)=", signif(DescriptiveStatisticsNoNA_NoBad[,17], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,19], 2), ", T(67)=", signif(DescriptiveStatisticsNoNA_NoBad[,21], 3), ", p=", signif(DescriptiveStatisticsNoNA_NoBad[,23], 3), sep=""),
  EffectOfSCHIZ=paste("Beta(+/-SE)=", signif(DescriptiveStatisticsNoNA_NoBad[,18], 3), " +/-", signif(DescriptiveStatisticsNoNA_NoBad[,20], 2), ", T(67)=", signif(DescriptiveStatisticsNoNA_NoBad[,22], 3), ", p=", signif(DescriptiveStatisticsNoNA_NoBad[,24], 3), sep=""))
  
write.csv(Demographics_TableFormat, "Demographics_TableFormat.csv")

  
##########################

#Outputting some more basic summary demographics & evaluating multicollinearity:

#Full sample size:
length(SubjectInfo$Diagnosis)
#[1] 72

#Sample size with subjects removed that failed QC:

length(SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F])
#[1] 69


table(SubjectInfo$Diagnosis)
# Control BP   Schiz 
# 27      21      24

table(SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F])
# Control BP   Schiz 
# 26      21      22 


#For the following analyses, subjects that failed basic QC or were missing the relevant variables are removed:

table(SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F], SubjectInfo$Gender[SubjectInfo$LowQualityRNA==F])
#         M  F
# Control 25  1
# BP      15  6
# Schiz   21  1

fisher.test(SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F], SubjectInfo$Gender[SubjectInfo$LowQualityRNA==F])
# Fisher's Exact Test for Count Data
# 
# data:  SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA == F] and SubjectInfo$Gender[SubjectInfo$LowQualityRNA == F]
# p-value = 0.01735
# alternative hypothesis: two.sided

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC/Variables_vsDiagnosis")

pdf("Boxplot_pH_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(pH~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,], xlab="Diagnosis", ylab="Brain pH", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(pH~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,], xlab="Diagnosis", ylab="Brain pH", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(pH~Diagnosis, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,]))
#Multiple R-squared:  0.04936

#Checking for collinearity with diagnosis in a larger model:

summary.lm(lm(pH~Diagnosis+Age+Gender+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.25362 -0.10536 -0.00621  0.10755  0.32286 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     6.8979804  0.0942541  73.185   <2e-16 ***
# DiagnosisBP    -0.0715345  0.0479693  -1.491    0.141    
# DiagnosisSchiz  0.0040476  0.0454300   0.089    0.929    
# Age            -0.0016265  0.0014125  -1.152    0.254    
# GenderF        -0.0296036  0.0621963  -0.476    0.636    
# Hours.Final     0.0009007  0.0026537   0.339    0.735    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1525 on 62 degrees of freedom
# Multiple R-squared:  0.07487,	Adjusted R-squared:  0.000262 
# F-statistic: 1.004 on 5 and 62 DF,  p-value: 0.4233

car::Anova(lm(pH~Diagnosis+Age+Gender+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: pH
# Sum Sq Df   F value Pr(>F)    
# (Intercept) 124.584  1 5356.0311 <2e-16 ***
#   Diagnosis     0.066  2    1.4285 0.2474    
# Age           0.031  1    1.3260 0.2539    
# Gender        0.005  1    0.2265 0.6358    
# Hours.Final   0.003  1    0.1152 0.7354    
# Residuals     1.442 62                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_Age_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(Age~Diagnosis, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,], xlab="Diagnosis", ylab="Age", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(Age~Diagnosis, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,], xlab="Diagnosis", ylab="Age", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(Age~Diagnosis, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,]))
#Multiple R-squared:  0.03222

summary.lm(lm(Age~Diagnosis+pH+Gender+Hours.Final, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,]))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -31.450  -6.869  -0.830   8.400  29.110 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    132.0878    76.9101   1.717   0.0908 .
# DiagnosisBP     -6.5852     4.2860  -1.536   0.1294  
# DiagnosisSchiz  -3.4494     3.9954  -0.863   0.3912  
# pH             -12.4669    11.2379  -1.109   0.2715  
# GenderF          0.5981     5.5737   0.107   0.9149  
# Hours.Final      0.1935     0.2357   0.821   0.4150  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13.64 on 63 degrees of freedom
# Multiple R-squared:  0.05998,	Adjusted R-squared:  -0.01463 
# F-statistic: 0.8039 on 5 and 63 DF,  p-value: 0.5512


car::Anova(lm(Age~Diagnosis+pH+Gender+Hours.Final, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: Age
# Sum Sq Df F value  Pr(>F)  
# (Intercept)   549.1  1  2.9496 0.09081 .
# Diagnosis     451.2  2  1.2118 0.30450  
# pH            229.1  1  1.2307 0.27149  
# Gender          2.1  1  0.0115 0.91489  
# Hours.Final   125.4  1  0.6734 0.41498  
# Residuals   11728.7 63                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_PMI_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(Hours.Final~Diagnosis, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,], xlab="Diagnosis", ylab="Post-Mortem Interval (PMI, hrs)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(Hours.Final~Diagnosis, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,], xlab="Diagnosis", ylab="Post-Mortem Interval (PMI, hrs)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


summary.lm(lm(Hours.Final~Diagnosis+Age+pH+Gender, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,]))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -13.902  -4.613  -0.690   4.031  23.920 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)     7.12782   41.82086   0.170    0.865
# DiagnosisBP     1.04105    2.31695   0.449    0.655
# DiagnosisSchiz -2.24101    2.11767  -1.058    0.294
# Age             0.05467    0.06662   0.821    0.415
# pH              1.89265    6.02727   0.314    0.755
# GenderF        -4.03383    2.91928  -1.382    0.172
# 
# Residual standard error: 7.253 on 63 degrees of freedom
# Multiple R-squared:  0.06283,	Adjusted R-squared:  -0.01155 
# F-statistic: 0.8447 on 5 and 63 DF,  p-value: 0.5232

car::Anova(lm(Hours.Final~Diagnosis+Age+pH+Gender, data=SubjectInfo[SubjectInfo$LowQualityRNA==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: Hours.Final
# Sum Sq Df F value Pr(>F)
# (Intercept)    1.5  1  0.0290 0.8652
# Diagnosis    112.8  2  1.0717 0.3486
# Age           35.4  1  0.6734 0.4150
# pH             5.2  1  0.0986 0.7545
# Gender       100.4  1  1.9093 0.1719
# Residuals   3314.3 63 


table(SubjectInfo$Cohort, SubjectInfo$Diagnosis)

#                 Control BP Schiz
# Cohort 11            6  1     3
# Cohort 12            9  4     2
# Cohort 13            5  0     2
# Cohort 7             0  6     6
# Cohort 8             0  2     4
# Dep Cohort 1         2  2     0
# Dep Cohort 5         1  4     0
# Dep Cohort 6         0  2     2
# Schiz Cohort 2       4  0     5

table(SubjectInfo$Cohort[SubjectInfo$LowQualityRNA==F], SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F])

#                 Control BP Schiz
# Cohort 11            5  1     3
# Cohort 12            9  4     2
# Cohort 13            5  0     2
# Cohort 7             0  6     5
# Cohort 8             0  2     4
# Dep Cohort 1         2  2     0
# Dep Cohort 5         1  4     0
# Dep Cohort 6         0  2     1
# Schiz Cohort 2       4  0     5

table(SubjectInfo$AFS, SubjectInfo$Diagnosis)

#     Control BP Schiz
# 0      27 21    24

table(SubjectInfo$Dissecton.Group, SubjectInfo$Diagnosis)

#         Control BP Schiz
# 1             3  1     2
# 2             2  2     1
# 3             2  2     2
# 4             2  2     2
# 5             2  2     2
# 6             2  0     2
# 7             3  1     2
# 8             2  2     2
# 9             2  2     2
# 10            2  2     2
# 11            2  2     2
# 12            3  2     2
# 6redone       0  1     1


table(SubjectInfo$Dissecton.Group[SubjectInfo$LowQualityRNA==F], SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F])
# Control BP Schiz
# 1             3  1     2
# 2             1  2     1
# 3             2  2     2
# 4             2  2     2
# 5             2  2     2
# 6             2  0     2
# 7             3  1     2
# 8             2  2     1
# 9             2  2     1
# 10            2  2     2
# 11            2  2     2
# 12            3  2     2
# 6redone       0  1     1

write.csv(table(SubjectInfo$Dissecton.Group[SubjectInfo$LowQualityRNA==F], SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F]), "Table_DissectionGroupVsDiagnosis_noBadRIN.csv")



pdf("Boxplot_BlockWeight_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(Block.Weight..g.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="Block Weight (g)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(Block.Weight..g.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="Block Weight (g)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


summary.lm(lm(Block.Weight..g.~Diagnosis+Age+Gender+pH+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.37707 -0.13287 -0.01944  0.10384  0.39341 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)    -0.0617824  0.9927578  -0.062  0.95058   
# DiagnosisBP     0.0375526  0.0558956   0.672  0.50427   
# DiagnosisSchiz  0.0606273  0.0512100   1.184  0.24112   
# Age             0.0007722  0.0016096   0.480  0.63315   
# GenderF         0.2057877  0.0704544   2.921  0.00491 **
#   pH              0.1120357  0.1430940   0.783  0.43673   
# Hours.Final    -0.0027875  0.0030520  -0.913  0.36472   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1718 on 60 degrees of freedom
# Multiple R-squared:  0.1875,	Adjusted R-squared:  0.1063 
# F-statistic: 2.308 on 6 and 60 DF,  p-value: 0.04527

summary.lm(lm(Block.Weight..g.~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32859 -0.10752 -0.01866  0.10733  0.45029 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        0.007428   1.084358   0.007   0.9946  
# DiagnosisBP        0.019271   0.060045   0.321   0.7496  
# DiagnosisSchiz     0.063579   0.053636   1.185   0.2416  
# Age                0.000633   0.001686   0.375   0.7089  
# GenderF            0.203745   0.087409   2.331   0.0239 *
# pH                 0.098325   0.156589   0.628   0.5330  
# Hours.Final       -0.002603   0.003481  -0.748   0.4582  
# Dissecton.Group2   0.134053   0.118045   1.136   0.2616  
# Dissecton.Group3  -0.063539   0.104956  -0.605   0.5477  
# Dissecton.Group4   0.064685   0.108250   0.598   0.5529  
# Dissecton.Group5   0.048939   0.105368   0.464   0.6444  
# Dissecton.Group6   0.018273   0.117414   0.156   0.8770  
# Dissecton.Group7  -0.035417   0.103568  -0.342   0.7338  
# Dissecton.Group8   0.124499   0.111279   1.119   0.2687  
# Dissecton.Group9   0.119036   0.110738   1.075   0.2877  
# Dissecton.Group10  0.009058   0.103395   0.088   0.9305  
# Dissecton.Group11  0.040326   0.105133   0.384   0.7030  
# Dissecton.Group12 -0.011711   0.099515  -0.118   0.9068  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1773 on 49 degrees of freedom
# Multiple R-squared:  0.2934,	Adjusted R-squared:  0.04825 
# F-statistic: 1.197 on 17 and 49 DF,  p-value: 0.3021

car::Anova(lm(Block.Weight..g.~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: Block.Weight..g.
# Sum Sq Df F value  Pr(>F)  
# (Intercept)     0.00000  1  0.0000 0.99456  
# Diagnosis       0.04478  2  0.7120 0.49565  
# Age             0.00443  1  0.1410 0.70892  
# Gender          0.17085  1  5.4333 0.02392 *
# pH              0.01240  1  0.3943 0.53297  
# Hours.Final     0.01758  1  0.5592 0.45816  
# Dissecton.Group 0.23090 11  0.6675 0.76152  
# Residuals       1.54082 49                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP_RIN_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TZP_BioAnalyzer_RIN~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TZP_BioAnalyzer_RIN~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Integrity Number (RIN)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_BioAnalyzer_RIN~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
#Multiple R-squared:  0.02291

summary.lm(lm(TZP_BioAnalyzer_RIN~Diagnosis+Age+Gender+pH+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.4477 -0.2235  0.1162  0.3298  0.9106 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     3.846187   3.158598   1.218   0.2281  
# DiagnosisBP    -0.025505   0.177840  -0.143   0.8864  
# DiagnosisSchiz -0.207421   0.162932  -1.273   0.2079  
# Age             0.005328   0.005121   1.040   0.3023  
# GenderF         0.483385   0.224161   2.156   0.0351 *
# pH              0.613641   0.455274   1.348   0.1828  
# Hours.Final    -0.020700   0.009710  -2.132   0.0371 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5467 on 60 degrees of freedom
# Multiple R-squared:  0.1902,	Adjusted R-squared:  0.1092 
# F-statistic: 2.349 on 6 and 60 DF,  p-value: 0.04198

summary.lm(lm(TZP_BioAnalyzer_RIN~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.42756 -0.26808  0.08467  0.24164  1.01001 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        3.471325   2.941516   1.180 0.243652    
# DiagnosisBP       -0.033457   0.162884  -0.205 0.838106    
# DiagnosisSchiz    -0.180757   0.145496  -1.242 0.220019    
# Age                0.007226   0.004573   1.580 0.120499    
# GenderF            0.250872   0.237113   1.058 0.295230    
# pH                 0.669853   0.424775   1.577 0.121240    
# Hours.Final       -0.018037   0.009442  -1.910 0.061942 .  
# Dissecton.Group2  -0.200505   0.320219  -0.626 0.534124    
# Dissecton.Group3  -0.066545   0.284712  -0.234 0.816171    
# Dissecton.Group4   0.354565   0.293649   1.207 0.233057    
# Dissecton.Group5  -0.209185   0.285831  -0.732 0.467746    
# Dissecton.Group6  -1.273176   0.318506  -3.997 0.000215 ***
# Dissecton.Group7  -0.056927   0.280946  -0.203 0.840266    
# Dissecton.Group8  -0.043233   0.301864  -0.143 0.886703    
# Dissecton.Group9  -0.234593   0.300397  -0.781 0.438591    
# Dissecton.Group10  0.036895   0.280477   0.132 0.895885    
# Dissecton.Group11 -0.123015   0.285193  -0.431 0.668114    
# Dissecton.Group12 -0.208580   0.269951  -0.773 0.443435    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


car::Anova(lm(TZP_BioAnalyzer_RIN~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_BioAnalyzer_RIN
# Sum Sq Df F value  Pr(>F)  
# (Intercept)      0.3223  1  1.3927 0.24365  
# Diagnosis        0.3760  2  0.8125 0.44963  
# Age              0.5778  1  2.4970 0.12050  
# Gender           0.2590  1  1.1194 0.29523  
# pH               0.5754  1  2.4868 0.12124  
# Hours.Final      0.8445  1  3.6496 0.06194 .
# Dissecton.Group  6.5965 11  2.5916 0.01110 *
# Residuals       11.3383 49                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP_rRNA28s18s_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Integrity (28s/18s rRNA ratio)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Integrity (28s/18s rRNA ratio)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


summary.lm(lm(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis+Age+Gender+pH+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.12294 -0.19584  0.02404  0.23533  0.77708 
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    -2.388689   2.255001  -1.059   0.2937  
# DiagnosisBP     0.069404   0.126964   0.547   0.5867  
# DiagnosisSchiz -0.065689   0.116321  -0.565   0.5744  
# Age            -0.003416   0.003656  -0.934   0.3538  
# GenderF         0.014569   0.160034   0.091   0.9278  
# pH              0.757980   0.325031   2.332   0.0231 *
# Hours.Final    -0.004283   0.006933  -0.618   0.5390  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3903 on 60 degrees of freedom
# Multiple R-squared:  0.1174,	Adjusted R-squared:  0.02911 
# F-statistic:  1.33 on 6 and 60 DF,  p-value: 0.2581

summary.lm(lm(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.87617 -0.15594 -0.00508  0.23984  0.51713 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)       -2.506780   2.036431  -1.231   0.2242  
# DiagnosisBP       -0.017166   0.112765  -0.152   0.8796  
# DiagnosisSchiz    -0.036550   0.100728  -0.363   0.7183  
# Age               -0.002442   0.003166  -0.771   0.4441  
# GenderF            0.141776   0.164155   0.864   0.3920  
# pH                 0.738354   0.294075   2.511   0.0154 *
# Hours.Final        0.002566   0.006537   0.393   0.6964  
# Dissecton.Group2  -0.271626   0.221690  -1.225   0.2263  
# Dissecton.Group3   0.286883   0.197108   1.455   0.1519  
# Dissecton.Group4  -0.065372   0.203295  -0.322   0.7492  
# Dissecton.Group5   0.143288   0.197882   0.724   0.4724  
# Dissecton.Group6  -0.528177   0.220504  -2.395   0.0205 *
# Dissecton.Group7  -0.160222   0.194501  -0.824   0.4141  
# Dissecton.Group8   0.505027   0.208982   2.417   0.0194 *
# Dissecton.Group9  -0.019568   0.207967  -0.094   0.9254  
# Dissecton.Group10  0.235502   0.194176   1.213   0.2310  
# Dissecton.Group11  0.252432   0.197441   1.279   0.2071  
# Dissecton.Group12  0.073004   0.186889   0.391   0.6978  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.333 on 49 degrees of freedom
# Multiple R-squared:  0.4753,	Adjusted R-squared:  0.2932 
# F-statistic: 2.611 on 17 and 49 DF,  p-value: 0.004527

car::Anova(lm(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TZP_BioAnalyzer_rRNA.Ratio..28s.18s.
# Sum Sq Df F value   Pr(>F)   
# (Intercept)     0.1681  1  1.5153 0.224211   
# Diagnosis       0.0146  2  0.0660 0.936245   
# Age             0.0660  1  0.5952 0.444137   
# Gender          0.0827  1  0.7459 0.391976   
# pH              0.6991  1  6.3040 0.015396 * 
# Hours.Final     0.0171  1  0.1541 0.696361   
# Dissecton.Group 3.7068 11  3.0385 0.003631 **
# Residuals       5.4343 49                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP_RNAConc_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TZP_Average.RNAConc..ng.uL.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TZP_Average.RNAConc..ng.uL.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Concentration (ng/uL)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
#Multiple R-squared:  0.02248

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Diagnosis+Age+Gender+pH+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -216.314  -29.978    5.199   35.635  105.844 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    526.0010   376.0338   1.399    0.167
# DiagnosisBP     15.1670    21.1720   0.716    0.477
# DiagnosisSchiz -18.5726    19.3972  -0.957    0.342
# Age             -0.2404     0.6097  -0.394    0.695
# GenderF        -30.2313    26.6865  -1.133    0.262
# pH             -10.0136    54.2007  -0.185    0.854
# Hours.Final     -1.7906     1.1560  -1.549    0.127
# 
# Residual standard error: 65.09 on 60 degrees of freedom
# Multiple R-squared:  0.07737,	Adjusted R-squared:  -0.0149 
# F-statistic: 0.8385 on 6 and 60 DF,  p-value: 0.5452

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -194.961  -22.483    5.044   32.990  126.766 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)       564.18056  382.92344   1.473    0.147
# DiagnosisBP         8.49676   21.20401   0.401    0.690
# DiagnosisSchiz    -16.41587   18.94052  -0.867    0.390
# Age                -0.05351    0.59532  -0.090    0.929
# GenderF           -21.55825   30.86705  -0.698    0.488
# pH                -16.13176   55.29676  -0.292    0.772
# Hours.Final        -1.29371    1.22911  -1.053    0.298
# Dissecton.Group2  -52.93246   41.68582  -1.270    0.210
# Dissecton.Group3  -32.54423   37.06355  -0.878    0.384
# Dissecton.Group4  -37.09684   38.22690  -0.970    0.337
# Dissecton.Group5   29.26204   37.20912   0.786    0.435
# Dissecton.Group6  -60.59907   41.46275  -1.462    0.150
# Dissecton.Group7  -53.44257   36.57330  -1.461    0.150
# Dissecton.Group8  -10.31860   39.29634  -0.263    0.794
# Dissecton.Group9   22.25489   39.10531   0.569    0.572
# Dissecton.Group10  26.56613   36.51219   0.728    0.470
# Dissecton.Group11   1.14417   37.12607   0.031    0.976
# Dissecton.Group12 -40.33924   35.14195  -1.148    0.257
# 
# Residual standard error: 62.62 on 49 degrees of freedom
# Multiple R-squared:  0.3026,	Adjusted R-squared:  0.06061 
# F-statistic:  1.25 on 17 and 49 DF,  p-value: 0.264

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value Pr(>F)
# (Intercept)       8512  1  2.1708 0.1471
# Diagnosis         5374  2  0.6852 0.5088
# Age                 32  1  0.0081 0.9287
# Gender            1913  1  0.4878 0.4882
# pH                 334  1  0.0851 0.7717
# Hours.Final       4344  1  1.1079 0.2977
# Dissecton.Group  62047 11  1.4384 0.1864
# Residuals       192145 49

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value  Pr(>F)  
# (Intercept)        8489  1  2.3261 0.13379  
# Diagnosis          8325  2  1.1406 0.32814  
# Age                 159  1  0.0435 0.83567  
# Gender             6832  1  1.8721 0.17761  
# pH                  890  1  0.2439 0.62365  
# Hours.Final        2673  1  0.7324 0.39635  
# Block.Weight..g.  16973  1  4.6509 0.03608 *
# Dissecton.Group   58353 11  1.4536 0.18080  
# Residuals        175172 48                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP_260280_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TZP_Average.260.280~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Purity (260/280)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TZP_Average.260.280~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Purity (260/280)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_Average.260.280~Diagnosis+Age+Gender+pH+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.052931 -0.013980 -0.000168  0.011889  0.054440 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     2.160e+00  1.268e-01  17.026   <2e-16 ***
# DiagnosisBP    -1.289e-02  7.142e-03  -1.804   0.0762 .  
# DiagnosisSchiz -6.823e-03  6.543e-03  -1.043   0.3013    
# Age             9.317e-05  2.057e-04   0.453   0.6522    
# GenderF         1.358e-02  9.002e-03   1.509   0.1366    
# pH             -2.054e-02  1.828e-02  -1.124   0.2656    
# Hours.Final     2.191e-04  3.900e-04   0.562   0.5763    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02196 on 60 degrees of freedom
# Multiple R-squared:  0.09794,	Adjusted R-squared:  0.007731 
# F-statistic: 1.086 on 6 and 60 DF,  p-value: 0.3813


summary.lm(lm(TZP_Average.260.280~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.038669 -0.012423 -0.000440  0.008134  0.044774 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.159e+00  1.335e-01  16.165   <2e-16 ***
# DiagnosisBP       -9.518e-03  7.395e-03  -1.287   0.2041    
# DiagnosisSchiz    -9.151e-03  6.605e-03  -1.385   0.1722    
# Age                6.966e-05  2.076e-04   0.336   0.7387    
# GenderF            2.819e-03  1.076e-02   0.262   0.7945    
# pH                -2.064e-02  1.928e-02  -1.070   0.2897    
# Hours.Final       -1.697e-04  4.286e-04  -0.396   0.6939    
# Dissecton.Group2   1.852e-02  1.454e-02   1.274   0.2088    
# Dissecton.Group3   1.253e-02  1.293e-02   0.969   0.3372    
# Dissecton.Group4   2.958e-02  1.333e-02   2.218   0.0312 *  
# Dissecton.Group5   9.896e-03  1.298e-02   0.763   0.4493    
# Dissecton.Group6   2.894e-02  1.446e-02   2.001   0.0509 .  
# Dissecton.Group7   2.108e-02  1.275e-02   1.653   0.1047    
# Dissecton.Group8   8.720e-04  1.370e-02   0.064   0.9495    
# Dissecton.Group9  -2.099e-03  1.364e-02  -0.154   0.8783    
# Dissecton.Group10  1.467e-02  1.273e-02   1.152   0.2548    
# Dissecton.Group11  3.413e-03  1.295e-02   0.264   0.7932    
# Dissecton.Group12  1.318e-02  1.226e-02   1.076   0.2874    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02184 on 49 degrees of freedom
# Multiple R-squared:  0.2712,	Adjusted R-squared:  0.01837 
# F-statistic: 1.073 on 17 and 49 DF,  p-value: 0.4051

car::Anova(lm(TZP_Average.260.280~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TZP_Average.260.280
# Sum Sq Df  F value Pr(>F)    
# (Intercept)     0.124616  1 261.2933 <2e-16 ***
# Diagnosis       0.001240  2   1.2996 0.2819    
# Age             0.000054  1   0.1126 0.7387    
# Gender          0.000033  1   0.0686 0.7945    
# pH              0.000546  1   1.1456 0.2897    
# Hours.Final     0.000075  1   0.1567 0.6939    
# Dissecton.Group 0.005556 11   1.0591 0.4124    
# Residuals       0.023369 49                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP_260230_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TZP_Average.260.230~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Purity (260/230)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TZP_Average.260.230~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TZP RNA Purity (260/230)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


summary.lm(lm(TZP_Average.260.230~Diagnosis+Age+Gender+pH+Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.43016 -0.11771  0.04608  0.15207  0.23348 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     2.5181109  1.0723170   2.348   0.0222 *
#   DiagnosisBP    -0.0243879  0.0603751  -0.404   0.6877  
# DiagnosisSchiz  0.0262787  0.0553139   0.475   0.6365  
# Age            -0.0010240  0.0017386  -0.589   0.5581  
# GenderF         0.0551801  0.0761006   0.725   0.4712  
# pH             -0.0592494  0.1545615  -0.383   0.7028  
# Hours.Final     0.0005882  0.0032966   0.178   0.8590  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1856 on 60 degrees of freedom
# Multiple R-squared:  0.02159,	Adjusted R-squared:  -0.07625 
# F-statistic: 0.2206 on 6 and 60 DF,  p-value: 0.9687

summary.lm(lm(TZP_Average.260.230~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49325 -0.08872  0.01388  0.10889  0.34681 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        1.9322932  1.0731601   1.801   0.0779 .
# DiagnosisBP       -0.0006422  0.0594252  -0.011   0.9914  
# DiagnosisSchiz     0.0106234  0.0530816   0.200   0.8422  
# Age               -0.0007403  0.0016684  -0.444   0.6592  
# GenderF           -0.0518322  0.0865063  -0.599   0.5518  
# pH                 0.0395542  0.1549716   0.255   0.7996  
# Hours.Final       -0.0028049  0.0034446  -0.814   0.4194  
# Dissecton.Group2   0.1970090  0.1168264   1.686   0.0981 .
# Dissecton.Group3  -0.1906488  0.1038722  -1.835   0.0725 .
# Dissecton.Group4   0.0806189  0.1071326   0.753   0.4553  
# Dissecton.Group5  -0.0325878  0.1042802  -0.313   0.7560  
# Dissecton.Group6   0.0189346  0.1162012   0.163   0.8712  
# Dissecton.Group7  -0.0936118  0.1024983  -0.913   0.3656  
# Dissecton.Group8  -0.2016586  0.1101297  -1.831   0.0732 .
# Dissecton.Group9  -0.0066860  0.1095944  -0.061   0.9516  
# Dissecton.Group10  0.0117836  0.1023270   0.115   0.9088  
# Dissecton.Group11  0.0346595  0.1040475   0.333   0.7405  
# Dissecton.Group12  0.0210007  0.0984869   0.213   0.8320  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1755 on 49 degrees of freedom
# Multiple R-squared:  0.2857,	Adjusted R-squared:  0.03783 
# F-statistic: 1.153 on 17 and 49 DF,  p-value: 0.3363

car::Anova(lm(TZP_Average.260.230~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TZP_Average.260.230
# Sum Sq Df F value  Pr(>F)  
# (Intercept)     0.09985  1  3.2420 0.07793 .
# Diagnosis       0.00150  2  0.0244 0.97592  
# Age             0.00606  1  0.1969 0.65922  
# Gender          0.01106  1  0.3590 0.55182  
# pH              0.00201  1  0.0651 0.79961  
# Hours.Final     0.02042  1  0.6630 0.41943  
# Dissecton.Group 0.55791 11  1.6468 0.11505  
# Residuals       1.50916 49                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

############


##Getting into the (more esoteric) newer TTP data and detailed analyses of relationships between RNA extraction related variables.


#Making some illustrative plots:


pdf("Boxplot_TTP_Conc_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TTP_Average.RNAConc..ng.uL.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TTP_Average.RNAConc..ng.uL.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP Average RNA Concentration (ng/uL)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

pdf("Boxplot_TTP_260280_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TTP_Average.260.280~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP Average RNA Purity (260/280)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TTP_Average.260.280~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP Average RNA Purity (260/280)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

pdf("Boxplot_TTP_260230_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TTP_Average.260.230~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP Average RNA Purity (260/230)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TTP_Average.260.230~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP Average RNA Purity (260/230)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

pdf("Boxplot_TTP_rRNA28s18s_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TTP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP RNA Integrity (28s/18s rRNA ratio)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TTP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP RNA Integrity (28s/18s rRNA ratio)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

pdf("Boxplot_TTP_RIN_ByDiagnosis_noBad.pdf", width=4, height=6)
boxplot(TTP_BioAnalyzer_RIN~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(TTP_BioAnalyzer_RIN~Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Diagnosis", ylab="TTP RNA Integrity Number (RIN)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()



#####

#Going down more of a rabbithole trying to interpret RNA concentration:  Quite a few of the RNA extraction related measurements relate to dissection group. Let's plot that:

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC")


pdf("Boxplot_Weight_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(Block.Weight..g.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="Tissue Block Weight (g)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(Block.Weight..g.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="Tissue Block Weight (g)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(Block.Weight..g.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))    

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32833 -0.11200 -0.02857  0.13383  0.42333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.715000   0.074831   9.555 2.78e-13 ***
#   Dissecton.Group2   0.135000   0.118319   1.141    0.259    
# Dissecton.Group3  -0.096667   0.105828  -0.913    0.365    
# Dissecton.Group4   0.123333   0.105828   1.165    0.249    
# Dissecton.Group5   0.001667   0.105828   0.016    0.987    
# Dissecton.Group6  -0.012500   0.118319  -0.106    0.916    
# Dissecton.Group7  -0.068333   0.105828  -0.646    0.521    
# Dissecton.Group8   0.099000   0.110993   0.892    0.376    
# Dissecton.Group9   0.089000   0.110993   0.802    0.426    
# Dissecton.Group10  0.048333   0.105828   0.457    0.650    
# Dissecton.Group11 -0.005000   0.105828  -0.047    0.962    
# Dissecton.Group12 -0.026429   0.101978  -0.259    0.796    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1833 on 55 degrees of freedom
# Multiple R-squared:  0.1526,	Adjusted R-squared:  -0.01693 
# F-statistic: 0.9001 on 11 and 55 DF,  p-value: 0.5461


summary.lm(lm(Block.Weight..g.~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))  

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32859 -0.10752 -0.01866  0.10733  0.45029 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        0.007428   1.084358   0.007   0.9946  
# GenderF            0.203745   0.087409   2.331   0.0239 *
# Age                0.000633   0.001686   0.375   0.7089  
# Hours.Final       -0.002603   0.003481  -0.748   0.4582  
# pH                 0.098325   0.156589   0.628   0.5330  
# DiagnosisBP        0.019271   0.060045   0.321   0.7496  
# DiagnosisSchiz     0.063579   0.053636   1.185   0.2416  
# Dissecton.Group2   0.134053   0.118045   1.136   0.2616  
# Dissecton.Group3  -0.063539   0.104956  -0.605   0.5477  
# Dissecton.Group4   0.064685   0.108250   0.598   0.5529  
# Dissecton.Group5   0.048939   0.105368   0.464   0.6444  
# Dissecton.Group6   0.018273   0.117414   0.156   0.8770  
# Dissecton.Group7  -0.035417   0.103568  -0.342   0.7338  
# Dissecton.Group8   0.124499   0.111279   1.119   0.2687  
# Dissecton.Group9   0.119036   0.110738   1.075   0.2877  
# Dissecton.Group10  0.009058   0.103395   0.088   0.9305  
# Dissecton.Group11  0.040326   0.105133   0.384   0.7030  
# Dissecton.Group12 -0.011711   0.099515  -0.118   0.9068  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1773 on 49 degrees of freedom
# Multiple R-squared:  0.2934,	Adjusted R-squared:  0.04825 
# F-statistic: 1.197 on 17 and 49 DF,  p-value: 0.3021

car::Anova(lm(Block.Weight..g.~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: Block.Weight..g.
# Sum Sq Df F value  Pr(>F)  
# (Intercept)     0.00000  1  0.0000 0.99456  
# Gender          0.17085  1  5.4333 0.02392 *
# Age             0.00443  1  0.1410 0.70892  
# Hours.Final     0.01758  1  0.5592 0.45816  
# pH              0.01240  1  0.3943 0.53297  
# Diagnosis       0.04478  2  0.7120 0.49565  
# Dissecton.Group 0.23398 12  0.6201 0.81478  
# Residuals       1.54082 49                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Interestingly, Block weight is *not* one of the variables that differs by dissection group.

summary.lm(lm(Block.Weight..g.~Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.192486 -0.076268 -0.007249  0.082738  0.288545 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           1.511471   0.860018   1.757  0.08535 .  
# GenderF                               0.178581   0.066449   2.688  0.00992 ** 
# Age                                   0.002613   0.001336   1.956  0.05647 .  
# Hours.Final                          -0.006082   0.002767  -2.198  0.03291 *  
# pH                                   -0.118005   0.125026  -0.944  0.35008    
# DiagnosisBP                           0.021905   0.045085   0.486  0.62933    
# DiagnosisSchiz                        0.054518   0.040916   1.332  0.18913    
# TZP_BioAnalyzer_RIN                  -0.133945   0.045396  -2.951  0.00493 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.414507   0.065572   6.321 8.77e-08 ***
# Dissecton.Group2                      0.219788   0.089941   2.444  0.01835 *  
# Dissecton.Group3                     -0.191368   0.081386  -2.351  0.02295 *  
# Dissecton.Group4                      0.139274   0.083341   1.671  0.10134    
# Dissecton.Group5                     -0.038474   0.080745  -0.476  0.63594    
# Dissecton.Group6                      0.066671   0.101665   0.656  0.51515    
# Dissecton.Group7                      0.023371   0.078307   0.298  0.76667    
# Dissecton.Group8                     -0.090629   0.090218  -1.005  0.32025    
# Dissecton.Group9                      0.095725   0.083718   1.143  0.25866    
# Dissecton.Group10                    -0.083618   0.078977  -1.059  0.29512    
# Dissecton.Group11                    -0.080786   0.081374  -0.993  0.32590    
# Dissecton.Group12                    -0.069910   0.075731  -0.923  0.36065    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1331 on 47 degrees of freedom
# Multiple R-squared:  0.6182,	Adjusted R-squared:  0.4639 
# F-statistic: 4.006 on 19 and 47 DF,  p-value: 5.5e-05

car::Anova(lm(Block.Weight..g.~Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: Block.Weight..g.
# Sum Sq Df F value    Pr(>F)    
# (Intercept)                          0.05471  1  3.0888  0.085347 .  
# Gender                               0.12793  1  7.2227  0.009924 ** 
# Age                                  0.06774  1  3.8245  0.056468 .  
# Hours.Final                          0.08557  1  4.8310  0.032914 *  
# pH                                   0.01578  1  0.8908  0.350075    
# Diagnosis                            0.03145  2  0.8879  0.418320    
# TZP_BioAnalyzer_RIN                  0.15420  1  8.7061  0.004933 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 0.70777  1 39.9606 8.765e-08 ***
# Dissecton.Group                      0.52215 12  2.4567  0.014061 *  
#   Residuals                            0.83245 47                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


summary.lm(lm(Block.Weight..g.~Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           1.0951744  0.8344738   1.312  0.19589    
# GenderF                               0.1963999  0.0635508   3.090  0.00339 ** 
# Age                                   0.0026527  0.0012698   2.089  0.04226 *  
# Hours.Final                          -0.0051988  0.0026538  -1.959  0.05619 .  
# pH                                   -0.0962210  0.1191239  -0.808  0.42340    
# DiagnosisBP                           0.0155330  0.0429165   0.362  0.71906    
# DiagnosisSchiz                        0.0651819  0.0391172   1.666  0.10245    
# TZP_BioAnalyzer_RIN                  -0.1376667  0.0431598  -3.190  0.00256 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.4039805  0.0624503   6.469 5.72e-08 ***
# TZP_Average.RNAConc..ng.uL.           0.0007140  0.0002901   2.462  0.01764 *  
# Dissecton.Group2                      0.2539763  0.0865801   2.933  0.00521 ** 
# Dissecton.Group3                     -0.1653587  0.0780482  -2.119  0.03955 *  
# Dissecton.Group4                      0.1663932  0.0799500   2.081  0.04301 *  
# Dissecton.Group5                     -0.0586376  0.0771572  -0.760  0.45115    
# Dissecton.Group6                      0.0996404  0.0975221   1.022  0.31226    
# Dissecton.Group7                      0.0596313  0.0758487   0.786  0.43579    
# Dissecton.Group8                     -0.0781064  0.0858724  -0.910  0.36779    
# Dissecton.Group9                      0.0787552  0.0798444   0.986  0.32912    
# Dissecton.Group10                    -0.0999699  0.0753342  -1.327  0.19105    
# Dissecton.Group11                    -0.0794038  0.0773210  -1.027  0.30982    
# Dissecton.Group12                    -0.0411155  0.0729014  -0.564  0.57550    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1265 on 46 degrees of freedom
# Multiple R-squared:  0.6627,	Adjusted R-squared:  0.516 
# F-statistic: 4.518 on 20 and 46 DF,  p-value: 1.184e-05

car::Anova(lm(Block.Weight..g.~Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: Block.Weight..g.
# Sum Sq Df F value    Pr(>F)    
# (Intercept)                          0.02754  1  1.7224  0.195892    
# Gender                               0.15272  1  9.5508  0.003388 ** 
# Age                                  0.06979  1  4.3642  0.042262 *  
# Hours.Final                          0.06137  1  3.8377  0.056189 .  
# pH                                   0.01043  1  0.6524  0.423398    
# Diagnosis                            0.04555  2  1.4243  0.251091    
# TZP_BioAnalyzer_RIN                  0.16269  1 10.1742  0.002565 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 0.66913  1 41.8457 5.723e-08 ***
# TZP_Average.RNAConc..ng.uL.          0.09690  1  6.0596  0.017641 *  
# Dissecton.Group                      0.55173 12  2.8753  0.004917 ** 
#   Residuals                            0.73556 46                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



pdf("Boxplot_TTPConc_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TTP_Average.RNAConc..ng.uL.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TTP_Average.RNAConc..ng.uL.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Concentration (ng/uL)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TTP_Average.RNAConc..ng.uL.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -229.329  -31.238    2.975   42.892  179.850 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         504.62      32.69  15.437  < 2e-16 ***
# Dissecton.Group2    -18.69      51.69  -0.362   0.7190    
# Dissecton.Group3      4.40      46.23   0.095   0.9245    
# Dissecton.Group4    -41.68      46.23  -0.902   0.3712    
# Dissecton.Group5    114.02      46.23   2.466   0.0168 *  
# Dissecton.Group6    -16.92      51.69  -0.327   0.7447    
# Dissecton.Group7     35.40      46.23   0.766   0.4471    
# Dissecton.Group8     94.56      48.48   1.950   0.0562 .  
# Dissecton.Group9    249.24      48.48   5.141 3.75e-06 ***
# Dissecton.Group10   340.73      46.23   7.371 9.32e-10 ***
# Dissecton.Group11    77.08      46.23   1.667   0.1011    
# Dissecton.Group12    23.11      44.55   0.519   0.6060    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 80.07 on 55 degrees of freedom
# Multiple R-squared:  0.7016,	Adjusted R-squared:  0.6419 
# F-statistic: 11.75 on 11 and 55 DF,  p-value: 6.893e-11

summary.lm(lm(TTP_Average.RNAConc..ng.uL.~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -196.292  -43.842    5.126   49.486  151.047 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       1168.7901   490.7820   2.381   0.0212 *  
# GenderF            -69.0343    39.5614  -1.745   0.0873 .  
# Age                 -0.3915     0.7630  -0.513   0.6102    
# Hours.Final         -1.2909     1.5753  -0.819   0.4165    
# pH                 -89.7916    70.8723  -1.267   0.2112    
# DiagnosisBP         30.0029    27.1766   1.104   0.2750    
# DiagnosisSchiz      -0.4577    24.2755  -0.019   0.9850    
# Dissecton.Group2   -29.2829    53.4275  -0.548   0.5861    
# Dissecton.Group3   -13.8679    47.5033  -0.292   0.7716    
# Dissecton.Group4   -17.4666    48.9943  -0.357   0.7230    
# Dissecton.Group5    98.7964    47.6899   2.072   0.0436 *  
# Dissecton.Group6    -7.3201    53.1416  -0.138   0.8910    
# Dissecton.Group7    29.4202    46.8750   0.628   0.5332    
# Dissecton.Group8    78.7084    50.3650   1.563   0.1245    
# Dissecton.Group9   229.8294    50.1201   4.586 3.15e-05 ***
# Dissecton.Group10  346.7522    46.7966   7.410 1.54e-09 ***
# Dissecton.Group11   62.6316    47.5834   1.316   0.1942    
# Dissecton.Group12   12.5114    45.0404   0.278   0.7823    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 80.26 on 49 degrees of freedom
# Multiple R-squared:  0.7329,	Adjusted R-squared:  0.6402 
# F-statistic: 7.907 on 17 and 49 DF,  p-value: 5.951e-09

summary.lm(lm(TTP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -192.090  -42.768    0.692   50.397  158.651 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       1168.5589   494.6942   2.362   0.0223 *  
#   Block.Weight..g.    31.1254    65.1727   0.478   0.6351    
# GenderF            -75.3759    42.0295  -1.793   0.0792 .  
# Age                 -0.4112     0.7702  -0.534   0.5959    
# Hours.Final         -1.2099     1.5969  -0.758   0.4524    
# pH                 -92.8521    71.7240  -1.295   0.2017    
# DiagnosisBP         29.4031    27.4220   1.072   0.2890    
# DiagnosisSchiz      -2.4366    24.8174  -0.098   0.9222    
# Dissecton.Group2   -33.4554    54.5575  -0.613   0.5426    
# Dissecton.Group3   -11.8902    48.0607  -0.247   0.8057    
# Dissecton.Group4   -19.4799    49.5645  -0.393   0.6960    
# Dissecton.Group5    97.2731    48.1757   2.019   0.0491 *  
#   Dissecton.Group6    -7.8889    53.5784  -0.147   0.8836    
# Dissecton.Group7    30.5226    47.3049   0.645   0.5219    
# Dissecton.Group8    74.8334    51.4108   1.456   0.1520    
# Dissecton.Group9   226.1244    51.1118   4.424 5.55e-05 ***
#   Dissecton.Group10  346.4703    47.1733   7.345 2.18e-09 ***
#   Dissecton.Group11   61.3765    48.0347   1.278   0.2075    
# Dissecton.Group12   12.8759    45.4059   0.284   0.7780    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 80.9 on 48 degrees of freedom
# Multiple R-squared:  0.7341,	Adjusted R-squared:  0.6344 
# F-statistic: 7.363 on 18 and 48 DF,  p-value: 1.559e-08

car::Anova(lm(TTP_Average.RNAConc..ng.uL.~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TTP_Average.RNAConc..ng.uL.
# Sum Sq Df F value  Pr(>F)    
# (Intercept)      36533  1  5.6715 0.02117 *  
# Gender           19614  1  3.0450 0.08726 .  
# Age               1696  1  0.2633 0.61016    
# Hours.Final       4325  1  0.6715 0.41650    
# pH               10340  1  1.6052 0.21116    
# Diagnosis         9282  2  0.7205 0.49160    
# Dissecton.Group 767930 12  9.9347 1.9e-09 ***
# Residuals       315633 49                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TTP260280_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TTP_Average.260.280~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Purity: Average 260/280", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TTP_Average.260.280~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Purity: Average 260/280", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TTP_Average.260.280~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))    
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.100833 -0.025500  0.001667  0.018750  0.104167 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.920000   0.020192  95.089   <2e-16 ***
#   Dissecton.Group2  -0.030000   0.031926  -0.940   0.3515    
# Dissecton.Group3   0.060833   0.028555   2.130   0.0376 *  
#   Dissecton.Group4  -0.014167   0.028555  -0.496   0.6218    
# Dissecton.Group5   0.059167   0.028555   2.072   0.0430 *  
#   Dissecton.Group6  -0.015000   0.031926  -0.470   0.6403    
# Dissecton.Group7   0.003333   0.028555   0.117   0.9075    
# Dissecton.Group8   0.071000   0.029949   2.371   0.0213 *  
#   Dissecton.Group9  -0.039000   0.029949  -1.302   0.1983    
# Dissecton.Group10 -0.060833   0.028555  -2.130   0.0376 *  
#   Dissecton.Group11  0.033333   0.028555   1.167   0.2481    
# Dissecton.Group12  0.025000   0.027516   0.909   0.3676    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04946 on 55 degrees of freedom
# Multiple R-squared:  0.4472,	Adjusted R-squared:  0.3366 
# F-statistic: 4.044 on 11 and 55 DF,  p-value: 0.0002398

summary.lm(lm(TTP_Average.260.280~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))   

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.11202 -0.02608 -0.00180  0.02099  0.11202 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.8800568  0.3043843   6.177 1.25e-07 ***
#   GenderF            0.0454929  0.0245361   1.854   0.0697 .  
# Age                0.0001653  0.0004732   0.349   0.7284    
# Hours.Final        0.0001283  0.0009770   0.131   0.8961    
# pH                 0.0045588  0.0439552   0.104   0.9178    
# DiagnosisBP       -0.0195947  0.0168550  -1.163   0.2506    
# DiagnosisSchiz    -0.0173856  0.0150557  -1.155   0.2538    
# Dissecton.Group2  -0.0294894  0.0331359  -0.890   0.3778    
# Dissecton.Group3   0.0710702  0.0294617   2.412   0.0196 *  
#   Dissecton.Group4  -0.0268793  0.0303864  -0.885   0.3807    
# Dissecton.Group5   0.0694711  0.0295774   2.349   0.0229 *  
#   Dissecton.Group6  -0.0099546  0.0329586  -0.302   0.7639    
# Dissecton.Group7   0.0101109  0.0290720   0.348   0.7295    
# Dissecton.Group8   0.0806779  0.0312365   2.583   0.0128 *  
#   Dissecton.Group9  -0.0292703  0.0310846  -0.942   0.3510    
# Dissecton.Group10 -0.0648609  0.0290234  -2.235   0.0300 *  
#   Dissecton.Group11  0.0435913  0.0295114   1.477   0.1460    
# Dissecton.Group12  0.0277598  0.0279342   0.994   0.3252    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04978 on 49 degrees of freedom
# Multiple R-squared:  0.5011,	Adjusted R-squared:  0.328 
# F-statistic: 2.895 on 17 and 49 DF,  p-value: 0.00188


#I added block weight, because I thought it might be related to Trizol contamination:
summary.lm(lm(TTP_Average.260.280~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))    

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.8796607  0.3019393   6.225 1.13e-07 ***
#   Block.Weight..g.   0.0533211  0.0397785   1.340   0.1864    
# GenderF            0.0346290  0.0256529   1.350   0.1834    
# Age                0.0001315  0.0004701   0.280   0.7808    
# Hours.Final        0.0002670  0.0009747   0.274   0.7853    
# pH                -0.0006840  0.0437772  -0.016   0.9876    
# DiagnosisBP       -0.0206223  0.0167371  -1.232   0.2239    
# DiagnosisSchiz    -0.0207757  0.0151474  -1.372   0.1766    
# Dissecton.Group2  -0.0366372  0.0332994  -1.100   0.2767    
# Dissecton.Group3   0.0744582  0.0293341   2.538   0.0144 *  
#   Dissecton.Group4  -0.0303284  0.0302519  -1.003   0.3211    
# Dissecton.Group5   0.0668616  0.0294043   2.274   0.0275 *  
#   Dissecton.Group6  -0.0109289  0.0327019  -0.334   0.7397    
# Dissecton.Group7   0.0119994  0.0288728   0.416   0.6796    
# Dissecton.Group8   0.0740395  0.0313789   2.360   0.0224 *  
#   Dissecton.Group9  -0.0356174  0.0311964  -1.142   0.2592    
# Dissecton.Group10 -0.0653439  0.0287925  -2.269   0.0278 *  
#   Dissecton.Group11  0.0414410  0.0293182   1.413   0.1640    
# Dissecton.Group12  0.0283843  0.0277137   1.024   0.3109    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04938 on 48 degrees of freedom
# Multiple R-squared:  0.5191,	Adjusted R-squared:  0.3388 
# F-statistic: 2.879 on 18 and 48 DF,  p-value: 0.0018


car::Anova(lm(TTP_Average.260.280~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TTP_Average.260.280
# Sum Sq Df F value    Pr(>F)    
# (Intercept)     0.094526  1 38.1503  1.25e-07 ***
# Gender          0.008518  1  3.4378 0.0697459 .  
# Age             0.000302  1  0.1220 0.7283892    
# Hours.Final     0.000043  1  0.0172 0.8960978    
# pH              0.000027  1  0.0108 0.9178195    
# Diagnosis       0.004827  2  0.9741 0.3847099    
# Dissecton.Group 0.116303 12  3.9116 0.0003213 ***
# Residuals       0.121409 49                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TTP260230_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TTP_Average.260.230~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Purity: Average 260/230", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TTP_Average.260.230~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Purity: Average 260/230", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TTP_Average.260.230~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))  

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.68500 -0.12508 -0.00333  0.11354  0.57750 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.97333    0.09976  19.780  < 2e-16 ***
#   Dissecton.Group2  -0.03833    0.15774  -0.243 0.808896    
# Dissecton.Group3  -0.07417    0.14109  -0.526 0.601222    
# Dissecton.Group4   0.06833    0.14109   0.484 0.630070    
# Dissecton.Group5   0.13500    0.14109   0.957 0.342825    
# Dissecton.Group6  -0.02708    0.15774  -0.172 0.864306    
# Dissecton.Group7  -0.52083    0.14109  -3.692 0.000513 ***
#   Dissecton.Group8  -0.29233    0.14797  -1.976 0.053227 .  
# Dissecton.Group9  -1.38433    0.14797  -9.355 5.74e-13 ***
#   Dissecton.Group10 -1.41167    0.14109 -10.006 5.47e-14 ***
#   Dissecton.Group11 -0.35583    0.14109  -2.522 0.014590 *  
#   Dissecton.Group12 -0.24333    0.13595  -1.790 0.078989 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2444 on 55 degrees of freedom
# Multiple R-squared:  0.8368,	Adjusted R-squared:  0.8042 
# F-statistic: 25.64 on 11 and 55 DF,  p-value: < 2.2e-16

#Almost *all* variation in 260/230 is predicted by dissection group.


summary.lm(lm(TTP_Average.260.230~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))   

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.6048 -0.1226 -0.0187  0.1356  0.5677 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        3.3308102  1.5600426   2.135   0.0378 *  
#   GenderF            0.0360952  0.1257534   0.287   0.7753    
# Age               -0.0002986  0.0024254  -0.123   0.9025    
# Hours.Final        0.0002808  0.0050074   0.056   0.9555    
# pH                -0.2010999  0.2252808  -0.893   0.3764    
# DiagnosisBP        0.0324725  0.0863858   0.376   0.7086    
# DiagnosisSchiz     0.0397889  0.0771643   0.516   0.6084    
# Dissecton.Group2  -0.0757520  0.1698294  -0.446   0.6575    
# Dissecton.Group3  -0.0813170  0.1509981  -0.539   0.5926    
# Dissecton.Group4   0.0485805  0.1557376   0.312   0.7564    
# Dissecton.Group5   0.1233333  0.1515912   0.814   0.4198    
# Dissecton.Group6  -0.0086252  0.1689206  -0.051   0.9595    
# Dissecton.Group7  -0.5086169  0.1490008  -3.414   0.0013 ** 
#   Dissecton.Group8  -0.2758741  0.1600946  -1.723   0.0912 .  
# Dissecton.Group9  -1.3769419  0.1593163  -8.643 2.02e-11 ***
#   Dissecton.Group10 -1.4209905  0.1487519  -9.553 8.99e-13 ***
#   Dissecton.Group11 -0.3650292  0.1512528  -2.413   0.0196 *  
#   Dissecton.Group12 -0.2592776  0.1431695  -1.811   0.0763 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2551 on 49 degrees of freedom
# Multiple R-squared:  0.8415,	Adjusted R-squared:  0.7866 
# F-statistic: 15.31 on 17 and 49 DF,  p-value: 4.406e-14

#I added block weight because I thought it might be related to Trizol contamination:
summary.lm(lm(TTP_Average.260.230~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))   

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.58715 -0.11240 -0.02059  0.11368  0.59204 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        3.3298397  1.5696970   2.121  0.03908 *  
#   Block.Weight..g.   0.1306520  0.2067974   0.632  0.53052    
# GenderF            0.0094755  0.1333623   0.071  0.94365    
# Age               -0.0003813  0.0024439  -0.156  0.87666    
# Hours.Final        0.0006209  0.0050671   0.123  0.90299    
# pH                -0.2139463  0.2275850  -0.940  0.35189    
# DiagnosisBP        0.0299548  0.0870117   0.344  0.73215    
# DiagnosisSchiz     0.0314821  0.0787472   0.400  0.69109    
# Dissecton.Group2  -0.0932664  0.1731144  -0.539  0.59255    
# Dissecton.Group3  -0.0730155  0.1524996  -0.479  0.63426    
# Dissecton.Group4   0.0401293  0.1572712   0.255  0.79969    
# Dissecton.Group5   0.1169393  0.1528646   0.765  0.44802    
# Dissecton.Group6  -0.0110126  0.1700079  -0.065  0.94862    
# Dissecton.Group7  -0.5039896  0.1501017  -3.358  0.00155 ** 
#   Dissecton.Group8  -0.2921401  0.1631297  -1.791  0.07963 .  
# Dissecton.Group9  -1.3924942  0.1621812  -8.586 2.91e-11 ***
#   Dissecton.Group10 -1.4221739  0.1496841  -9.501 1.32e-12 ***
#   Dissecton.Group11 -0.3702979  0.1524171  -2.430  0.01891 *  
#   Dissecton.Group12 -0.2577474  0.1440758  -1.789  0.07993 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2567 on 48 degrees of freedom
# Multiple R-squared:  0.8428,	Adjusted R-squared:  0.7839 
# F-statistic:  14.3 on 18 and 48 DF,  p-value: 1.466e-13


car::Anova(lm(TTP_Average.260.230~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TTP_Average.260.230
# Sum Sq Df F value    Pr(>F)    
# (Intercept)      0.2967  1  4.5586   0.03778 *  
# Gender           0.0054  1  0.0824   0.77530    
# Age              0.0010  1  0.0152   0.90251    
# Hours.Final      0.0002  1  0.0031   0.95551    
# pH               0.0519  1  0.7968   0.37640    
# Diagnosis        0.0198  2  0.1520   0.85936    
# Dissecton.Group 15.5406 12 19.8978 7.762e-15 ***
# Residuals        3.1892 49                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TTP_RIN_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TTP_BioAnalyzer_RIN~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TTP_BioAnalyzer_RIN~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TTP RNA Integrity Number (RIN)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TTP_BioAnalyzer_RIN~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))  

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8714 -0.2000  0.0500  0.2500  0.9286 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.750e+00  2.040e-01  37.989  < 2e-16 ***
#   Dissecton.Group2  -2.750e-01  3.226e-01  -0.853   0.3976    
# Dissecton.Group3  -2.667e-01  2.885e-01  -0.924   0.3594    
# Dissecton.Group4   2.000e-01  2.885e-01   0.693   0.4911    
# Dissecton.Group5  -5.333e-01  2.885e-01  -1.849   0.0699 .  
# Dissecton.Group6  -1.625e+00  3.226e-01  -5.038 5.41e-06 ***
#   Dissecton.Group7   5.743e-15  2.885e-01   0.000   1.0000    
# Dissecton.Group8   9.000e-02  3.026e-01   0.297   0.7673    
# Dissecton.Group9  -2.700e-01  3.026e-01  -0.892   0.3761    
# Dissecton.Group10  6.083e-15  2.885e-01   0.000   1.0000    
# Dissecton.Group11 -7.667e-01  2.885e-01  -2.657   0.0103 *  
#   Dissecton.Group12 -4.786e-01  2.780e-01  -1.721   0.0908 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4997 on 55 degrees of freedom
# Multiple R-squared:  0.4819,	Adjusted R-squared:  0.3782 
# F-statistic:  4.65 on 11 and 55 DF,  p-value: 5.492e-05

summary.lm(lm(TTP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.47522 -0.20998  0.05379  0.27395  1.09424 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.106100   2.821083   0.747   0.4589    
# GenderF            0.247547   0.227405   1.089   0.2817    
# Age                0.006101   0.004386   1.391   0.1705    
# Hours.Final       -0.022494   0.009055  -2.484   0.0165 *  
# pH                 0.858926   0.407384   2.108   0.0401 *  
# DiagnosisBP       -0.072653   0.156215  -0.465   0.6439    
# DiagnosisSchiz    -0.263868   0.139539  -1.891   0.0645 .  
# Dissecton.Group2  -0.147832   0.307109  -0.481   0.6324    
# Dissecton.Group3  -0.195245   0.273055  -0.715   0.4780    
# Dissecton.Group4   0.236657   0.281626   0.840   0.4048    
# Dissecton.Group5  -0.337454   0.274128  -1.231   0.2242    
# Dissecton.Group6  -1.520483   0.305465  -4.978 8.35e-06 ***
# Dissecton.Group7   0.030112   0.269444   0.112   0.9115    
# Dissecton.Group8   0.004273   0.289505   0.015   0.9883    
# Dissecton.Group9  -0.316598   0.288098  -1.099   0.2772    
# Dissecton.Group10 -0.046515   0.268993  -0.173   0.8634    
# Dissecton.Group11 -0.587824   0.273516  -2.149   0.0366 *  
# Dissecton.Group12 -0.397205   0.258899  -1.534   0.1314    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4613 on 49 degrees of freedom
# Multiple R-squared:  0.6066,	Adjusted R-squared:  0.4701 
# F-statistic: 4.444 on 17 and 49 DF,  p-value: 2.078e-05



car::Anova(lm(TTP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingPH)==F,]), type=3)
# Anova Table (Type III tests)
# 
# Response: TTP_BioAnalyzer_RIN
# Sum Sq Df F value    Pr(>F)    
# (Intercept)      0.1186  1  0.5573 0.4588979    
# Gender           0.2522  1  1.1850 0.2816658    
# Age              0.4118  1  1.9348 0.1705222    
# Hours.Final      1.3134  1  6.1709 0.0164533 *  
# pH               0.9461  1  4.4453 0.0401335 *  
# Diagnosis        0.7766  2  1.8244 0.1721033    
# Dissecton.Group  9.6348 12  3.7724 0.0004538 ***
#   Residuals       10.4289 49                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Ooh, that's pretty. :)

pdf("Boxplot_TZPConc_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TZP_Average.RNAConc..ng.uL.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TZP_Average.RNAConc..ng.uL.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Concentration (ng/uL)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))    

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -201.12  -25.06    4.61   33.23  140.69 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       417.9533    24.6912  16.927   <2e-16 ***
#   Dissecton.Group2  -52.3421    39.0402  -1.341   0.1855    
# Dissecton.Group3  -27.5117    34.9186  -0.788   0.4341    
# Dissecton.Group4  -49.1267    34.9186  -1.407   0.1651    
# Dissecton.Group5   29.2492    34.9186   0.838   0.4059    
# Dissecton.Group6  -70.0408    39.0402  -1.794   0.0783 .  
# Dissecton.Group7  -52.2783    34.9186  -1.497   0.1401    
# Dissecton.Group8   -0.5803    36.6230  -0.016   0.9874    
# Dissecton.Group9   32.4337    36.6230   0.886   0.3797    
# Dissecton.Group10  25.2300    34.9186   0.723   0.4730    
# Dissecton.Group11   1.2792    34.9186   0.037   0.9709    
# Dissecton.Group12 -37.9433    33.6484  -1.128   0.2644    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 60.48 on 55 degrees of freedom
# Multiple R-squared:  0.2698,	Adjusted R-squared:  0.1237 
# F-statistic: 1.847 on 11 and 55 DF,  p-value: 0.06781

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))  

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -194.961  -22.483    5.044   32.990  126.766 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)       564.18056  382.92344   1.473    0.147
# GenderF           -21.55825   30.86705  -0.698    0.488
# Age                -0.05351    0.59532  -0.090    0.929
# Hours.Final        -1.29371    1.22911  -1.053    0.298
# pH                -16.13176   55.29676  -0.292    0.772
# DiagnosisBP         8.49676   21.20401   0.401    0.690
# DiagnosisSchiz    -16.41587   18.94052  -0.867    0.390
# Dissecton.Group2  -52.93246   41.68582  -1.270    0.210
# Dissecton.Group3  -32.54423   37.06355  -0.878    0.384
# Dissecton.Group4  -37.09684   38.22690  -0.970    0.337
# Dissecton.Group5   29.26204   37.20912   0.786    0.435
# Dissecton.Group6  -60.59907   41.46275  -1.462    0.150
# Dissecton.Group7  -53.44257   36.57330  -1.461    0.150
# Dissecton.Group8  -10.31860   39.29634  -0.263    0.794
# Dissecton.Group9   22.25489   39.10531   0.569    0.572
# Dissecton.Group10  26.56613   36.51219   0.728    0.470
# Dissecton.Group11   1.14417   37.12607   0.031    0.976
# Dissecton.Group12 -40.33924   35.14195  -1.148    0.257
# 
# Residual standard error: 62.62 on 49 degrees of freedom
# Multiple R-squared:  0.3026,	Adjusted R-squared:  0.06061 
# F-statistic:  1.25 on 17 and 49 DF,  p-value: 0.264

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -180.792  -20.185   -6.594   35.939  148.753 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)       563.4010   369.4089   1.525   0.1338  
# Block.Weight..g.  104.9553    48.6672   2.157   0.0361 *
# GenderF           -42.9424    31.3852  -1.368   0.1776  
# Age                -0.1200     0.5751  -0.209   0.8357  
# Hours.Final        -1.0205     1.1925  -0.856   0.3964  
# pH                -26.4515    53.5593  -0.494   0.6236  
# DiagnosisBP         6.4742    20.4771   0.316   0.7532  
# DiagnosisSchiz    -23.0889    18.5322  -1.246   0.2189  
# Dissecton.Group2  -67.0021    40.7403  -1.645   0.1066  
# Dissecton.Group3  -25.8754    35.8889  -0.721   0.4744  
# Dissecton.Group4  -43.8858    37.0119  -1.186   0.2416  
# Dissecton.Group5   24.1256    35.9748   0.671   0.5057  
# Dissecton.Group6  -62.5169    40.0093  -1.563   0.1247  
# Dissecton.Group7  -49.7254    35.3246  -1.408   0.1657  
# Dissecton.Group8  -23.3854    38.3906  -0.609   0.5453  
# Dissecton.Group9    9.7614    38.1674   0.256   0.7992  
# Dissecton.Group10  25.6155    35.2263   0.727   0.4707  
# Dissecton.Group11  -3.0882    35.8695  -0.086   0.9317  
# Dissecton.Group12 -39.1101    33.9065  -1.153   0.2544  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 60.41 on 48 degrees of freedom
# Multiple R-squared:  0.3642,	Adjusted R-squared:  0.1257 
# F-statistic: 1.527 on 18 and 48 DF,  p-value: 0.1218

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value  Pr(>F)  
# (Intercept)        8489  1  2.3261 0.13379  
# Block.Weight..g.  16973  1  4.6509 0.03608 *
# Gender             6832  1  1.8721 0.17761  
# Age                 159  1  0.0435 0.83567  
# Hours.Final        2673  1  0.7324 0.39635  
# pH                  890  1  0.2439 0.62365  
# Diagnosis          8325  2  1.1406 0.32814  
# Dissecton.Group   58353 11  1.4536 0.18080  
# Residuals        175172 48                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s. +Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value  Pr(>F)  
# (Intercept)                            2546  1  0.6975 0.40794  
# Block.Weight..g.                      22123  1  6.0596 0.01764 *
# Gender                                10165  1  2.7842 0.10199  
# Age                                    2124  1  0.5817 0.44952  
# Hours.Final                             127  1  0.0348 0.85293  
# pH                                      141  1  0.0387 0.84491  
# Diagnosis                              8140  2  1.1149 0.33666  
# TZP_BioAnalyzer_RIN                    5305  1  1.4532 0.23418  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   6214  1  1.7020 0.19851  
# Dissecton.Group                       60842 11  1.5150 0.15880  
# Residuals                            167940 46                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value    Pr(>F)    
# (Intercept)                           45726  1 16.2748 0.0002150 ***
# Block.Weight..g.                       4969  1  1.7685 0.1904242    
# Gender                                 4940  1  1.7584 0.1916669    
# Age                                     235  1  0.0837 0.7737598    
# Hours.Final                            2849  1  1.0142 0.3194147    
# pH                                     2992  1  1.0648 0.3077580    
# Diagnosis                             11224  2  1.9974 0.1477978    
# TZP_BioAnalyzer_RIN                    1260  1  0.4483 0.5066435    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.    274  1  0.0975 0.7563148    
# TZP_Average.260.280                   43720  1 15.5608 0.0002834 ***
# TZP_Average.260.230                    3191  1  1.1358 0.2923652    
# Dissecton.Group                       36120 11  1.1687 0.3356493    
# Residuals                            123623 44                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Gender+Age+Hours.Final+pH+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           3780.9284   937.2180   4.034 0.000215 ***
#   Block.Weight..g.                        82.3634    61.9349   1.330 0.190424    
# GenderF                                -38.3069    28.8880  -1.326 0.191667    
# Age                                     -0.1630     0.5635  -0.289 0.773760    
# Hours.Final                             -1.1907     1.1824  -1.007 0.319415    
# pH                                     -53.3181    51.6698  -1.032 0.307758    
# DiagnosisBP                             -6.9929    18.2693  -0.383 0.703735    
# DiagnosisSchiz                         -32.7282    16.7594  -1.953 0.057219 .  
# TZP_BioAnalyzer_RIN                     13.4370    20.0688   0.670 0.506643    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   -11.7409    37.5994  -0.312 0.756315    
# TZP_Average.260.280                  -1483.5088   376.0744  -3.945 0.000283 ***
# TZP_Average.260.230                    -47.1247    44.2188  -1.066 0.292365    
# Dissecton.Group2                       -27.7158    40.8927  -0.678 0.501468    
# Dissecton.Group3                       -13.4469    34.8826  -0.385 0.701731    
# Dissecton.Group4                        -0.2819    37.4020  -0.008 0.994021    
# Dissecton.Group5                        42.8697    32.3157   1.327 0.191489    
# Dissecton.Group6                        -7.3788    42.5895  -0.173 0.863247    
# Dissecton.Group7                       -24.7761    32.5972  -0.760 0.451266    
# Dissecton.Group8                       -22.2718    37.1401  -0.600 0.551801    
# Dissecton.Group9                        11.9438    33.8251   0.353 0.725695    
# Dissecton.Group10                       50.4089    32.0457   1.573 0.122875    
# Dissecton.Group11                        9.1367    32.8864   0.278 0.782448    
# Dissecton.Group12                      -15.1686    30.6580  -0.495 0.623225    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 53.01 on 44 degrees of freedom
# Multiple R-squared:  0.5513,	Adjusted R-squared:  0.3269 
# F-statistic: 2.457 on 22 and 44 DF,  p-value: 0.005521

#There is a trend towards RNA conc being lower in Schiz.... in an extremely unwieldy model that may or may not make sense.

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Diagnosis+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))


pdf("Boxplot_TZP260280_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TZP_Average.260.280~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Purity: Average 260/280", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TZP_Average.260.280~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Purity: Average 260/280", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_Average.260.280~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))    

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.04333 -0.01333 -0.00250  0.01000  0.04500 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.013333   0.008767 229.646   <2e-16 ***
#   Dissecton.Group2   0.019167   0.013862   1.383   0.1724    
# Dissecton.Group3   0.011667   0.012399   0.941   0.3508    
# Dissecton.Group4   0.028333   0.012399   2.285   0.0262 *  
# Dissecton.Group5   0.008333   0.012399   0.672   0.5043    
# Dissecton.Group6   0.026667   0.013862   1.924   0.0596 .  
# Dissecton.Group7   0.020000   0.012399   1.613   0.1124    
# Dissecton.Group8  -0.001333   0.013004  -0.103   0.9187    
# Dissecton.Group9  -0.003333   0.013004  -0.256   0.7986    
# Dissecton.Group10  0.013333   0.012399   1.075   0.2869    
# Dissecton.Group11  0.001667   0.012399   0.134   0.8936    
# Dissecton.Group12  0.013810   0.011948   1.156   0.2527    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02147 on 55 degrees of freedom
# Multiple R-squared:  0.209,	Adjusted R-squared:  0.05078 
# F-statistic: 1.321 on 11 and 55 DF,  p-value: 0.2381

summary.lm(lm(TZP_Average.260.280~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.038669 -0.012423 -0.000440  0.008134  0.044774 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.159e+00  1.335e-01  16.165   <2e-16 ***
# GenderF            2.819e-03  1.076e-02   0.262   0.7945    
# Age                6.966e-05  2.076e-04   0.336   0.7387    
# Hours.Final       -1.697e-04  4.286e-04  -0.396   0.6939    
# pH                -2.064e-02  1.928e-02  -1.070   0.2897    
# DiagnosisBP       -9.518e-03  7.395e-03  -1.287   0.2041    
# DiagnosisSchiz    -9.151e-03  6.605e-03  -1.385   0.1722    
# Dissecton.Group2   1.852e-02  1.454e-02   1.274   0.2088    
# Dissecton.Group3   1.253e-02  1.293e-02   0.969   0.3372    
# Dissecton.Group4   2.958e-02  1.333e-02   2.218   0.0312 *  
# Dissecton.Group5   9.896e-03  1.298e-02   0.763   0.4493    
# Dissecton.Group6   2.894e-02  1.446e-02   2.001   0.0509 .  
# Dissecton.Group7   2.108e-02  1.275e-02   1.653   0.1047    
# Dissecton.Group8   8.720e-04  1.370e-02   0.064   0.9495    
# Dissecton.Group9  -2.099e-03  1.364e-02  -0.154   0.8783    
# Dissecton.Group10  1.467e-02  1.273e-02   1.152   0.2548    
# Dissecton.Group11  3.413e-03  1.295e-02   0.264   0.7932    
# Dissecton.Group12  1.318e-02  1.226e-02   1.076   0.2874    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02184 on 49 degrees of freedom
# Multiple R-squared:  0.2712,	Adjusted R-squared:  0.01837 
# F-statistic: 1.073 on 17 and 49 DF,  p-value: 0.4051


car::Anova(lm(TZP_Average.260.280~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.260.280
# Sum Sq Df  F value Pr(>F)    
# (Intercept)     0.124616  1 261.2933 <2e-16 ***
# Gender          0.000033  1   0.0686 0.7945    
# Age             0.000054  1   0.1126 0.7387    
# Hours.Final     0.000075  1   0.1567 0.6939    
# pH              0.000546  1   1.1456 0.2897    
# Diagnosis       0.001240  2   1.2996 0.2819    
# Dissecton.Group 0.005556 11   1.0591 0.4124    
# Residuals       0.023369 49                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP260230_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TZP_Average.260.230~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Purity: Average 260/230", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TZP_Average.260.230~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Purity: Average 260/230", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_Average.260.230~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))  

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49200 -0.09583  0.02000  0.10857  0.38167 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.106667   0.068707  30.662   <2e-16 ***
# Dissecton.Group2   0.173333   0.108635   1.596   0.1163    
# Dissecton.Group3  -0.188333   0.097166  -1.938   0.0577 .  
# Dissecton.Group4   0.048333   0.097166   0.497   0.6209    
# Dissecton.Group5  -0.040000   0.097166  -0.412   0.6822    
# Dissecton.Group6   0.010833   0.108635   0.100   0.9209    
# Dissecton.Group7  -0.090000   0.097166  -0.926   0.3584    
# Dissecton.Group8  -0.186667   0.101909  -1.832   0.0724 .  
# Dissecton.Group9   0.005333   0.101909   0.052   0.9585    
# Dissecton.Group10  0.006667   0.097166   0.069   0.9455    
# Dissecton.Group11  0.028333   0.097166   0.292   0.7717    
# Dissecton.Group12  0.016190   0.093632   0.173   0.8634    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1683 on 55 degrees of freedom
# Multiple R-squared:  0.2626,	Adjusted R-squared:  0.1152 
# F-statistic: 1.781 on 11 and 55 DF,  p-value: 0.08001


summary.lm(lm(TZP_Average.260.230~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))   

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49325 -0.08872  0.01388  0.10889  0.34681 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        1.9322932  1.0731601   1.801   0.0779 .
# GenderF           -0.0518322  0.0865063  -0.599   0.5518  
# Age               -0.0007403  0.0016684  -0.444   0.6592  
# Hours.Final       -0.0028049  0.0034446  -0.814   0.4194  
# pH                 0.0395542  0.1549716   0.255   0.7996  
# DiagnosisBP       -0.0006422  0.0594252  -0.011   0.9914  
# DiagnosisSchiz     0.0106234  0.0530816   0.200   0.8422  
# Dissecton.Group2   0.1970090  0.1168264   1.686   0.0981 .
# Dissecton.Group3  -0.1906488  0.1038722  -1.835   0.0725 .
# Dissecton.Group4   0.0806189  0.1071326   0.753   0.4553  
# Dissecton.Group5  -0.0325878  0.1042802  -0.313   0.7560  
# Dissecton.Group6   0.0189346  0.1162012   0.163   0.8712  
# Dissecton.Group7  -0.0936118  0.1024983  -0.913   0.3656  
# Dissecton.Group8  -0.2016586  0.1101297  -1.831   0.0732 .
# Dissecton.Group9  -0.0066860  0.1095944  -0.061   0.9516  
# Dissecton.Group10  0.0117836  0.1023270   0.115   0.9088  
# Dissecton.Group11  0.0346595  0.1040475   0.333   0.7405  
# Dissecton.Group12  0.0210007  0.0984869   0.213   0.8320  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1755 on 49 degrees of freedom
# Multiple R-squared:  0.2857,	Adjusted R-squared:  0.03783 
# F-statistic: 1.153 on 17 and 49 DF,  p-value: 0.3363

car::Anova(lm(TZP_Average.260.230~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.260.230
# Sum Sq Df F value  Pr(>F)  
# (Intercept)     0.09985  1  3.2420 0.07793 .
# Gender          0.01106  1  0.3590 0.55182  
# Age             0.00606  1  0.1969 0.65922  
# Hours.Final     0.02042  1  0.6630 0.41943  
# pH              0.00201  1  0.0651 0.79961  
# Diagnosis       0.00150  2  0.0244 0.97592  
# Dissecton.Group 0.55791 11  1.6468 0.11505  
# Residuals       1.50916 49                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("Boxplot_TZP_RIN_ByDissectionGroup_noBad.pdf", width=12, height=6)
boxplot(TZP_BioAnalyzer_RIN~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(TZP_BioAnalyzer_RIN~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Dissection/Extraction Group", ylab="TZP RNA Integrity Number (RIN)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

summary.lm(lm(TZP_BioAnalyzer_RIN~Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))  

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.72857 -0.20167  0.06667  0.26250  0.87500 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        8.00000    0.20410  39.197  < 2e-16 ***
#   Dissecton.Group2  -0.27500    0.32271  -0.852 0.397817    
# Dissecton.Group3  -0.11667    0.28864  -0.404 0.687635    
# Dissecton.Group4   0.35000    0.28864   1.213 0.230469    
# Dissecton.Group5  -0.36667    0.28864  -1.270 0.209313    
# Dissecton.Group6  -1.35000    0.32271  -4.183 0.000104 ***
# Dissecton.Group7  -0.08333    0.28864  -0.289 0.773887    
# Dissecton.Group8   0.02000    0.30273   0.066 0.947565    
# Dissecton.Group9  -0.20000    0.30273  -0.661 0.511586    
# Dissecton.Group10  0.08333    0.28864   0.289 0.773887    
# Dissecton.Group11 -0.26667    0.28864  -0.924 0.359585    
# Dissecton.Group12 -0.27143    0.27814  -0.976 0.333399    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4999 on 55 degrees of freedom
# Multiple R-squared:  0.3793,	Adjusted R-squared:  0.2552 
# F-statistic: 3.055 on 11 and 55 DF,  p-value: 0.002961

summary.lm(lm(TZP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.42756 -0.26808  0.08467  0.24164  1.01001 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        3.471325   2.941516   1.180 0.243652    
# GenderF            0.250872   0.237113   1.058 0.295230    
# Age                0.007226   0.004573   1.580 0.120499    
# Hours.Final       -0.018037   0.009442  -1.910 0.061942 .  
# pH                 0.669853   0.424775   1.577 0.121240    
# DiagnosisBP       -0.033457   0.162884  -0.205 0.838106    
# DiagnosisSchiz    -0.180757   0.145496  -1.242 0.220019    
# Dissecton.Group2  -0.200505   0.320219  -0.626 0.534124    
# Dissecton.Group3  -0.066545   0.284712  -0.234 0.816171    
# Dissecton.Group4   0.354565   0.293649   1.207 0.233057    
# Dissecton.Group5  -0.209185   0.285831  -0.732 0.467746    
# Dissecton.Group6  -1.273176   0.318506  -3.997 0.000215 ***
# Dissecton.Group7  -0.056927   0.280946  -0.203 0.840266    
# Dissecton.Group8  -0.043233   0.301864  -0.143 0.886703    
# Dissecton.Group9  -0.234593   0.300397  -0.781 0.438591    
# Dissecton.Group10  0.036895   0.280477   0.132 0.895885    
# Dissecton.Group11 -0.123015   0.285193  -0.431 0.668114    
# Dissecton.Group12 -0.208580   0.269951  -0.773 0.443435    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.481 on 49 degrees of freedom
# Multiple R-squared:  0.488,	Adjusted R-squared:  0.3104 
# F-statistic: 2.748 on 17 and 49 DF,  p-value: 0.002963

car::Anova(lm(TZP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_BioAnalyzer_RIN
# Sum Sq Df F value  Pr(>F)  
# (Intercept)      0.3223  1  1.3927 0.24365  
# Gender           0.2590  1  1.1194 0.29523  
# Age              0.5778  1  2.4970 0.12050  
# Hours.Final      0.8445  1  3.6496 0.06194 .
# pH               0.5754  1  2.4868 0.12124  
# Diagnosis        0.3760  2  0.8125 0.44963  
# Dissecton.Group  6.5965 11  2.5916 0.01110 *
# Residuals       11.3383 49                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##So it seems like Dissection/Extraction Group is a particularly big deal at the TTP level, where we had some pretty profound contamination in 2 dissection groups that ended up strongly influencing RNA concentration values.

##By the time we are at the TZP level, there is still a little bit of evidence that one of the groups differed for contamination or RIN, but it is a much weaker effect and doesn't seem to be influencing concentration.

########

#What about cohort?  Cohort for the most part should be meaningless (because there is a secondary dissection), but it does partially reflect tissue age, so maybe I should look at it in respect to RIN? Unfortunately, it is also unbalanced based on diagnosis, so I'm not sure if I can interpret any results I get, but let's take a peek...

table(SubjectInfo$Cohort[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F], SubjectInfo$Dissecton.Group[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F])
#                1 2 3 4 5 6 7 8 9 10 11 12 6redone
# Cohort 11      1 0 1 0 1 1 1 1 1  0  1  1       0
# Cohort 12      1 1 1 1 2 1 1 2 0  1  2  1       0
# Cohort 13      1 0 0 1 0 0 1 1 2  0  0  1       0
# Cohort 7       0 1 1 2 2 0 1 0 0  2  0  1       0
# Cohort 8       0 1 1 0 0 0 1 0 1  1  0  1       0
# Dep Cohort 1   1 0 0 0 0 0 0 1 1  0  0  1       0
# Dep Cohort 5   0 1 1 1 0 0 0 0 0  1  1  0       0
# Dep Cohort 6   1 0 0 1 0 0 0 0 0  0  0  1       0
# Schiz Cohort 2 1 0 1 0 1 2 1 0 0  1  2  0       0

summary.lm(lm(TTP_BioAnalyzer_RIN~Cohort, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.10000 -0.26000  0.08571  0.38452  0.85000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           7.36667    0.20623  35.720   <2e-16 ***
# CohortCohort 12       0.08333    0.26434   0.315    0.754    
# CohortCohort 13       0.24762    0.31180   0.794    0.430    
# CohortCohort 7        0.19333    0.28427   0.680    0.499    
# CohortCohort 8       -0.55000    0.32609  -1.687    0.097 .  
# CohortDep Cohort 1    0.35833    0.37179   0.964    0.339    
# CohortDep Cohort 5    0.29333    0.34510   0.850    0.399    
# CohortDep Cohort 6    0.60000    0.41247   1.455    0.151    
# CohortSchiz Cohort 2 -0.06667    0.29166  -0.229    0.820    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6187 on 58 degrees of freedom
# Multiple R-squared:  0.1624,	Adjusted R-squared:  0.04689 
# F-statistic: 1.406 on 8 and 58 DF,  p-value: 0.2135

summary.lm(lm(TTP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.39777 -0.33016  0.08459  0.42797  0.79471 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           2.127712   3.864142   0.551  0.58425   
# GenderF               0.274189   0.363547   0.754  0.45413   
# Age                   0.005174   0.006001   0.862  0.39257   
# Hours.Final          -0.032578   0.010839  -3.006  0.00407 **
# pH                    0.833080   0.560942   1.485  0.14354   
# DiagnosisBP          -0.066100   0.217611  -0.304  0.76253   
# DiagnosisSchiz       -0.247242   0.199725  -1.238  0.22131   
# CohortCohort 12       0.162476   0.254012   0.640  0.52522   
# CohortCohort 13       0.409673   0.300467   1.363  0.17861   
# CohortCohort 7        0.357527   0.307122   1.164  0.24969   
# CohortCohort 8       -0.318543   0.331055  -0.962  0.34040   
# CohortDep Cohort 1    0.249836   0.387841   0.644  0.52230   
# CohortDep Cohort 5    0.251160   0.416046   0.604  0.54868   
# CohortDep Cohort 6    0.477585   0.547706   0.872  0.38723   
# CohortSchiz Cohort 2  0.039461   0.281285   0.140  0.88897   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5794 on 52 degrees of freedom
# Multiple R-squared:  0.3414,	Adjusted R-squared:  0.1641 
# F-statistic: 1.926 on 14 and 52 DF,  p-value: 0.04478

car::Anova(lm(TTP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TTP_BioAnalyzer_RIN
# Sum Sq Df F value   Pr(>F)   
# (Intercept)  0.1018  1  0.3032 0.584245   
# Gender       0.1910  1  0.5688 0.454129   
# Age          0.2495  1  0.7433 0.392574   
# Hours.Final  3.0327  1  9.0338 0.004073 **
# pH           0.7404  1  2.2057 0.143543   
# Diagnosis    0.5433  2  0.8092 0.450752   
# Cohort       2.4408  8  0.9088 0.516516   
# Residuals   17.4566 52                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#This is pushing it as far as df:
summary.lm(lm(TTP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.91739 -0.21893  0.05646  0.21422  0.88715 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           5.404476   3.118536   1.733  0.09061 .  
# GenderF              -0.091985   0.318519  -0.289  0.77420    
# Age                   0.004596   0.004551   1.010  0.31845    
# Hours.Final          -0.024519   0.009033  -2.714  0.00967 ** 
# pH                    0.372534   0.453999   0.821  0.41664    
# DiagnosisBP           0.022352   0.172982   0.129  0.89782    
# DiagnosisSchiz       -0.130187   0.152192  -0.855  0.39729    
# CohortCohort 12       0.025945   0.194280   0.134  0.89441    
# CohortCohort 13       0.052584   0.241251   0.218  0.82854    
# CohortCohort 7       -0.028413   0.249519  -0.114  0.90990    
# CohortCohort 8       -0.744399   0.267243  -2.785  0.00806 ** 
# CohortDep Cohort 1    0.112323   0.299705   0.375  0.70976    
# CohortDep Cohort 5    0.116252   0.326216   0.356  0.72339    
# CohortDep Cohort 6    0.397029   0.424691   0.935  0.35533    
# CohortSchiz Cohort 2  0.200485   0.216475   0.926  0.35980    
# Dissecton.Group2      0.096642   0.310097   0.312  0.75689    
# Dissecton.Group3     -0.078742   0.270496  -0.291  0.77244    
# Dissecton.Group4      0.384847   0.282889   1.360  0.18113    
# Dissecton.Group5     -0.327452   0.270713  -1.210  0.23336    
# Dissecton.Group6     -1.504883   0.300924  -5.001 1.12e-05 ***
# Dissecton.Group7      0.208347   0.263662   0.790  0.43396    
# Dissecton.Group8      0.056797   0.276251   0.206  0.83812    
# Dissecton.Group9     -0.136098   0.282922  -0.481  0.63304    
# Dissecton.Group10     0.203991   0.281831   0.724  0.47329    
# Dissecton.Group11    -0.634697   0.273094  -2.324  0.02515 *  
# Dissecton.Group12    -0.285795   0.245440  -1.164  0.25099    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4316 on 41 degrees of freedom
# Multiple R-squared:  0.7119,	Adjusted R-squared:  0.5362 
# F-statistic: 4.053 on 25 and 41 DF,  p-value: 3.671e-05

car::Anova(lm(TTP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TTP_BioAnalyzer_RIN
# Sum Sq Df F value    Pr(>F)    
# (Intercept)     0.5594  1  3.0033 0.0906091 .  
# Gender          0.0155  1  0.0834 0.7741976    
# Age             0.1900  1  1.0200 0.3184521    
# Hours.Final     1.3722  1  7.3674 0.0096715 ** 
# pH              0.1254  1  0.6733 0.4166418    
# Diagnosis       0.1956  2  0.5250 0.5954720    
# Cohort          2.7922  8  1.8739 0.0908026 .  
# Dissecton.Group 9.8199 11  4.7929 0.0001036 ***
#   Residuals       7.6366 41                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Just double-checking with a few less df:
summary.lm(lm(TTP_BioAnalyzer_RIN~Hours.Final+pH+Cohort+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 
#Looks basically the same

summary.lm(lm(TZP_BioAnalyzer_RIN~Cohort, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.7000 -0.2536  0.1000  0.3000  0.7667 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            7.6444     0.1831  41.751   <2e-16 ***
#   CohortCohort 12        0.2556     0.2347   1.089   0.2807    
# CohortCohort 13        0.3127     0.2768   1.130   0.2633    
# CohortCohort 7         0.2056     0.2524   0.814   0.4187    
# CohortCohort 8        -0.5111     0.2895  -1.765   0.0827 .  
# CohortDep Cohort 1     0.4306     0.3301   1.304   0.1973    
# CohortDep Cohort 5     0.4756     0.3064   1.552   0.1261    
# CohortDep Cohort 6     0.6222     0.3662   1.699   0.0946 .  
# CohortSchiz Cohort 2   0.1556     0.2589   0.601   0.5504    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5493 on 58 degrees of freedom
# Multiple R-squared:  0.2098,	Adjusted R-squared:  0.1008 
# F-statistic: 1.925 on 8 and 58 DF,  p-value: 0.07327

summary.lm(lm(TZP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.15706 -0.16197  0.07424  0.30471  0.82191 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           6.241539   3.555608   1.755   0.0851 .
# GenderF               0.331007   0.334520   0.989   0.3270  
# Age                   0.004939   0.005522   0.894   0.3752  
# Hours.Final          -0.024799   0.009974  -2.486   0.0162 *
# pH                    0.240471   0.516153   0.466   0.6432  
# DiagnosisBP           0.043503   0.200236   0.217   0.8289  
# DiagnosisSchiz       -0.038534   0.183778  -0.210   0.8347  
# CohortCohort 12       0.305836   0.233730   1.308   0.1965  
# CohortCohort 13       0.431930   0.276476   1.562   0.1243  
# CohortCohort 7        0.219419   0.282600   0.776   0.4410  
# CohortCohort 8       -0.455814   0.304622  -1.496   0.1406  
# CohortDep Cohort 1    0.443442   0.356873   1.243   0.2196  
# CohortDep Cohort 5    0.357440   0.382826   0.934   0.3548  
# CohortDep Cohort 6    0.299459   0.503974   0.594   0.5550  
# CohortSchiz Cohort 2  0.240173   0.258826   0.928   0.3577  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5331 on 52 degrees of freedom
# Multiple R-squared:  0.3326,	Adjusted R-squared:  0.1529 
# F-statistic: 1.851 on 14 and 52 DF,  p-value: 0.05529

car::Anova(lm(TZP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_BioAnalyzer_RIN
# Sum Sq Df F value  Pr(>F)  
# (Intercept)  0.8759  1  3.0815 0.08508 .
# Gender       0.2783  1  0.9791 0.32700  
# Age          0.2274  1  0.8000 0.37522  
# Hours.Final  1.7572  1  6.1823 0.01615 *
# pH           0.0617  1  0.2171 0.64324  
# Diagnosis    0.0476  2  0.0838 0.91975  
# Cohort       3.1546  8  1.3873 0.22406  
# Residuals   14.7802 52                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(TZP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])) 

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.05050 -0.19901  0.02366  0.20085  0.86237 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           8.125690   3.078438   2.640 0.011687 *  
#   GenderF               0.021881   0.314424   0.070 0.944858    
# Age                   0.005262   0.004492   1.171 0.248228    
# Hours.Final          -0.020308   0.008917  -2.277 0.028047 *  
# pH                   -0.031250   0.448162  -0.070 0.944747    
# DiagnosisBP           0.134125   0.170758   0.785 0.436694    
# DiagnosisSchiz        0.047789   0.150235   0.318 0.752027    
# CohortCohort 12       0.181990   0.191782   0.949 0.348210    
# CohortCohort 13       0.164992   0.238149   0.693 0.492335    
# CohortCohort 7       -0.129766   0.246310  -0.527 0.601143    
# CohortCohort 8       -0.797474   0.263807  -3.023 0.004303 ** 
# CohortDep Cohort 1    0.354486   0.295852   1.198 0.237724    
# CohortDep Cohort 5    0.161241   0.322021   0.501 0.619249    
# CohortDep Cohort 6    0.190976   0.419230   0.456 0.651124    
# CohortSchiz Cohort 2  0.346367   0.213692   1.621 0.112712    
# Dissecton.Group2      0.055632   0.306110   0.182 0.856684    
# Dissecton.Group3      0.091847   0.267018   0.344 0.732624    
# Dissecton.Group4      0.545150   0.279251   1.952 0.057768 .  
# Dissecton.Group5     -0.170766   0.267232  -0.639 0.526365    
# Dissecton.Group6     -1.261052   0.297055  -4.245 0.000122 ***
# Dissecton.Group7      0.179294   0.260272   0.689 0.494783    
# Dissecton.Group8     -0.022348   0.272699  -0.082 0.935083    
# Dissecton.Group9     -0.044642   0.279284  -0.160 0.873789    
# Dissecton.Group10     0.318401   0.278207   1.144 0.259065    
# Dissecton.Group11    -0.204551   0.269583  -0.759 0.452331    
# Dissecton.Group12    -0.060372   0.242284  -0.249 0.804468    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.426 on 41 degrees of freedom
# Multiple R-squared:  0.664,	Adjusted R-squared:  0.4591 
# F-statistic: 3.241 on 25 and 41 DF,  p-value: 0.0004119

car::Anova(lm(TZP_BioAnalyzer_RIN~Gender+Age+Hours.Final+pH+Diagnosis+Cohort+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_BioAnalyzer_RIN
# Sum Sq Df F value   Pr(>F)   
# (Intercept)     1.2646  1  6.9672 0.011687 * 
# Gender          0.0009  1  0.0048 0.944858   
# Age             0.2490  1  1.3720 0.248228   
# Hours.Final     0.9414  1  5.1865 0.028047 * 
# pH              0.0009  1  0.0049 0.944747   
# Diagnosis       0.1122  2  0.3090 0.735893   
# Cohort          3.8968  8  2.6837 0.018026 * 
# Dissecton.Group 7.3387 11  3.6758 0.001112 **
#   Residuals       7.4415 41                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Very similar to TTP


########

#Going down more of a rabbithole trying to interpret RNA concentration:

setwd("~/Documents/Microarray Gen/FrontalPole/Output/CorrelationsBetweenVariables/AfterBasicQC")


pdf("Scatterplot_TTPConcbyWeight_noBad.pdf", width=6, height=6)
plot(TTP_Average.RNAConc..ng.uL.~Block.Weight..g., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Tissue Block Weight (g)", ylab="TTP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_Average.RNAConc..ng.uL.~Block.Weight..g., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#No correlation - but maybe because of the bad Conc measurements from batches 9 & 10?

pdf("Scatterplot_TZPConcbyWeight_noBad.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL.~Block.Weight..g., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Tissue Block Weight (g)", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

cor(SubjectInfo$TZP_Average.RNAConc..ng.uL.[SubjectInfo$LowQualityRNA==F], SubjectInfo$Block.Weight..g.[SubjectInfo$LowQualityRNA==F])
#[1] 0.1788493


pdf("Scatterplot_TZPConcbyTZP260280_noBad.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL.~TZP_Average.260.280, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="TZP RNA Purity: Average 260/280", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL.~TZP_Average.260.280, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TZPConcbyTZP260230_noBad.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL.~TZP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="TZP RNA Purity: Average 260/230", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL.~TZP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TTPConcbyTTP260280_noBad.pdf", width=6, height=6)
plot(TTP_Average.RNAConc..ng.uL.~TTP_Average.260.280, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="TTP RNA Purity: Average 260/280", ylab="TTP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_Average.RNAConc..ng.uL.~TTP_Average.260.280, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TTPConcbyTTP260230_noBad.pdf", width=6, height=6)
plot(TTP_Average.RNAConc..ng.uL.~TTP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="TTP RNA Purity: Average 260/230", ylab="TTP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_Average.RNAConc..ng.uL.~TTP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

#Oh that is so much fun: you can see how we corrected for tissue weight during the initial extraction, how purified RNA and unpurified RNA have similar concentration... to a point (probably because we're maxing out the spin columns), and how phenol contamination is present in the unpurified RNA (TTP) and is inflating the concentration estimates, whereas in the purified RNA there is some other contamination that is also inflating concentration estimates.


pdf("Scatterplot_TZPConcbyTTP_RIN_noBad.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL.~TTP_BioAnalyzer_RIN, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="TTP RNA Integrity Number (RIN)", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL.~TTP_BioAnalyzer_RIN, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()


#So hypothetically TTP RNA integrity would probably be more closely related to pre- and post-mortem factors, right?  Let's see.


pdf("Scatterplot_Weight_byPMI_noBad.pdf", width=6, height=6)
plot(Block.Weight..g.~Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Post-Mortem Interval (PMI)", ylab="Tissue Block Weight (g)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(Block.Weight..g.~Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TTP_RINbyPMI_noBad.pdf", width=6, height=6)
plot(TTP_BioAnalyzer_RIN~Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Post-Mortem Interval (PMI)", ylab="TTP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TTP_BioAnalyzer_RIN~Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()

pdf("Scatterplot_TZP_RINbyPMI_noBad.pdf", width=6, height=6)
plot(TZP_BioAnalyzer_RIN~Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,], xlab="Post-Mortem Interval (PMI)", ylab="TZP RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_BioAnalyzer_RIN~Hours.Final, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])
abline(ForBestFitLine, col="black", lwd=3)
dev.off()


#################


#Going further down a rabbithole:  This code was all written later to try to double-check that RNA concentration really doesn't relate to diagnosis.

#First, figuring out whether it might be better to correct RNA concentration for other variables influencing its measurement: Block weight, purity, integrity, dissection/extraction group:

#TZP concentration:

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -137.555  -25.028   -6.041   12.189  133.843 
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           3498.552    650.909   5.375 1.27e-06 ***
# Block.Weight..g.                        41.910     40.959   1.023    0.310    
# TZP_Average.260.280                  -1559.379    310.021  -5.030 4.60e-06 ***
# TZP_Average.260.230                    -31.237     39.755  -0.786    0.435    
# TZP_BioAnalyzer_RIN                      9.515     13.388   0.711    0.480    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.     8.341     22.226   0.375    0.709    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 54.25 on 61 degrees of freedom
# Multiple R-squared:  0.3483,	Adjusted R-squared:  0.2949 
# F-statistic:  6.52 on 5 and 61 DF,  p-value: 6.436e-05

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value    Pr(>F)    
# (Intercept)                           48611  1 16.9637 0.0001428 ***
# Block.Weight..g.                       4891  1  1.7067 0.1973902    
# TZP_Average.260.280                   41070  1 14.3321 0.0004115 ***
# TZP_Average.260.230                    2640  1  0.9213 0.3417634    
# TZP_BioAnalyzer_RIN                    2749  1  0.9593 0.3320801    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   1070  1  0.3732 0.5440115    
# Dissecton.Group                       36272 11  1.1507 0.3444648    
# Residuals                            143279 50                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value    Pr(>F)    
# (Intercept)                           75415  1 25.9655 4.401e-06 ***
# Block.Weight..g.                       4636  1  1.5961    0.2118    
# TZP_Average.260.280                   68888  1 23.7185 9.816e-06 ***
# TZP_Average.260.230                    1336  1  0.4599    0.5005    
# TZP_BioAnalyzer_RIN                     730  1  0.2514    0.6181    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.    423  1  0.1456    0.7042    
# Age                                     311  1  0.1070    0.7448    
# pH                                     5174  1  1.7813    0.1875    
# Hours.Final                            3164  1  1.0895    0.3011    
# Gender                                 2300  1  0.7918    0.3774    
# Diagnosis                              9812  2  1.6892    0.1941    
# Residuals                            159743 55                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -129.276  -26.896    0.197   16.683  132.632 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           4006.2791   786.2173   5.096 4.40e-06 ***
# Block.Weight..g.                        60.7906    48.1183   1.263   0.2118    
# TZP_Average.260.280                  -1581.1513   324.6606  -4.870 9.82e-06 ***
# TZP_Average.260.230                    -27.2363    40.1633  -0.678   0.5005    
# TZP_BioAnalyzer_RIN                      7.8075    15.5699   0.501   0.6181    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.     9.6732    25.3479   0.382   0.7042    
# Age                                     -0.1764     0.5393  -0.327   0.7448    
# pH                                     -63.0442    47.2370  -1.335   0.1875    
# Hours.Final                             -1.0557     1.0114  -1.044   0.3011    
# GenderF                                -23.6772    26.6079  -0.890   0.3774    
# DiagnosisBP                             -8.6286    18.0870  -0.477   0.6352    
# DiagnosisSchiz                         -30.0752    16.6081  -1.811   0.0756 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 53.89 on 55 degrees of freedom
# Multiple R-squared:  0.4202,	Adjusted R-squared:  0.3042 
# F-statistic: 3.623 on 11 and 55 DF,  p-value: 0.0006891

#There is a trend towards Schiz being lower.

car::Anova(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]), type=3)

# Anova Table (Type III tests)
# 
# Response: TZP_Average.RNAConc..ng.uL.
# Sum Sq Df F value    Pr(>F)    
# (Intercept)                           45726  1 16.2748 0.0002150 ***
# Block.Weight..g.                       4969  1  1.7685 0.1904242    
# TZP_Average.260.280                   43720  1 15.5608 0.0002834 ***
# TZP_Average.260.230                    3191  1  1.1358 0.2923652    
# TZP_BioAnalyzer_RIN                    1260  1  0.4483 0.5066435    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.    274  1  0.0975 0.7563148    
# Age                                     235  1  0.0837 0.7737598    
# pH                                     2992  1  1.0648 0.3077580    
# Hours.Final                            2849  1  1.0142 0.3194147    
# Gender                                 4940  1  1.7584 0.1916669    
# Diagnosis                             11224  2  1.9974 0.1477978    
# Dissecton.Group                       36120 11  1.1687 0.3356493    
# Residuals                            123623 44                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -117.446  -25.643   -0.789   28.207  127.767 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           3780.9284   937.2180   4.034 0.000215 ***
# Block.Weight..g.                        82.3634    61.9349   1.330 0.190424    
# TZP_Average.260.280                  -1483.5088   376.0744  -3.945 0.000283 ***
# TZP_Average.260.230                    -47.1247    44.2188  -1.066 0.292365    
# TZP_BioAnalyzer_RIN                     13.4370    20.0688   0.670 0.506643    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   -11.7409    37.5994  -0.312 0.756315    
# Age                                     -0.1630     0.5635  -0.289 0.773760    
# pH                                     -53.3181    51.6698  -1.032 0.307758    
# Hours.Final                             -1.1907     1.1824  -1.007 0.319415    
# GenderF                                -38.3069    28.8880  -1.326 0.191667    
# DiagnosisBP                             -6.9929    18.2693  -0.383 0.703735    
# DiagnosisSchiz                         -32.7282    16.7594  -1.953 0.057219 .  
# Dissecton.Group2                       -27.7158    40.8927  -0.678 0.501468    
# Dissecton.Group3                       -13.4469    34.8826  -0.385 0.701731    
# Dissecton.Group4                        -0.2819    37.4020  -0.008 0.994021    
# Dissecton.Group5                        42.8697    32.3157   1.327 0.191489    
# Dissecton.Group6                        -7.3788    42.5895  -0.173 0.863247    
# Dissecton.Group7                       -24.7761    32.5972  -0.760 0.451266    
# Dissecton.Group8                       -22.2718    37.1401  -0.600 0.551801    
# Dissecton.Group9                        11.9438    33.8251   0.353 0.725695    
# Dissecton.Group10                       50.4089    32.0457   1.573 0.122875    
# Dissecton.Group11                        9.1367    32.8864   0.278 0.782448    
# Dissecton.Group12                      -15.1686    30.6580  -0.495 0.623225    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 53.01 on 44 degrees of freedom
# Multiple R-squared:  0.5513,	Adjusted R-squared:  0.3269 
# F-statistic: 2.457 on 22 and 44 DF,  p-value: 0.005521


#There is still a trend towards Schiz being lower, similar to later analyses of housekeeping gene expression.


#Since block weight could reflect genuine biological differences, I'm going to remove it from the analysis and see if anything changes

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -118.815  -26.925    3.934   23.825  115.310 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           4.274e+03  8.681e+02   4.923 1.19e-05 ***
# TZP_Average.260.280                  -1.653e+03  3.569e+02  -4.632 3.10e-05 ***
# TZP_Average.260.230                  -4.594e+01  4.459e+01  -1.030   0.3084    
# TZP_BioAnalyzer_RIN                   1.954e+00  1.827e+01   0.107   0.9153    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  2.381e+01  2.667e+01   0.893   0.3767    
# Age                                   7.165e-02  5.397e-01   0.133   0.8950    
# pH                                   -6.732e+01  5.102e+01  -1.320   0.1937    
# Hours.Final                          -1.729e+00  1.120e+00  -1.543   0.1298    
# GenderF                              -2.315e+01  2.677e+01  -0.865   0.3918    
# DiagnosisBP                          -6.791e+00  1.842e+01  -0.369   0.7141    
# DiagnosisSchiz                       -2.983e+01  1.676e+01  -1.780   0.0818 .  
# Dissecton.Group2                     -6.419e+00  3.795e+01  -0.169   0.8664    
# Dissecton.Group3                     -2.729e+01  3.358e+01  -0.813   0.4206    
# Dissecton.Group4                      1.636e+01  3.555e+01   0.460   0.6477    
# Dissecton.Group5                      4.112e+01  3.256e+01   1.263   0.2132    
# Dissecton.Group6                      3.162e+00  4.220e+01   0.075   0.9406    
# Dissecton.Group7                     -1.897e+01  3.258e+01  -0.582   0.5633    
# Dissecton.Group8                     -3.008e+01  3.698e+01  -0.813   0.4203    
# Dissecton.Group9                      1.940e+01  3.364e+01   0.577   0.5670    
# Dissecton.Group10                     4.568e+01  3.212e+01   1.422   0.1619    
# Dissecton.Group11                     2.609e+00  3.279e+01   0.080   0.9369    
# Dissecton.Group12                    -1.892e+01  3.079e+01  -0.614   0.5421    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 53.46 on 45 degrees of freedom
# Multiple R-squared:  0.5333,	Adjusted R-squared:  0.3154 
# F-statistic: 2.448 on 21 and 45 DF,  p-value: 0.005911


#Same results - so having block weight in the model just cleans up the noise a little bit, it does not mediate the diagnosis effect.

#If we trim the model down just to things that seem to matter (note: this is completely exploratory):
summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+pH+Hours.Final+Gender+Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -130.831  -30.506    1.983   17.258  139.288 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          4056.050    765.756   5.297 1.90e-06 ***
# Block.Weight..g.       69.378     40.450   1.715   0.0917 .  
# TZP_Average.260.280 -1604.654    318.046  -5.045 4.77e-06 ***
# TZP_Average.260.230   -33.401     36.941  -0.904    0.3696    
# pH                    -50.432     44.329  -1.138   0.2599    
# Hours.Final            -1.266      0.942  -1.344   0.1842    
# GenderF               -20.938     24.050  -0.871   0.3876    
# DiagnosisBP            -7.799     17.495  -0.446   0.6574    
# DiagnosisSchiz        -32.395     16.029  -2.021   0.0479 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 52.98 on 58 degrees of freedom
# Multiple R-squared:  0.4092,	Adjusted R-squared:  0.3277 
# F-statistic: 5.021 on 8 and 58 DF,  p-value: 9.278e-05

#Interesting. It's possible that there may be something there.

#If we trim the model down just to things that seem to matter (note: this is completely exploratory) & remove the highly collinear purity measurements which may or may not be spurious:
summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Age+pH+Hours.Final+Gender+Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      532.6737   363.4770   1.465   0.1481  
# Block.Weight..g. 108.0026    47.2655   2.285   0.0259 *
# Age               -0.3238     0.5904  -0.548   0.5855  
# pH               -22.1137    52.6561  -0.420   0.6760  
# Hours.Final       -1.4896     1.1251  -1.324   0.1906  
# GenderF          -52.4569    27.5675  -1.903   0.0619 .
# DiagnosisBP       11.1113    20.5412   0.541   0.5906  
# DiagnosisSchiz   -25.1206    18.9666  -1.324   0.1905  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 62.91 on 59 degrees of freedom
# Multiple R-squared:  0.1524,	Adjusted R-squared:  0.05181 
# F-statistic: 1.515 on 7 and 59 DF,  p-value: 0.1799

#Not as related 


#What about TTP?

#Just the technical variables first:
summary.lm(lm(TTP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.260.280+TTP_Average.260.230+TTP_BioAnalyzer_RIN+TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -359.99  -42.97   -3.78   51.80  320.83 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            -37.37     470.83  -0.079   0.9370    
# Block.Weight..g.                        26.37      73.39   0.359   0.7206    
# TTP_Average.260.280                    482.25     244.68   1.971   0.0533 .  
# TTP_Average.260.230                   -165.23      30.07  -5.495 8.09e-07 ***
# TTP_BioAnalyzer_RIN                    -17.23      23.37  -0.737   0.4639    
# TTP_BioAnalyzer_rRNA.Ratio..28s.18s.    25.49      36.11   0.706   0.4830    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 104.9 on 61 degrees of freedom
# Multiple R-squared:  0.4314,	Adjusted R-squared:  0.3848 
# F-statistic: 9.258 on 5 and 61 DF,  p-value: 1.346e-06


summary.lm(lm(TTP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.260.280+TTP_Average.260.230+TTP_BioAnalyzer_RIN+TTP_BioAnalyzer_rRNA.Ratio..28s.18s.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -193.534  -36.653   -1.893   35.758  200.296 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          -617.423    439.354  -1.405   0.1661    
# Block.Weight..g.                      -55.391     64.484  -0.859   0.3944    
# TTP_Average.260.280                   543.327    220.649   2.462   0.0173 *  
# TTP_Average.260.230                    65.358     44.484   1.469   0.1480    
# TTP_BioAnalyzer_RIN                    -8.709     21.335  -0.408   0.6849    
# TTP_BioAnalyzer_rRNA.Ratio..28s.18s.   24.075     37.866   0.636   0.5278    
# Dissecton.Group2                        9.209     51.906   0.177   0.8599    
# Dissecton.Group3                      -37.902     49.568  -0.765   0.4481    
# Dissecton.Group4                      -27.472     45.817  -0.600   0.5515    
# Dissecton.Group5                       72.507     47.314   1.532   0.1317    
# Dissecton.Group6                        3.839     68.775   0.056   0.9557    
# Dissecton.Group7                       57.424     50.468   1.138   0.2606    
# Dissecton.Group8                       71.892     52.571   1.368   0.1776    
# Dissecton.Group9                      360.279     76.843   4.688 2.16e-05 ***
# Dissecton.Group10                     452.677     75.923   5.962 2.51e-07 ***
# Dissecton.Group11                      69.658     52.181   1.335   0.1879    
# Dissecton.Group12                      16.935     46.271   0.366   0.7159    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 76.37 on 50 degrees of freedom
# Multiple R-squared:  0.7532,	Adjusted R-squared:  0.6742 
# F-statistic: 9.536 on 16 and 50 DF,  p-value: 3.272e-10


summary.lm(lm(TTP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.260.280+TTP_Average.260.230+TTP_BioAnalyzer_RIN+TTP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -126.556  -44.975   -2.012   34.713  177.457 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          -307.2256   625.7502  -0.491 0.625887    
# Block.Weight..g.                      -10.5496    70.7841  -0.149 0.882204    
# TTP_Average.260.280                   700.0894   224.3829   3.120 0.003188 ** 
# TTP_Average.260.230                    44.8320    43.7945   1.024 0.311577    
# TTP_BioAnalyzer_RIN                     8.0010    23.6432   0.338 0.736666    
# TTP_BioAnalyzer_rRNA.Ratio..28s.18s.   -3.0187    39.5593  -0.076 0.939519    
# Age                                    -0.5445     0.7375  -0.738 0.464271    
# pH                                    -89.5559    69.1280  -1.296 0.201901    
# Hours.Final                            -1.2564     1.5450  -0.813 0.420504    
# GenderF                              -102.5333    40.8652  -2.509 0.015862 *  
# DiagnosisBP                            43.1062    25.6396   1.681 0.099803 .  
# DiagnosisSchiz                         12.4178    24.5688   0.505 0.615780    
# Dissecton.Group2                       -3.0412    51.7307  -0.059 0.953387    
# Dissecton.Group3                      -58.2674    49.6701  -1.173 0.247072    
# Dissecton.Group4                       -2.1932    46.5134  -0.047 0.962605    
# Dissecton.Group5                       47.3975    47.0214   1.008 0.318967    
# Dissecton.Group6                        9.3309    68.6109   0.136 0.892444    
# Dissecton.Group7                       45.3414    49.3813   0.918 0.363526    
# Dissecton.Group8                       36.9348    53.0786   0.696 0.490181    
# Dissecton.Group9                      316.1385    76.2439   4.146 0.000152 ***
# Dissecton.Group10                     458.3436    74.7001   6.136 2.13e-07 ***
# Dissecton.Group11                      54.3594    51.2887   1.060 0.294986    
# Dissecton.Group12                       8.1232    45.2042   0.180 0.858214    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 73.91 on 44 degrees of freedom
# Multiple R-squared:  0.7965,	Adjusted R-squared:  0.6948 
# F-statistic:  7.83 on 22 and 44 DF,  p-value: 4.553e-09


#There is actually a trend towards BP having *higher* concentration in the TTP


#So is something happening during the purification process?  e.g., maybe there is a higher percentage of small RNA getting filtered out for the Schiz samples?

#This model is waaaay too large but maybe still interesting:
summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.RNAConc..ng.uL.+TZP_Average.260.280+TZP_Average.260.230+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -120.694  -18.647    0.721   14.703  106.364 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           772.46198  919.49415   0.840  0.40550    
# Block.Weight..g.                      108.12321   48.58799   2.225  0.03136 *  
# TTP_Average.RNAConc..ng.uL.             0.48638    0.09003   5.403 2.68e-06 ***
# TZP_Average.260.280                  -393.55229  356.23908  -1.105  0.27541    
# TZP_Average.260.230                   -26.49229   34.73288  -0.763  0.44978    
# TZP_BioAnalyzer_RIN                    16.89181   15.68104   1.077  0.28739    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  -25.74467   29.46861  -0.874  0.38717    
# Age                                    -0.10867    0.44002  -0.247  0.80611    
# pH                                     17.52829   42.41729   0.413  0.68149    
# Hours.Final                            -0.15477    0.94281  -0.164  0.87038    
# GenderF                               -10.86258   23.11832  -0.470  0.64082    
# DiagnosisBP                           -11.81933   14.29104  -0.827  0.41278    
# DiagnosisSchiz                        -24.27540   13.17749  -1.842  0.07235 .  
# Dissecton.Group2                      -44.28390   32.07254  -1.381  0.17449    
# Dissecton.Group3                      -10.53951   27.23867  -0.387  0.70071    
# Dissecton.Group4                      -29.49238   29.69667  -0.993  0.32621    
# Dissecton.Group5                      -13.82816   27.32500  -0.506  0.61540    
# Dissecton.Group6                      -39.21653   33.76841  -1.161  0.25191    
# Dissecton.Group7                      -61.26850   26.33024  -2.327  0.02475 *  
# Dissecton.Group8                      -53.32913   29.56020  -1.804  0.07823 .  
# Dissecton.Group9                      -99.94421   33.55993  -2.978  0.00475 ** 
# Dissecton.Group10                    -131.54060   41.95392  -3.135  0.00309 ** 
# Dissecton.Group11                     -22.84040   26.34833  -0.867  0.39083    
# Dissecton.Group12                     -34.01109   24.18796  -1.406  0.16688    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 41.38 on 43 degrees of freedom
# Multiple R-squared:  0.7327,	Adjusted R-squared:  0.5898 
# F-statistic: 5.125 on 23 and 43 DF,  p-value: 2.032e-06

#What if instead we correct for impurity/integrity in the TTP measurements?

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -94.884 -16.449  -0.171  15.860 121.079 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          -663.56490  304.57512  -2.179  0.03488 *  
# Block.Weight..g.                       96.33365   34.36784   2.803  0.00756 ** 
# TTP_Average.RNAConc..ng.uL.             0.41703    0.07318   5.699    1e-06 ***
# TTP_Average.260.280                   354.64748  120.36448   2.946  0.00517 ** 
# TTP_Average.260.230                    30.05601   21.50987   1.397  0.16949    
# TTP_BioAnalyzer_RIN                    20.72433   11.49154   1.803  0.07833 .  
# TTP_BioAnalyzer_rRNA.Ratio..28s.18s.  -34.01999   19.20369  -1.772  0.08356 .  
# Age                                    -0.22338    0.36022  -0.620  0.53845    
# pH                                      1.24794   34.18927   0.037  0.97105    
# Hours.Final                            -0.26822    0.75558  -0.355  0.72433    
# GenderF                               -37.00569   21.20793  -1.745  0.08815 .  
# DiagnosisBP                             2.24279   12.83918   0.175  0.86215    
# DiagnosisSchiz                        -15.22630   11.96045  -1.273  0.20984    
# Dissecton.Group2                      -42.30451   25.11146  -1.685  0.09930 .  
# Dissecton.Group3                      -30.13001   24.48441  -1.231  0.22517    
# Dissecton.Group4                      -34.62830   22.57853  -1.534  0.13243    
# Dissecton.Group5                      -43.07428   23.08662  -1.866  0.06891 .  
# Dissecton.Group6                      -58.52536   33.31130  -1.757  0.08605 .  
# Dissecton.Group7                      -42.07161   24.19865  -1.739  0.08926 .  
# Dissecton.Group8                      -63.58811   25.90618  -2.455  0.01823 *  
# Dissecton.Group9                      -23.37888   43.64512  -0.536  0.59495    
# Dissecton.Group10                     -29.59035   49.39399  -0.599  0.55227    
# Dissecton.Group11                     -12.69080   25.21173  -0.503  0.61727    
# Dissecton.Group12                     -34.10539   21.95055  -1.554  0.12758    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 35.88 on 43 degrees of freedom
# Multiple R-squared:  0.7991,	Adjusted R-squared:  0.6916 
# F-statistic: 7.436 on 23 and 43 DF,  p-value: 1.03e-08

#That's model produces a better fit, but TZP has less of a relationship with Schiz.



#Or both TZP and TTP purity?

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.RNAConc..ng.uL.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TZP_Average.260.230+TTP_BioAnalyzer_RIN+TTP_BioAnalyzer_rRNA.Ratio..28s.18s.+Age+pH+Hours.Final+Gender+Diagnosis+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -92.980 -17.898  -0.693  17.539 120.346 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          -365.86643  980.23100  -0.373   0.7109    
# Block.Weight..g.                       95.61423   35.49171   2.694   0.0102 *  
# TTP_Average.RNAConc..ng.uL.             0.40390    0.08333   4.847 1.83e-05 ***
# TTP_Average.260.280                   347.12614  137.15292   2.531   0.0153 *  
# TTP_Average.260.230                    26.79336   23.90164   1.121   0.2688    
# TZP_Average.260.280                  -117.82582  345.35070  -0.341   0.7347    
# TZP_Average.260.230                    -0.50462   34.87309  -0.014   0.9885    
# TTP_BioAnalyzer_RIN                    20.53778   11.80989   1.739   0.0895 .  
# TTP_BioAnalyzer_rRNA.Ratio..28s.18s.  -34.80135   20.54271  -1.694   0.0978 .  
# Age                                    -0.22103    0.37000  -0.597   0.5535    
# pH                                     -2.66988   36.67276  -0.073   0.9423    
# Hours.Final                            -0.31480    0.78974  -0.399   0.6922    
# GenderF                               -37.00525   21.68775  -1.706   0.0955 .  
# DiagnosisBP                             1.48835   13.44624   0.111   0.9124    
# DiagnosisSchiz                        -16.38577   12.81560  -1.279   0.2082    
# Dissecton.Group2                      -40.91059   26.61256  -1.537   0.1319    
# Dissecton.Group3                      -28.53316   25.43589  -1.122   0.2685    
# Dissecton.Group4                      -31.32540   25.13067  -1.247   0.2197    
# Dissecton.Group5                      -39.84696   25.48155  -1.564   0.1256    
# Dissecton.Group6                      -56.36885   34.60530  -1.629   0.1110    
# Dissecton.Group7                      -40.64155   25.06347  -1.622   0.1126    
# Dissecton.Group8                      -62.48213   26.71828  -2.339   0.0243 *  
# Dissecton.Group9                      -25.22144   45.19671  -0.558   0.5799    
# Dissecton.Group10                     -27.90966   50.74279  -0.550   0.5853    
# Dissecton.Group11                     -12.19798   26.30206  -0.464   0.6453    
# Dissecton.Group12                     -33.00178   23.05680  -1.431   0.1599    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 36.69 on 41 degrees of freedom
# Multiple R-squared:  0.7997,	Adjusted R-squared:  0.6776 
# F-statistic: 6.548 on 25 and 41 DF,  p-value: 7.382e-08

#Looks like only technical variables matter.

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.RNAConc..ng.uL.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TZP_Average.260.230+TTP_BioAnalyzer_RIN+TTP_BioAnalyzer_rRNA.Ratio..28s.18s.+Dissecton.Group, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residual standard error: 36.16 on 47 degrees of freedom
# Multiple R-squared:  0.7769,	Adjusted R-squared:  0.6867 
# F-statistic: 8.613 on 19 and 47 DF,  p-value: 1.113e-09

#...And those are largely captured by a smaller model.

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.RNAConc..ng.uL.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TZP_Average.260.230+TTP_BioAnalyzer_RIN+TTP_BioAnalyzer_rRNA.Ratio..28s.18s., data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residual standard error: 35.96 on 58 degrees of freedom
# Multiple R-squared:  0.7278,	Adjusted R-squared:  0.6902 
# F-statistic: 19.38 on 8 and 58 DF,  p-value: 8.159e-14

#...Or an even smaller model.

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.RNAConc..ng.uL.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TTP_BioAnalyzer_RIN, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residual standard error: 36.14 on 60 degrees of freedom
# Multiple R-squared:  0.7156,	Adjusted R-squared:  0.6871 
# F-statistic: 25.16 on 6 and 60 DF,  p-value: 1.083e-14

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TTP_BioAnalyzer_RIN, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residual standard error: 50.77 on 61 degrees of freedom
# Multiple R-squared:  0.4293,	Adjusted R-squared:  0.3825 
# F-statistic: 9.178 on 5 and 61 DF,  p-value: 1.497e-06

#For which TTP concentration is a major contributor (although purity still matters - esp. once TTP concentration is out of the model)

#But neither reduced model produces a relationship with diagnosis:

summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.RNAConc..ng.uL.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TTP_BioAnalyzer_RIN+Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -103.908  -17.751   -1.097   14.793  126.301 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  322.73449  602.38090   0.536  0.59417    
# Block.Weight..g.              53.64489   25.67097   2.090  0.04104 *  
# TTP_Average.RNAConc..ng.uL.    0.39245    0.05082   7.723 1.79e-10 ***
# TTP_Average.260.280          245.73077   89.40671   2.748  0.00797 ** 
# TTP_Average.260.230           36.46041   12.41523   2.937  0.00475 ** 
# TZP_Average.260.280         -389.40674  252.50699  -1.542  0.12847    
# TTP_BioAnalyzer_RIN           10.25764    7.37082   1.392  0.16934    
# DiagnosisBP                   -9.58712   11.20716  -0.855  0.39582    
# DiagnosisSchiz               -16.07853   11.36353  -1.415  0.16244    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 36.12 on 58 degrees of freedom
# Multiple R-squared:  0.7253,	Adjusted R-squared:  0.6874 
# F-statistic: 19.14 on 8 and 58 DF,  p-value: 1.052e-13

#No relationship with diagnosis still...


summary.lm(lm(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+TTP_Average.260.280+TTP_Average.260.230+TZP_Average.260.280+TTP_BioAnalyzer_RIN+Diagnosis, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,]))

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -157.444  -23.590   -1.352   16.579  146.026 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          2429.677    758.360   3.204  0.00219 ** 
# Block.Weight..g.       50.212     36.244   1.385  0.17114    
# TTP_Average.260.280   328.812    125.331   2.624  0.01106 *  
# TTP_Average.260.230   -27.979     12.981  -2.155  0.03523 *  
# TZP_Average.260.280 -1328.941    312.444  -4.253 7.62e-05 ***
# TTP_BioAnalyzer_RIN     5.841     10.377   0.563  0.57564    
# DiagnosisBP            -1.520     15.756  -0.096  0.92348    
# DiagnosisSchiz        -17.479     16.044  -1.089  0.28039    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 51.01 on 59 degrees of freedom
# Multiple R-squared:  0.4428,	Adjusted R-squared:  0.3766 
# F-statistic: 6.697 on 7 and 59 DF,  p-value: 7.658e-06

#Nope.

#... so what if we corrected our RNA concentration estimates for impurity, and then looked at relationships with later gene expression (esp. housekeeping gene expression)?  Is RNA concentration more related to gene expression if we correct the measurement for purity?

#I removed the one subject with super low 260/230 for this analysis: 

TTP_Average.RNAConc..ng.uL._PurityCorrected<-lm(TTP_Average.RNAConc..ng.uL.~TTP_Average.260.280+TTP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])$residuals+mean(SubjectInfo$TTP_Average.RNAConc..ng.uL.[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F], na.rm=T)

TZP_Average.RNAConc..ng.uL._PurityCorrected<-lm(TZP_Average.RNAConc..ng.uL.~TZP_Average.260.280+TZP_Average.260.230, data=SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F,])$residuals+mean(SubjectInfo$TZP_Average.RNAConc..ng.uL.[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F], na.rm=T)


pdf("Scatterplot_TZPConcbyTTPConc_PurityCorrected.pdf", width=6, height=6)
plot(TZP_Average.RNAConc..ng.uL._PurityCorrected~TTP_Average.RNAConc..ng.uL._PurityCorrected, xlab="TTP Average RNA Concentration (ng/uL, purity corrected)", ylab="TZP Average RNA Concentration (ng/uL, purity corrected)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(TZP_Average.RNAConc..ng.uL._PurityCorrected~TTP_Average.RNAConc..ng.uL._PurityCorrected)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

#hmm... that is actually a worse correlation than before.  Maybe only correct the TTP concentration?

pdf("Scatterplot_TZPConcbyTTPConc_OnlyTTPPurityCorrected.pdf", width=6, height=6)
plot(SubjectInfo$TZP_Average.RNAConc..ng.uL.[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F]~TTP_Average.RNAConc..ng.uL._PurityCorrected, xlab="TTP Average RNA Concentration (ng/uL, purity corrected)", ylab="TZP Average RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(SubjectInfo$TZP_Average.RNAConc..ng.uL.[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==F]~TTP_Average.RNAConc..ng.uL._PurityCorrected)
abline(ForBestFitLine, col="black", lwd=3)
abline(a=0, b=1, col="red", lwd=3)
dev.off()

#This is old code - I don't think I'm seeing much of a diagnosis effect, so I didn't bother updating it.

# pdf("Boxplot_TZPConcPurityCorrected_ByDiagnosis_noBadRIN.pdf", width=4, height=6)
# boxplot(TZP_Average.RNAConc..ng.uL._PurityCorrected~SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F & SubjectInfo$TZP_Average.260.230>1], xlab="Diagnosis", ylab="TZP Average RNA Concentration (ng/uL, purity corrected)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
# stripchart(TZP_Average.RNAConc..ng.uL._PurityCorrected~SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F & SubjectInfo$TZP_Average.260.230>1], xlab="Diagnosis", ylab="TZP Average RNA Concentration (ng/uL, purity corrected)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
# dev.off()
# 
# pdf("Boxplot_TTPConcPurityCorrected_ByDiagnosis_noBadRIN.pdf", width=4, height=6)
# boxplot(TTP_Average.RNAConc..ng.uL._PurityCorrected~SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F & SubjectInfo$TZP_Average.260.230>1], xlab="Diagnosis", ylab="TTP Average RNA Concentration (ng/uL, purity corrected)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
# stripchart(TTP_Average.RNAConc..ng.uL._PurityCorrected~SubjectInfo$Diagnosis[SubjectInfo$LowQualityRNA==F & SubjectInfo$TZP_Average.260.230>1], xlab="Diagnosis", ylab="TTP Average RNA Concentration (ng/uL, purity corrected)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
# dev.off()


#... so if TZP RNA concentration is inflated by contamination (picked up by with 260/280 ratio) and also related to contamination of the TTP samples (which artificially inflated the RNAconc)...which values should we use as a control in our gene expression measurements?  RNAconcentration is traditional - hypothetically the contamination in the TTP that decreased the final TZP concentration should be accounted for in the final RNA concentration value.
#But maybe correcting for the effects of contamination in the final RNA concentration value would make sense?

#That said, TZP concentration is most likely to relate to gene expression, and it seems like the TTP concentraton measurements is what is most affected by impurity (given the small variability in TZP 260/280, I suspect the correlation that I am finding with TZP concentration is spurious)... so maybe correcting for purity doesn't make sense?

#I should come back to this later.


#########################