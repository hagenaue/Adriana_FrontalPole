#Code for calculating NegDeltaCq & then doing QC/PCA again:
#Now with the final 3 outliers removed.

##############################

#Scaling -DeltaCq and Running PCA again:
#Note - I'm reusing code and objects for this analysis for efficiency sake - I have just changed the output directory and saved the workspace separately:

#While scaling I'm going to also change the orientation back to match my previous code:
GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed<-t(scale(GabaGlu_NegDeltaCq_AllSubjects_QCed2, center=T, scale=T))
head(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
#Sanity check:
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, 1, function(y) mean(y, na.rm=T))
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, 1, function(y) sd(y, na.rm=T))
#Looks good.

GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed<-cor(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")

boxplot(GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed)

setwd("~/Documents/Microarray Gen/FrontalPole/Output/GabaGlu/GabaGlu_NegDeltaPCA_NoOutliers")

pdf("GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed.pdf", height=14, width=14)
heatmap(GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed)
dev.off()
#Again - some of the replicates most closely resemble each other, but certainly not all. There are correlation blocks (batch effects??).


#Dealing with the NA values and looking at overall variation in -DeltaCq values across samples:

#Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed[is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)]<-0
sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed))
#[1] 0

pdf("GabaGlu_Boxplot_NegDeltaCqZscore_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq Z-scores per sample for all genes", xlab="Sample ID", ylab="-DeltaCq Z-score")
dev.off()

pdf("GabaGlu_Boxplot_NegDeltaCqZscore_perSample_AllGenes_ylim6.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq Z-scores per sample for all genes", xlab="Sample ID", ylab="-DeltaCq Z-score", ylim=c(-6.5,6.5))
dev.off()

pdf("Boxplot_NegDeltaCqZscore_perSample_HousekeepingGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed[83:93,]), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq Z-scores per sample for housekeeping genes", xlab="Sample ID", ylab="-DeltaCq Z-score")
dev.off()



#Trying PCA:

pca_GabaGlu<-prcomp(t(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed))
tmp<-pca_GabaGlu$x[,1:10]
dim(tmp)
#[1]  133  10
rownames(tmp)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_GabaGlu_NegaDeltaCq.csv")

tmp<-pca_GabaGlu$rotation[,1:10]
write.csv(tmp, "pca_GabaGlu_NegDeltaCq_Eigenvectors.csv")

png("PCA_ScreePlot_GabaGlu_NegDeltaCq.png")
plot(summary(pca_GabaGlu)$importance[2,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 explains a smaller portion of variance now (~27%)


png("PCA_ScreePlot2_GabaGlu_NegDeltaCq.png")
plot(summary(pca_GabaGlu)$importance[3,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()

png("PC1vsPC2_GabaGlu_NegDeltaCq_byDiagnosis.png")
plot(pca_GabaGlu$x[,1]~pca_GabaGlu$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis))
dev.off()
#PC1 Not diagnosis
#No outliers now.
hist(pca_GabaGlu$x[,1], breaks=20)

pdf("PC1_GabaGlu_NegDeltaCq_byMeanHousekeeping.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_MeanHousekeeping_QCed2)
dev.off()
#Again, outliers are gone, so is correlation with HK expression.

colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Temp<-cbind(pca_GabaGlu$x[,c(1:4)],GabaGlu_MeanHousekeeping_QCed2, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3[,c(7,9,11,16,24:28,38:42)] )
cor(Temp)
write.csv(cor(Temp), "GabaGlu_CorMatrix_PCA_GabaGlu_NegDeltaCq_VsSubjectVar.csv")

#Potentially worth exploring:
#PC1 vs. GabaGlu_MeanHousekeeping_QCed2, Age, pH, PMI, TZP_BioAnalyzer_rRNA.Ratio..28s.18s.,TZP_Average.260.280
#PC2 vs. GabaGlu_MeanHousekeeping_QCed2, pH, PMI, RIN, TZP_BioAnalyzer_rRNA.Ratio..28s.18s.,TZP_Average.260.280
#PC3 vs. GabaGlu_MeanHousekeeping_QCed2, Age, pH,  TZP_BioAnalyzer_rRNA.Ratio..28s.18s.


#To be thorough, I came back and looped the basic statistical output for each variable in relationship to PC1-4:
colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats<-matrix(0, 20, 4)
row.names(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3[,c(2,5:9,11,14,16,24:28,38:42,50)])
colnames(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats)<-c("PC1_ANOVA_Pval","PC2_ANOVA_Pval","PC3_ANOVA_Pval","PC4_ANOVA_Pval")

Temp<-data.frame(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3[,c(2,5:9,11,14,16,24:28,38:42,50)], pca_GabaGlu$x[,c(1:4)], ID=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3[,1])

for(i in c(1:20)){
  
  Model<-lmer(Temp$PC1~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats[i,1]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp$PC2~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats[i,2]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp$PC3~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats[i,3]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp$PC4~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats[i,4]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
}

setwd("~/Documents/Microarray Gen/FrontalPole/Output/GabaGlu/GabaGlu_NegDeltaPCA_NoOutliers")
write.csv(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats, "GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats.csv")


#I came back and looped a statistical output for all of these

#Technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)
# (Intercept)                          0.0016  1     0.9684
# Diagnosis                            0.4606  2     0.7943
# Block.Weight..g.                     0.4398  1     0.5072
# TZP_BioAnalyzer_RIN                  0.0730  1     0.7871
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 2.6114  1     0.1061
# TZP_Average.RNAConc..ng.uL.          0.6769  1     0.4106
# TZP_Average.260.230                  0.4164  1     0.5187
# TZP_Average.260.280                  0.0060  1     0.9385

#Basically no relationship.

#Biological variables:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)  5.3714  1  0.0204699 *  
#   Diagnosis    0.8872  2  0.6417116    
# Age         14.4563  1  0.0001434 ***
#   Gender       0.1385  1  0.7097992    
# pH           6.1323  1  0.0132735 *  
#   Hours.Final  7.0943  1  0.0077328 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Biological and technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           1.6535  1  0.1984803    
# Diagnosis                             0.3451  2  0.8415235    
# Age                                  18.9823  1  1.319e-05 ***
# Gender                                0.0040  1  0.9498084    
# pH                                    5.0026  1  0.0253090 *  
# Hours.Final                          14.4961  1  0.0001405 ***
# Block.Weight..g.                      0.1518  1  0.6968143    
# TZP_BioAnalyzer_RIN                   4.4325  1  0.0352609 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.1928  1  0.6605921    
# TZP_Average.RNAConc..ng.uL.           2.8477  1  0.0915029 .  
# TZP_Average.260.230                   2.3593  1  0.1245395    
# TZP_Average.260.280                   0.1535  1  0.6951698    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - RNA metrics less after controlling for HK expression, but other biological variables matter more.

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   9.5620  1   0.001987 ** 
# Diagnosis                     3.2166  2   0.200231    
# Age                          21.2015  1  4.134e-06 ***
# Gender                        0.3318  1   0.564601    
# pH                            5.7439  1   0.016546 *  
# Hours.Final                   7.7398  1   0.005402 ** 
# TZP_BioAnalyzer_RIN           2.9564  1   0.085539 .  
# TZP_Average.RNAConc..ng.uL.   8.4148  1   0.003722 ** 
# Card                        240.4470 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting -  RNAconc matters more now.

#Biological and more technical variables & card
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            0.2580  1    0.61148    
# Diagnosis                              2.3014  2    0.31641    
# Age                                   18.1558  1  2.035e-05 ***
# Gender                                 0.0681  1    0.79406    
# pH                                     3.2378  1    0.07196 .  
# Hours.Final                            5.8413  1    0.01565 *  
# TZP_BioAnalyzer_RIN                    1.1790  1    0.27757    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   1.0487  1    0.30580    
# TZP_Average.RNAConc..ng.uL.            3.5084  1    0.06106 .  
# TZP_Average.260.230                    4.1677  1    0.04120 *  
# TZP_Average.260.280                    0.0242  1    0.87649    
# Card                                 246.9030 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - purity seems to matter now (???).

#With dissection group instead of card:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           0.2400  1   0.624193    
# Diagnosis                             1.3534  2   0.508290    
# Age                                  16.0079  1  6.308e-05 ***
# Gender                                1.7429  1   0.186772    
# pH                                    3.6140  1   0.057294 .  
# Hours.Final                           6.8457  1   0.008885 ** 
# TZP_BioAnalyzer_RIN                   2.1572  1   0.141905    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.7268  1   0.393923    
# TZP_Average.RNAConc..ng.uL.           3.4962  1   0.061511 .  
# TZP_Average.260.230                   5.1061  1   0.023841 *  
# TZP_Average.260.280                   0.0778  1   0.780237    
# Dissecton.Group                      14.4673 12   0.271866    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Again, Dissection group doesn't seem to matter, but purity might.



#PC1 varies with Age, pH, PMI and Card. It no longer varies with RNA concentration (maybe? or maybe it does?), RIN, and dissection - maybe due to normalizing with the housekeeping genes?  

png("PC1_GabaGlu_NegDeltaCq_byDiagnosis.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, ylab="PC1")
dev.off()
#pretty unconvincing
#Notably, those outlier samples are outliers within the control distribution too

pdf("PC1_GabaGlu_NegDeltaCq_byDissection.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Dissecton.Group)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Hours.Final)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byPh.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$pH)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byAge.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Age)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byGender.pdf", width=5, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Gender)
dev.off()

#I'm going to check on some other metrics that didn't matter when looking at Cq:
pdf("PC1_GabaGlu_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC1_GabaGlu_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC1_GabaGlu_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.230)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.)
dev.off()
#skewed distribution - careful

pdf("PC1_GabaGlu_NegDeltaCq_byCard.pdf", width=20, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card, las=3, par(cex.axis=0.75))
dev.off()
#Card is the most convincing so far, but that may just be because there are so few samples on a card the distributions can be tight.
#Although replicate cards sometimes do not look so similar...


#hmmm... controlling for card explicitly would use up a whole bunch of degrees of freedom. :(
#It also seems somewhat redundant with just controlling for ID.

#But it does seem to clean things up a lot...

#I'm nervous about including it though, because there are so few samples per card that making the card effect essentially 0 seems like it would be removing a lot of the subject-level variability, although I suppose if we're simultaneously controlling for diagnosis it is less of an issue?
#Alternatively, we could make the gene values average out to the mean for both replicate cards.


#PC2:


#Technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)                          0.0203  1     0.8866  
# Diagnosis                            2.1863  2     0.3352  
# Block.Weight..g.                     0.0939  1     0.7593  
# TZP_BioAnalyzer_RIN                  4.2011  1     0.0404 *
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 0.2162  1     0.6420  
# TZP_Average.RNAConc..ng.uL.          2.1200  1     0.1454  
# TZP_Average.260.230                  0.1525  1     0.6962  
# TZP_Average.260.280                  0.0026  1     0.9595  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#RIN may matter a little bit for PC2

#Biological variables:
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 11.3983  1  0.0007351 ***
# Diagnosis    1.3447  2  0.5105004    
# Age          2.9256  1  0.0871854 .  
# Gender       0.6603  1  0.4164692    
# pH          10.4225  1  0.0012449 ** 
# Hours.Final  2.0618  1  0.1510307    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#pH is strongly related to PC2, maybe Age

#Biological and technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)                           4.3180  1   0.037712 * 
# Diagnosis                             0.9155  2   0.632710   
# Age                                   0.8694  1   0.351111   
# Gender                                0.2600  1   0.610108   
# pH                                   10.8104  1   0.001009 **
# Hours.Final                           3.6522  1   0.055995 . 
# Block.Weight..g.                      0.6208  1   0.430732   
# TZP_BioAnalyzer_RIN                   6.5065  1   0.010748 * 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  2.0460  1   0.152611   
# TZP_Average.RNAConc..ng.uL.           0.5035  1   0.477965   
# TZP_Average.260.230                   0.1578  1   0.691169   
# TZP_Average.260.280                   0.7316  1   0.392380   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                  12.4907  1   0.000409 ***
# Diagnosis                     0.2544  2   0.880537    
# Age                           4.0842  1   0.043286 *  
# Gender                        0.8389  1   0.359710    
# pH                            7.7670  1   0.005321 ** 
# Hours.Final                   2.8940  1   0.088911 .  
# TZP_BioAnalyzer_RIN           8.5571  1   0.003442 ** 
# TZP_Average.RNAConc..ng.uL.   7.2067  1   0.007263 ** 
# Card                        134.3669 35  1.444e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#PC2 is related to RIN, pH, PMI, CARD

#Biological and more technical variables & card
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            5.9258  1   0.014921 *  
#   Diagnosis                              0.4868  2   0.783952    
# Age                                    3.0048  1   0.083017 .  
# Gender                                 0.3140  1   0.575228    
# pH                                    10.1446  1   0.001447 ** 
# Hours.Final                            4.2870  1   0.038405 *  
# TZP_BioAnalyzer_RIN                    6.9740  1   0.008270 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   0.0272  1   0.869032    
# TZP_Average.RNAConc..ng.uL.            1.2954  1   0.255055    
# TZP_Average.260.230                    4.2648  1   0.038909 *  
# TZP_Average.260.280                    1.3905  1   0.238328    
# Card                                 138.7427 35  2.721e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Maybe purity again???

#With dissection group instead of card:
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)                           2.6907  1   0.100936   
# Diagnosis                             0.5368  2   0.764604   
# Age                                   0.3704  1   0.542784   
# Gender                                0.0187  1   0.891183   
# pH                                    4.8505  1   0.027638 * 
# Hours.Final                           6.8363  1   0.008932 **
# TZP_BioAnalyzer_RIN                   9.1447  1   0.002494 **
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  1.6222  1   0.202785   
# TZP_Average.RNAConc..ng.uL.           2.0010  1   0.157196   
# TZP_Average.260.230                   0.3055  1   0.580458   
# TZP_Average.260.280                   0.4301  1   0.511940   
# Dissecton.Group                      25.0394 12   0.014637 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Dissection group might matter a little, but not as much as card.

png("PC2_GabaGlu_NegDeltaCq_byDiagnosis.png")
boxplot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, ylab="PC2")
dev.off()
#Nope

pdf("PC2_GabaGlu_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$pH)
dev.off()

pdf("PC2_GabaGlu_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Hours.Final)
dev.off()

pdf("PC2_GabaGlu_NegDeltaCq_byRIN.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_RIN)
dev.off()
#Pretty.

pdf("PC2_GabaGlu_NegDeltaCq_byRNAconc.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.)
dev.off()
#RNA concentration looks like a nice relationship, but with a few high leverage points.


#Double-checking some of the other metrics that didn't matter for Cq:

pdf("PC2_GabaGlu_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC2_GabaGlu_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC2_GabaGlu_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.230)
dev.off()

pdf("PC2_GabaGlu_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()

#Technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           4.5195  1   0.033511 *  
# Diagnosis                             0.3842  2   0.825241    
# Block.Weight..g.                      0.9160  1   0.338535    
# TZP_BioAnalyzer_RIN                   9.2629  1   0.002338 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 17.7638  1  2.501e-05 ***
# TZP_Average.RNAConc..ng.uL.           3.8723  1   0.049091 *  
# TZP_Average.260.230                   1.8169  1   0.177678    
# TZP_Average.260.280                   4.4439  1   0.035025 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Huh. That's interesting - strong relationship with TZP_BioAnalyzer_rRNA.Ratio..28s.18s, slight relationship with purity.

#Biological variables:
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)  3.4038  1    0.06505 .  
# Diagnosis    0.4256  2    0.80831    
# Age         20.8419  1  4.988e-06 ***
# Gender       0.0714  1    0.78928    
# pH           4.4999  1    0.03390 *  
# Hours.Final  3.9680  1    0.04637 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Strong relationship with age.

#Biological and technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           2.7423  1  0.0977227 .  
# Diagnosis                             1.4324  2  0.4886106    
# Age                                  13.7455  1  0.0002093 ***
# Gender                                1.4598  1  0.2269636    
# pH                                    1.1740  1  0.2785776    
# Hours.Final                           1.6712  1  0.1960939    
# Block.Weight..g.                      0.5943  1  0.4407799    
# TZP_BioAnalyzer_RIN                   4.9671  1  0.0258335 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  8.4037  1  0.0037447 ** 
# TZP_Average.RNAConc..ng.uL.           3.1547  1  0.0757077 .  
# TZP_Average.260.230                   0.6438  1  0.4223501    
# TZP_Average.260.280                   4.7174  1  0.0298586 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   1.0075  1    0.31551    
# Diagnosis                     1.0355  2    0.59586    
# Age                          25.5637  1   4.28e-07 ***
# Gender                        0.3922  1    0.53116    
# pH                            4.2466  1    0.03933 *  
# Hours.Final                   2.8009  1    0.09421 .  
# TZP_BioAnalyzer_RIN           0.8971  1    0.34356    
# TZP_Average.RNAConc..ng.uL.   0.5528  1    0.45716    
# Card                        350.6892 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#No relationship with RIN or RNAconc by itself. Huh.


#Biological and more technical variables & card
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            2.7703  1   0.096031 .  
# Diagnosis                              2.8249  2   0.243541    
# Age                                   16.1514  1  5.848e-05 ***
#   Gender                                 0.0075  1   0.931049    
# pH                                     0.3791  1   0.538072    
# Hours.Final                            0.4374  1   0.508356    
# TZP_BioAnalyzer_RIN                    8.1942  1   0.004202 ** 
#   TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  15.8181  1  6.973e-05 ***
#   TZP_Average.RNAConc..ng.uL.            4.8326  1   0.027926 *  
#   TZP_Average.260.230                    0.0023  1   0.961572    
# TZP_Average.260.280                    3.3334  1   0.067886 .  
# Card                                 350.7672 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#With dissection group instead of card:
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#Dissection group definitely still doesn't matter.

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)                           1.9601  1   0.161500   
# Diagnosis                             1.1144  2   0.572800   
# Age                                  10.3171  1   0.001318 **
# Gender                                0.3070  1   0.579555   
# pH                                    1.3419  1   0.246698   
# Hours.Final                           0.7878  1   0.374762   
# TZP_BioAnalyzer_RIN                   6.2154  1   0.012665 * 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  7.6532  1   0.005667 **
# TZP_Average.RNAConc..ng.uL.           3.2043  1   0.073446 . 
# TZP_Average.260.230                   0.0065  1   0.935754   
# TZP_Average.260.280                   3.0436  1   0.081054 . 
# Dissecton.Group                       5.8479 12   0.923555   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


png("PC3_GabaGlu_NegDeltaCq_byDiagnosis.png")
boxplot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, ylab="PC2")
dev.off()
#Nope

pdf("PC3_GabaGlu_NegDeltaCq_byAge.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Age)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$pH)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Hours.Final)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byRIN.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_RIN)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.)
dev.off()
#The relationship with RNAconc seems driven by a few subjects with low concentration

#Double-checking some of the other metrics that didn't matter for Cq:

pdf("PC3_GabaGlu_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC3_GabaGlu_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC3_GabaGlu_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.230)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()
#Reasonably convincing


#Is there some way to collapse card? (so that it doesn't take up so many df in the model???)

dim(pca_GabaGlu$x)

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$DissectionGroupCard<-paste(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Dissecton.Group,SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card)

pca_GabaGlu_byDissectionCard<-matrix(0, 52, 10)
row.names(pca_GabaGlu_byDissectionCard)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$DissectionGroupCard))

for(i in c(1:10)){
  pca_GabaGlu_byDissectionCard[,i]<-tapply(X=pca_GabaGlu$x[,i], INDEX=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$DissectionGroupCard, FUN=mean)
}

#There are some very strong blocks - let's take a peek
pca_GabaGlu_byDissectionCard_CorMatrix<-cor(t(pca_GabaGlu_byDissectionCard))
dim(pca_GabaGlu_byDissectionCard_CorMatrix)
row.names(pca_GabaGlu_byDissectionCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$DissectionGroupCard))
colnames(pca_GabaGlu_byDissectionCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$DissectionGroupCard))

pdf("Heatmap_pca_GabaGlu_byDissectionCard_CorMatrix.pdf", width=10, height=10)
heatmap(pca_GabaGlu_byDissectionCard_CorMatrix)
dev.off()

pca_GabaGlu_byDissectionCard_CorMatrix

write.csv(pca_GabaGlu_byDissectionCard_CorMatrix, "pca_GabaGlu_byDissectionCard_NegDeltaCq_CorMatrix.csv")

pca_GabaGlu_byDissection<-matrix(0, 13, 10)
row.names(pca_GabaGlu_byDissection)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Dissecton.Group))

for(i in c(1:10)){
  pca_GabaGlu_byDissection[,i]<-tapply(X=pca_GabaGlu$x[,i], INDEX=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Dissecton.Group, FUN=mean)
}

#There are some very strong blocks - let's take a peek
pca_GabaGlu_byDissection_CorMatrix<-cor(t(pca_GabaGlu_byDissection))
dim(pca_GabaGlu_byDissection_CorMatrix)
row.names(pca_GabaGlu_byDissection_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Dissecton.Group))
colnames(pca_GabaGlu_byDissection_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Dissecton.Group))

pdf("Heatmap_pca_GabaGlu_byDissection_CorMatrix.pdf", width=10, height=10)
heatmap(pca_GabaGlu_byDissection_CorMatrix)
dev.off()

pca_GabaGlu_byDissection_CorMatrix
#The correlation blocks here are not as profound as for the other matrix.

row.names(pca_GabaGlu_byDissection)#Matches numeric order
plot(pca_GabaGlu_byDissection[4,]~pca_GabaGlu_byDissection[6,])
plot(pca_GabaGlu_byDissection[7,]~pca_GabaGlu_byDissection[9,])
#I'm feeling... underwhelmed.

pca_GabaGlu_byCard<-matrix(0, 36, 10)
row.names(pca_GabaGlu_byCard)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card))

for(i in c(1:10)){
  pca_GabaGlu_byCard[,i]<-tapply(X=pca_GabaGlu$x[,i], INDEX=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card, FUN=mean)
}


pca_GabaGlu_byCard_CorMatrix<-cor(t(pca_GabaGlu_byCard))
dim(pca_GabaGlu_byCard_CorMatrix)
row.names(pca_GabaGlu_byCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card))
colnames(pca_GabaGlu_byCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card))

pdf("Heatmap_pca_GabaGlu_byCard_CorMatrix.pdf", width=10, height=10)
heatmap(pca_GabaGlu_byCard_CorMatrix)
dev.off()
#There are two really clear blocks, maybe a few additional cards that don't fit in the larger blocks
#The blocks are quite notably not defined by replicate cards - actually, the two blocks are quite clearly defined by odd vs. even cards.

pca_GabaGlu_byCard_CorMatrix
#E.g., Card 18 and Card 12 have R=0.86 but are not replicate cards.

#Example plot:
row.names(pca_GabaGlu_byCard)
plot(pca_GabaGlu_byCard[4,]~pca_GabaGlu_byCard[10,])
#the correlation is almost entirely driven by what I'm assuming are PC1 and PC2 or PC2.

write.csv(pca_GabaGlu_byCard_CorMatrix, "pca_GabaGlu_byCard_NegDeltaCq_CorMatrix.csv")
#Alright I peaked at the cor matrix, and the correlations are strong and definitely in blocks.
#Maybe we can effectively reduce this down.


SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2[SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2%in%c("Card 24.eds", "Card 4.eds", "Card 2.eds", "Card 28.eds", "Card 30.eds", "Card 22.eds", "Card 26.eds", "Card 12.eds", "Card 18.eds", "Card 20.eds", "Card 34.eds", "Card 32.eds", "Card 36.eds", "Card 8.eds", "Card 14.eds", "Card 16.eds", "Card 10.eds", "Card 6.eds")]<-"EvenCards"

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2[SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2%in%c("Card 27.eds", "Card 9.eds", "Card 13.eds", "Card 29.eds", "Card 31.eds", "Card 21.eds", "Card 1.eds", "Card 15.eds", "Card 23.eds", "Card 25.eds", "Card 7.eds", "Card 19.eds", "Card 5.eds", "Card 33.eds", "Card 35.eds", "Card 17.eds", "Card 11.eds", "Card 3.eds")]<-"OddCards"

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2)
# EvenCards  OddCards 
# 67        69 


#Let's look at the results:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                  9.5242  1   0.002028 ** 
# Diagnosis                    0.5779  2   0.749052    
# Age                         20.5435  1  5.829e-06 ***
# Gender                       0.0051  1   0.943254    
# pH                           5.4789  1   0.019247 *  
# Hours.Final                 15.1954  1  9.694e-05 ***
# TZP_BioAnalyzer_RIN          5.3279  1   0.020986 *  
# TZP_Average.RNAConc..ng.uL.  3.7556  1   0.052632 .  
# CardBlock2                  75.0648  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           1.4088  1    0.23526    
# Diagnosis                             0.6562  2    0.72030    
# Age                                  21.8398  1  2.964e-06 ***
# Gender                                0.0236  1    0.87792    
# pH                                    5.4936  1    0.01909 *  
# Hours.Final                          15.8511  1  6.853e-05 ***
# TZP_BioAnalyzer_RIN                   4.6002  1    0.03197 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.1437  1    0.70465    
# TZP_Average.RNAConc..ng.uL.           3.1389  1    0.07645 .  
# TZP_Average.260.230                   2.3157  1    0.12808    
# TZP_Average.260.280                   0.0825  1    0.77392    
# CardBlock2                           75.0465  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 12.2538  1  0.0004643 ***
# Diagnosis                    1.1819  2  0.5538139    
# Age                          2.1174  1  0.1456321    
# Gender                       0.0013  1  0.9714122    
# pH                           8.4669  1  0.0036166 ** 
# Hours.Final                  2.7202  1  0.0990836 .  
# TZP_BioAnalyzer_RIN          4.4748  1  0.0343991 *  
# TZP_Average.RNAConc..ng.uL.  2.3065  1  0.1288344    
# CardBlock2                   7.1066  1  0.0076802 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)                           4.2573  1   0.039081 * 
# Diagnosis                             0.7291  2   0.694515   
# Age                                   1.3591  1   0.243697   
# Gender                                0.0285  1   0.865852   
# pH                                   10.5935  1   0.001135 **
# Hours.Final                           3.2349  1   0.072085 . 
# TZP_BioAnalyzer_RIN                   5.8537  1   0.015544 * 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  1.3708  1   0.241683   
# TZP_Average.RNAConc..ng.uL.           0.3428  1   0.558236   
# TZP_Average.260.230                   0.0594  1   0.807374   
# TZP_Average.260.280                   0.7255  1   0.394336   
# CardBlock2                            7.0957  1   0.007727 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   3.4428  1    0.06353 .  
# Diagnosis                     0.9562  2    0.61997    
# Age                          21.8020  1  3.023e-06 ***
# Gender                        0.1690  1    0.68097    
# pH                            5.3388  1    0.02086 *  
# Hours.Final                   2.4906  1    0.11453    
# TZP_BioAnalyzer_RIN           0.6454  1    0.42176    
# TZP_Average.RNAConc..ng.uL.   0.1085  1    0.74185    
# CardBlock2                  215.3465  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            2.8395  1  0.0919720 .  
# Diagnosis                              3.4620  2  0.1771114    
# Age                                   18.8014  1  1.451e-05 ***
# Gender                                 1.1351  1  0.2866919    
# pH                                     1.2999  1  0.2542363    
# Hours.Final                            1.9300  1  0.1647565    
# TZP_BioAnalyzer_RIN                    6.4316  1  0.0112109 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  13.2352  1  0.0002747 ***
# TZP_Average.RNAConc..ng.uL.            4.7345  1  0.0295642 *  
# TZP_Average.260.230                    0.6282  1  0.4280125    
# TZP_Average.260.280                    5.2454  1  0.0220045 *  
# CardBlock2                           215.3483  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#So going back to this question:
#I've included our normal biological co-variates. As expected, Age, pH, and PMI all seem to matter quite a bit for PC1-PC3 in this cleaned/normalized dataset (unlike for the noisier Cq). 
#Gender *never* matters (which is true in our microarray and RNA-Seq data too!)...but it is also the only variable that is significantly correlated with diagnosis and therefore seems important to include to prevent confound... and reviewers always ask about it. :(
#RIN still matters, although not as much as before.
#RNA concentration no longer seems to matter as much - I'm guessing the normalization using average housekeeping gene took care of these effects in the dataset (as one might hope they would...). Remove from model? Still matters a little for PC3, but maybe just due to a few outlier subjects.
#Also other more esoteric RNA measurements seem to potentially matter for PC3: TZP_BioAnalyzer_rRNA.Ratio..28s.18s., TZP_Average.260.280. Interestingly, RNA concentration only seems to matter when these other co-variates are included.
#Including all of those with Card is a pretty unwieldy model. With cardblock, maybe doable.

#Trying stepwise to see if that helps increase my confidence regarding which variables to keep/toss from the model:


library(lmerTest)
#This is backward elimination for mixed-effects models.

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+CardBlock2+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)

#using more severe alpha because I'm worried about overfitting:

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
# y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + CardBlock2 + (1 | ID)


Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
# Model found:
#   y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_Average.RNAConc..ng.uL. + TZP_Average.260.230 + Card + (1 | ID)
#hmm... interesting. A little different, but not much - contains purity when card is in the model instead of dissection group.

#With more severe alpha:
step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))
# Model found:
#   y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_Average.RNAConc..ng.uL. + TZP_Average.260.230 + Card + (1 | ID)



Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+CardBlock2+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_Average.RNAConc..ng.uL. + Dissecton.Group + CardBlock2 + (1 | ID)

#with more severe alpha:
step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))
# Model found:
#   y ~ Diagnosis + Gender + Hours.Final + TZP_BioAnalyzer_RIN + TZP_Average.RNAConc..ng.uL. + Dissecton.Group + CardBlock2 + (1 | ID)


Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_Average.260.230 + TZP_Average.260.280 + Card + (1 | ID)
#Interesting - RNAconc is replaced with Purity when Card is in the model (and not dissection group)

#More severe alpha:
step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))
# Model found:
#   y ~ Diagnosis + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_Average.260.230 + TZP_Average.260.280 + Card + (1 | ID)


Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+CardBlock2+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Age + Gender + TZP_BioAnalyzer_RIN + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + TZP_Average.RNAConc..ng.uL. + TZP_Average.260.280 + CardBlock2 + (1 | ID)

#Interesting - purity matters for PC3

#More severe alpha:
step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))
# Model found:
#   y ~ Diagnosis + Age + Gender + TZP_BioAnalyzer_RIN + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + TZP_Average.RNAConc..ng.uL. + TZP_Average.260.280 + CardBlock2 + (1 | ID)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Age + Gender + TZP_BioAnalyzer_RIN + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + TZP_Average.RNAConc..ng.uL. + TZP_Average.260.280 + Card + (1 | ID)

#more severe alpha:
step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))
# y ~ Diagnosis + Age + Gender + TZP_BioAnalyzer_RIN + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + TZP_Average.RNAConc..ng.uL. + TZP_Average.260.280 + Card + (1 | ID)

#Interesting - purity matters for PC3


#...but I'm still nervous about the amount of Df used up by some of these variables, and there hasn't been a direct comparison of models with Card to models with dissection group (w/ cardblock2?)

#PC1 - which model fits best?
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ Dissecton.Group + (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ Dissecton.Group + CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# Model  11 788.08 819.87 -383.04   766.08                            
# Model1 46 738.88 871.84 -323.44   646.88 119.2     35  4.095e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The model with card has better AIC but worse BIC
#According to ANOVA, card is still a better model fit.

anova(Model, Model2)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model2: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model2:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 788.08 819.87 -383.04   766.08                         
# Model2 23 802.12 868.60 -378.06   756.12 9.9564     12     0.6198

#Dissection group by itself is not an improvement.

anova(Model, Model3)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 788.08 819.87 -383.04   766.08                             
# Model3 12 738.59 773.28 -357.30   714.59 51.488      1  7.203e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#CardBlock2 -  better than no technical co-variate on all counts.
#And better in terms of both AIC and BIC.

anova(Model3, Model1)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model3 12 738.59 773.28 -357.30   714.59                             
# Model1 46 738.88 871.84 -323.44   646.88 67.711     34  0.0005129 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#CardBlock2 has better AIC and BIC... but ANOVA says that Card is still a better model???
#That is confusing.

anova(Model3, Model4)

# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + CardBlock2 + 
#   Model4:     (1 | ID)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# Model3 12 738.59 773.28 -357.30   714.59                        
# Model4 24 752.95 822.32 -352.47   704.95 9.645     12     0.6471

#Nope - adding dissection group is still not an improvement.



#Let's take a peek at PC2 also: 

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ Dissecton.Group + (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ Dissecton.Group + CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
# Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
# Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
# Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 564.22 596.02 -271.11   542.22                             
# Model1 46 544.27 677.22 -226.13   452.27 89.956     35   9.97e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#AIC is better with Card & BIC is worse with Card 
#... but Anova seems to suggest using card is better than no technical co-variate. (?!)

anova(Model, Model3)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Model  11 564.22 596.02 -271.11   542.22                            
# Model3 12 559.43 594.12 -267.72   535.43 6.7903      1   0.009165 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#card block 2 is an improvement over no covariates.

anova(Model3, Model1)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model3 12 559.43 594.12 -267.72   535.43                             
# Model1 46 544.27 677.22 -226.13   452.27 83.166     34  5.301e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Card has better AIC and worse BIC than CardBlock2
#But the anova still suggests a significant improvement??

anova(Model, Model2)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model2: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model2:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Model  11 564.22 596.02 -271.11   542.22                           
# Model2 23 566.93 633.41 -260.46   520.93 21.297     12     0.0462 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Dissection group is worse for AIC and BIC, and a very slight improvement with ANOVA.

anova(Model3, Model4)

# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + CardBlock2 + 
#   Model4:     (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Model3 12 559.43 594.12 -267.72   535.43                           
# Model4 24 562.08 631.45 -257.04   514.08 21.351     12    0.04547 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Similarly, when cardblock2 is in the model, dissection is a very slight improvement with ANOVA, but worse for both AIC and BIC.


#One more check: PC3

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ Dissecton.Group + (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ Dissecton.Group + CardBlock2+ (1 | ID), data = Temp, REML=F)


anova(Model, Model1)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 654.40 686.20 -316.20   632.40                             
# Model1 46 585.48 718.43 -246.74   493.48 138.92     35  2.538e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Card has better AIC but worse BIC than no technical covariates
#Anova indicates Card is much better fit.

anova(Model, Model3)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 654.40 686.20 -316.20   632.40                             
# Model3 12 558.07 592.75 -267.03   534.07 98.337      1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Card block2 is a major improvement on all counts than a model with no technical co-variates

anova(Model3, Model1)

# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model3 12 558.07 592.75 -267.03   534.07                         
# Model1 46 585.48 718.43 -246.74   493.48 40.588     34     0.2027

#Card is not significantly better than CardBlock2, and has better AIC and BIC.

#So the most efficient model is just odd/even cards for PC3

anova(Model, Model2)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model2: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model2:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 654.40 686.20 -316.20   632.40                         
# Model2 23 669.54 736.02 -311.77   623.54 8.8599     12     0.7149

#Dissection group is not an improvement.

anova(Model3, Model4)

# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + CardBlock2 + 
#   Model4:     (1 | ID)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# Model3 12 558.07 592.75 -267.03   534.07                        
# Model4 24 570.61 639.97 -261.30   522.61 11.46     12     0.4899

#Dissection group is not an improvement after controlling for cardblock2 either.

#So it seems like a toss up regarding which choice of co-variates is best for the data as a whole.
#Also a toss-up: I've included our normal biological co-variates, but some don't seem to matter with this dataset nearly as much as the technical co-variates. Toss out?  And which technical co-variates should we include - there were some non-traditional ones that seemed related to the PCs in this cleaned up analysis (purity, 28s/18s)



#*********

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+Card+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+ CardBlock2+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     Card + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Model  45 816.99 948.06 -363.50   726.99                           
# Model1 46 816.14 950.13 -362.07   724.14 2.8461      1    0.09159 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(Model3, Model4)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     CardBlock2 + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model3 11 806.64 838.68 -392.32   784.64                         
# Model4 12 807.00 841.95 -391.50   783.00 1.6385      1     0.2005

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+Card+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+ CardBlock2+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
# Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
# Model:     Card + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
# Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  45 600.34 731.41 -255.17   510.34                         
# Model1 46 600.85 734.83 -254.43   508.85 1.4911      1     0.2221

anova(Model3, Model4)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     CardBlock2 + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model3 11 595.01 627.04 -286.50   573.01                         
# Model4 12 596.53 631.48 -286.26   572.53 0.4787      1      0.489

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+Card+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+ CardBlock2+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     Card + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  45 572.76 703.83 -241.38   482.76                         
# Model1 46 572.37 706.35 -240.19   480.37 2.3908      1      0.122

anova(Model3, Model4)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     CardBlock2 + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model3 11 545.40 577.44 -261.7   523.40                         
# Model4 12 546.41 581.36 -261.2   522.41 0.9924      1     0.3192


#Hmmm....  I'm tempted to throw RNAConcentration out given the small sample size.
#On the other hand, if we want to compare results before and after -DeltaCq normalization, it might help to keep it in the model to make things equal.





#**************************