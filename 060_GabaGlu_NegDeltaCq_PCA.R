#Code for calculating NegDeltaCq & then doing QC/PCA again:

#########################################################

#Previously there was an additional qc step here (now placed earlier), so now this section of code simply exists to make sure that the objects are renamed so that I can reuse later code.

dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#[1] 136  48

dim(GabaGlu_Cq_AllSubjects_QCed2)
#[1]   93 136

length(GabaGlu_MeanHousekeeping_QCed)
#[1] 136

dim(GabaGlu_Cq_AllSubjects_QCed3_Imputed)
#[1]  93 136

GabaGlu_Cq_AllSubjects_QCed3<-GabaGlu_Cq_AllSubjects_QCed2
dim(GabaGlu_Cq_AllSubjects_QCed3)
#[1]  93 136


#####################


#Calculating negative delta Cq:
str(GabaGlu_Cq_AllSubjects_QCed3)
# num [1:93, 1:136] 20.5 25.7 23.7 27.3 22.8 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:93] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# ..$ : chr [1:136] "1" "2" "3" "4" ...
head(GabaGlu_Cq_AllSubjects_QCed3)
# 1      2      3      4     17     18     19     20     21     22     23     24     21     22     23     24     25     26
# ABAT    20.513 20.655 20.957 20.553 21.166 21.669 20.939 21.598 21.161 20.989 21.125 20.875 20.683 22.012 21.678 20.426 21.716 21.090
# ADCY7   25.734 24.956 25.178 25.612 26.201 26.061 25.861 26.895 25.576 24.902 25.002 25.145 25.469 25.750 25.353 25.210 26.383 25.021
# ADORA1  23.737 23.150 23.728 23.443 23.945 24.597 24.181 24.652 24.189 23.408 23.796 23.617 24.012 24.702 24.492 23.849 24.586 23.241
# ADORA2A 27.270 26.826 26.960 26.980 27.683 27.646 27.191 28.035 27.333 26.829 26.842 26.795 27.148 28.093 27.573 26.906 27.727 26.831

#Sanity check:
GabaGlu_Cq_AllSubjects_QCed3[1,]
#20.513 20.655 20.957 20.553...
GabaGlu_MeanHousekeeping_QCed

GabaGlu_Cq_AllSubjects_QCed3[1,]-GabaGlu_MeanHousekeeping_QCed


GabaGlu_NegDeltaCq_AllSubjects_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed3, 1, function(y) ((y-GabaGlu_MeanHousekeeping_QCed)*(-1)))
str(GabaGlu_NegDeltaCq_AllSubjects_QCed)

# num [1:138, 1:93] 0.934 1.447 1.575 1.674 1.626 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:138] "1" "2" "3" "4" ...
# ..$ : chr [1:93] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...

head(GabaGlu_NegDeltaCq_AllSubjects_QCed)
#The apply function transposed things. Huh. I wonder why.

sum(is.na(GabaGlu_Cq_AllSubjects_QCed3))
#[1] 19

pdf("GabaGlu_Boxplot_Cq_perSample_AllGenes_QCed3.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_QCed3), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq per sample for all genes", xlab="Sample ID", ylab="Cq")
dev.off()

sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_QCed))
#[1] 19

pdf("GabaGlu_Boxplot_NegDeltaCq_perSample_AllGenes_QCed.pdf", width=20, height=5)
boxplot(data.frame(t(GabaGlu_NegDeltaCq_AllSubjects_QCed)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq per sample for all genes", xlab="Sample ID", ylab="-DeltaCq")
dev.off()
#That seems to have reduced the variability between most samples and increased the variability (or perhaps the ability to detect the variability?) between a few others.

#***************************

#Scaling -DeltaCq and Running PCA again:

#While scaling I'm going to also change the orientation back to match my previous code:
GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed<-t(scale(GabaGlu_NegDeltaCq_AllSubjects_QCed, center=T, scale=T))
head(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
#Sanity check:
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, 1, function(y) mean(y, na.rm=T))
apply(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, 1, function(y) sd(y, na.rm=T))
#Looks good.

GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed<-cor(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")

boxplot(GabaGlu_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed)

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

#Huh. It centered the distribution near zero for most samples, but then shifted a few samples away from 0 (x1, x37 - notably x1.1 and x37.1 look normal)

pdf("Boxplot_NegDeltaCqZscore_perSample_HousekeepingGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed[83:93,]), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq Z-scores per sample for housekeeping genes", xlab="Sample ID", ylab="-DeltaCq Z-score")
dev.off()
#Yeah, the subjects with the weirdest distributions for all genes also have super variable housekeeping gene distributions (Sample 1 and 37. 25 also has a variable housekeeping gene distribution, but doesn't have values that seem weird for all genes.)


#Trying PCA:

pca_GabaGlu<-prcomp(t(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed))
tmp<-pca_GabaGlu$x[,1:10]
dim(tmp)
#[1]  136  10
rownames(tmp)<-colnames(GabaGlu_NegDeltaCq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_GabaGlu_NegaDeltaCq.csv")

tmp<-pca_GabaGlu$rotation[,1:10]
write.csv(tmp, "pca_GabaGlu_NegDeltaCq_Eigenvectors.csv")

png("PCA_ScreePlot_GabaGlu_NegDeltaCq.png")
plot(summary(pca_GabaGlu)$importance[2,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 explains a smaller portion of variance now (35%)

png("PCA_ScreePlot2_GabaGlu_NegDeltaCq.png")
plot(summary(pca_GabaGlu)$importance[3,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()

png("PC1vsPC2_GabaGlu_NegDeltaCq_byDiagnosis.png")
plot(pca_GabaGlu$x[,1]~pca_GabaGlu$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis))
dev.off()
#PC1 Not diagnosis
#Although there are several outliers that are all from the same diagnosis group (controls?)
hist(pca_GabaGlu$x[,1], breaks=20)
pca_GabaGlu$x[pca_GabaGlu$x[,1]<(-15),1]
# 1        37        16 
# -19.95283 -32.47732 -21.75890 
pca_GabaGlu$x[pca_GabaGlu$x[,1]<(-8),1]
#the replicate samples are not nearby.

pdf("PC1_GabaGlu_NegDeltaCq_byMeanHousekeeping.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_MeanHousekeeping_QCed)
dev.off()
#The general correlation with mean housekeeping gene Cq has gone away, but the three samples with PC1<(-15) are all also samples with lower Mean Housekeeping Cq (higher expression levels)
#Sample 1 and 37 have been flagged before for other reasons.
#All 3 samples had relatively low replicate sample-sample correlations in the Cq data prior to QC (<0.975)

colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Temp<-cbind(pca_GabaGlu$x[,c(1:4)],GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[,c(7,9,11,16,24:28,38:42)] )
cor(Temp)
write.csv(cor(Temp), "GabaGlu_CorMatrix_PCA_GabaGlu_NegDeltaCq_VsSubjectVar.csv")

#Potentially worth exploring:
#PC1 vs. GabaGlu_MeanHousekeeping_QCed, Age, pH, PMI, TZP_BioAnalyzer_rRNA.Ratio..28s.18s.,TZP_Average.260.280
#PC2 vs. GabaGlu_MeanHousekeeping_QCed, pH, PMI, RIN, TZP_BioAnalyzer_rRNA.Ratio..28s.18s.,TZP_Average.260.280
#PC3 vs. GabaGlu_MeanHousekeeping_QCed, Age, pH,  TZP_BioAnalyzer_rRNA.Ratio..28s.18s.


#I came back and looped a statistical output for all of these

#Technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)                          1.6512  1    0.19879  
# Diagnosis                            3.6188  2    0.16376  
# Block.Weight..g.                     1.8130  1    0.17815  
# TZP_BioAnalyzer_RIN                  0.1871  1    0.66534  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 2.9412  1    0.08635 .
# TZP_Average.RNAConc..ng.uL.          0.0645  1    0.79953  
# TZP_Average.260.230                  0.2395  1    0.62459  
# TZP_Average.260.280                  1.7581  1    0.18486  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Basically no relationship.

#Biological variables:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)  6.1618  1   0.013054 *  
# Diagnosis    3.2734  2   0.194620    
# Age         19.7408  1  8.869e-06 ***
# Gender       0.2991  1   0.584434    
# pH           6.9062  1   0.008590 ** 
# Hours.Final  9.1305  1   0.002514 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Biological and technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           0.0004  1   0.984972    
# Diagnosis                             2.9263  2   0.231501    
# Age                                  19.7744  1  8.714e-06 ***
#   Gender                                0.6982  1   0.403383    
# pH                                    4.9556  1   0.026006 *  
#   Hours.Final                          11.3834  1   0.000741 ***
#   Block.Weight..g.                      0.4395  1   0.507359    
# TZP_BioAnalyzer_RIN                   0.9602  1   0.327137    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.0259  1   0.872030    
# TZP_Average.RNAConc..ng.uL.           0.3392  1   0.560290    
# TZP_Average.260.230                   1.8776  1   0.170604    
# TZP_Average.260.280                   0.8608  1   0.353525    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - RNA metrics less after controlling for HK expression, but other biological variables matter more.

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                  13.1033  1  0.0002948 ***
# Diagnosis                     8.3063  2  0.0157145 *  
# Age                          32.3583  1  1.282e-08 ***
# Gender                        0.2613  1  0.6092096    
# pH                           10.5384  1  0.0011692 ** 
# Hours.Final                   9.6199  1  0.0019248 ** 
# TZP_BioAnalyzer_RIN           0.7934  1  0.3730726    
# TZP_Average.RNAConc..ng.uL.   2.9013  1  0.0885096 .  
# Card                        136.3578 35  6.771e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - RIN and RNAconc matter less after controlling for HK expression, but other biological variables matter more.

#Biological and more technical variables & card
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            0.3459  1   0.556416    
# Diagnosis                              4.1000  2   0.128733    
# Age                                   30.2507  1  3.797e-08 ***
# Gender                                 0.0144  1   0.904446    
# pH                                     6.8754  1   0.008739 ** 
# Hours.Final                            7.4725  1   0.006265 ** 
# TZP_BioAnalyzer_RIN                    1.2727  1   0.259264    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   0.2001  1   0.654606    
# TZP_Average.RNAConc..ng.uL.            0.0203  1   0.886592    
# TZP_Average.260.230                    6.7231  1   0.009517 ** 
# TZP_Average.260.280                    2.7403  1   0.097847 .  
# Card                                 141.3161 35  1.012e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - purity seems to matter now (???).

#With dissection group instead of card:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           0.5965  1   0.439912    
# Diagnosis                             2.5821  2   0.274979    
# Age                                  22.3879  1  2.228e-06 ***
# Gender                                0.7755  1   0.378515    
# pH                                    4.5380  1   0.033150 *  
# Hours.Final                           5.6529  1   0.017427 *  
# TZP_BioAnalyzer_RIN                   2.0316  1   0.154056    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.0622  1   0.803077    
# TZP_Average.RNAConc..ng.uL.           0.0008  1   0.977735    
# TZP_Average.260.230                   8.6818  1   0.003214 ** 
# TZP_Average.260.280                   2.6970  1   0.100536    
# Dissecton.Group                      15.7131 12   0.204733    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Again, Dissection group doesn't seem to matter, but purity might.



#PC1 varies with Age, pH, PMI and Card, as well as diagnosis. It no longer varies with RNA concentration, RIN, and dissection - maybe due to normalizing with the housekeeping genes?  

png("PC1_GabaGlu_NegDeltaCq_byDiagnosis.png")
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis, ylab="PC1")
dev.off()
#pretty unconvincing
#Notably, those outlier samples are outliers within the control distribution too

pdf("PC1_GabaGlu_NegDeltaCq_byDissection.pdf", width=10, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Hours.Final)
dev.off()

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Hours.Final + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept) 3.9229  1    0.04763 *
#   Hours.Final 4.3397  1    0.03723 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("PC1_GabaGlu_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$pH)
dev.off()

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~pH + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept) 6.1899  1    0.01285 *
#   pH          6.1921  1    0.01283 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("PC1_GabaGlu_NegDeltaCq_byAge.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Age)
dev.off()
#Maybe Age...the outliers are interfering though.

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Age + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 19.773  1  8.722e-06 ***
#   Age         21.290  1  3.948e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf("PC1_GabaGlu_NegDeltaCq_byGender.pdf", width=5, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Gender)
dev.off()
#The outliers are all male, but otherwise the distributions are noticeably different.


#I'm going to check on some other metrics that didn't matter when looking at Cq:
pdf("PC1_GabaGlu_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC1_GabaGlu_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC1_GabaGlu_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.230)
dev.off()
#Contains big outliers - not so helpful.

pdf("PC1_GabaGlu_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()
#Contains big outliers, not so helpful.

pdf("PC1_GabaGlu_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
dev.off()

pdf("PC1_GabaGlu_NegDeltaCq_byCard.pdf", width=20, height=5)
boxplot(pca_GabaGlu$x[,1]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, las=3, par(cex.axis=0.75))
dev.off()
#Card is the most convincing so far, but that may just be because there are so few samples on a card the distributions can be tight.
#Although replicate cards sometimes do not look so similar...


Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Card + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)   3.9816  1      0.046 *  
#   Card        127.6107 35  1.837e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#hmmm... controlling for card explicitly would use up a whole bunch of degrees of freedom. :(
#It also seems somewhat redundant with just controlling for ID.

#But it does seem to clean things up a lot...

#I'm nervous about including it though, because there are so few samples per card that making the card effect essentially 0 seems like it would be removing a lot of the subject-level variability, although I suppose if we're simultaneously controlling for diagnosis it is less of an issue?
#Alternatively, we could make the gene values average out to the mean for both replicate cards.


#PC2:


#Technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)                          0.2929  1    0.58839  
# Diagnosis                            3.3550  2    0.18684  
# Block.Weight..g.                     0.1741  1    0.67646  
# TZP_BioAnalyzer_RIN                  4.1750  1    0.04103 *
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 0.0069  1    0.93365  
# TZP_Average.RNAConc..ng.uL.          0.5570  1    0.45548  
# TZP_Average.260.230                  0.2100  1    0.64680  
# TZP_Average.260.280                  0.1640  1    0.68551  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#RIN may matter a little bit for PC2

#Biological variables:
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 11.5642  1  0.0006724 ***
#   Diagnosis    3.0561  2  0.2169558    
# Age          0.9667  1  0.3254968    
# Gender       0.6536  1  0.4188143    
# pH          10.9409  1  0.0009407 ***
#   Hours.Final  2.9473  1  0.0860192 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#pH is strongly related to PC2, maybe PMI

#Biological and technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           6.6648  1  0.0098336 ** 
# Diagnosis                             1.7628  2  0.4141991    
# Age                                   0.0187  1  0.8912514    
# Gender                                0.5025  1  0.4784000    
# pH                                   11.4261  1  0.0007242 ***
# Hours.Final                           6.5249  1  0.0106374 *  
# Block.Weight..g.                      1.3786  1  0.2403431    
# TZP_BioAnalyzer_RIN                   8.7870  1  0.0030338 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  2.0016  1  0.1571313    
# TZP_Average.RNAConc..ng.uL.           0.0067  1  0.9349673    
# TZP_Average.260.230                   0.3872  1  0.5337597    
# TZP_Average.260.280                   1.7621  1  0.1843666    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 11.8942  1  0.0005631 ***
# Diagnosis                    0.6546  2  0.7208566    
# Age                          1.4671  1  0.2258068    
# Gender                       0.5665  1  0.4516333    
# pH                           6.8448  1  0.0088900 ** 
# Hours.Final                  3.7738  1  0.0520605 .  
# TZP_BioAnalyzer_RIN          9.3493  1  0.0022307 ** 
# TZP_Average.RNAConc..ng.uL.  1.5078  1  0.2194759    
# Card                        85.1012 35  4.661e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#PC2 is related to RIN, pH, PMI, CARD

#Biological and more technical variables & card
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           7.0702  1   0.007838 ** 
# Diagnosis                             0.4511  2   0.798071    
# Age                                   1.2019  1   0.272947    
# Gender                                0.2352  1   0.627679    
# pH                                    8.0487  1   0.004554 ** 
# Hours.Final                           4.7458  1   0.029369 *  
# TZP_BioAnalyzer_RIN                   5.2243  1   0.022274 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.5717  1   0.449603    
# TZP_Average.RNAConc..ng.uL.           0.0123  1   0.911654    
# TZP_Average.260.230                   2.7850  1   0.095153 .  
# TZP_Average.260.280                   2.5732  1   0.108687    
# Card                                 87.4409 35   2.23e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Maybe purity again???

#With dissection group instead of card:
Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#Dissection group doesn't matter.


png("PC2_GabaGlu_NegDeltaCq_byDiagnosis.png")
boxplot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis, ylab="PC2")
dev.off()
#Nope

pdf("PC2_GabaGlu_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$pH)
dev.off()

pdf("PC2_GabaGlu_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Hours.Final)
dev.off()

pdf("PC2_GabaGlu_NegDeltaCq_byRIN.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_RIN)
dev.off()

#Double-checking some of the other metrics that didn't matter for Cq:

pdf("PC2_GabaGlu_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC2_GabaGlu_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC2_GabaGlu_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.230)
dev.off()

pdf("PC2_GabaGlu_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,2]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()

#Technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           5.1387  1    0.02340 *  
# Diagnosis                             0.6539  2    0.72112    
# Block.Weight..g.                      0.8681  1    0.35149    
# TZP_BioAnalyzer_RIN                   6.5784  1    0.01032 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 16.4114  1  5.098e-05 ***
# TZP_Average.RNAConc..ng.uL.           6.1631  1    0.01304 *  
# TZP_Average.260.230                   1.9197  1    0.16589    
# TZP_Average.260.280                   5.2131  1    0.02242 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Huh. That's interesting - strong relationship with TZP_BioAnalyzer_rRNA.Ratio..28s.18s, slight relationship with purity.

#Biological variables:
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)  3.3950  1    0.06540 .  
# Diagnosis    0.6060  2    0.73860    
# Age         15.9636  1  6.457e-05 ***
# Gender       0.1141  1    0.73548    
# pH           4.4277  1    0.03536 *  
# Hours.Final  2.5270  1    0.11191    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Strong relationship with age.

#Biological and technical variables
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)                          2.6533  1   0.103335   
# Diagnosis                            1.8696  2   0.392673   
# Age                                  9.5446  1   0.002005 **
# Gender                               1.2192  1   0.269527   
# pH                                   1.0690  1   0.301162   
# Hours.Final                          0.7252  1   0.394434   
# Block.Weight..g.                     0.6028  1   0.437506   
# TZP_BioAnalyzer_RIN                  3.7115  1   0.054038 . 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 8.1415  1   0.004326 **
# TZP_Average.RNAConc..ng.uL.          4.7105  1   0.029980 * 
# TZP_Average.260.230                  0.8227  1   0.364384   
# TZP_Average.260.280                  4.6764  1   0.030580 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   1.7234  1    0.18926    
# Diagnosis                     0.9950  2    0.60804    
# Age                          24.1694  1  8.822e-07 ***
# Gender                        0.2970  1    0.58578    
# pH                            5.0966  1    0.02397 *  
# Hours.Final                   2.7092  1    0.09977 .  
# TZP_BioAnalyzer_RIN           0.2265  1    0.63411    
# TZP_Average.RNAConc..ng.uL.   2.4367  1    0.11853    
# Card                        373.9835 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#No relationship with RIN by itself.

#Biological and more technical variables & card
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            3.7167  1  0.0538712 .  
# Diagnosis                              3.2430  2  0.1976067    
# Age                                   14.5167  1  0.0001389 ***
# Gender                                 0.0297  1  0.8631631    
# pH                                     0.6677  1  0.4138532    
# Hours.Final                            0.3400  1  0.5598579    
# TZP_BioAnalyzer_RIN                    5.2268  1  0.0222408 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  13.8372  1  0.0001993 ***
# TZP_Average.RNAConc..ng.uL.           10.1907  1  0.0014115 ** 
# TZP_Average.260.230                    0.0100  1  0.9205366    
# TZP_Average.260.280                    5.2639  1  0.0217719 *  
# Card                                 371.0424 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#With dissection group instead of card:
Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#Dissection group definitely still doesn't matter.


png("PC3_GabaGlu_NegDeltaCq_byDiagnosis.png")
boxplot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis, ylab="PC2")
dev.off()
#Nope

pdf("PC3_GabaGlu_NegDeltaCq_byAge.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Age)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$pH)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Hours.Final)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byRIN.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_RIN)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
dev.off()
#The relationship with RNAconc seems driven by a few subjects with low concentration

#Double-checking some of the other metrics that didn't matter for Cq:

pdf("PC3_GabaGlu_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC3_GabaGlu_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC3_GabaGlu_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.230)
dev.off()

pdf("PC3_GabaGlu_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,3]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()
#Reasonably convincing


#Is there some way to collapse card? (so that it doesn't take up so many df in the model???)

dim(pca_GabaGlu$x)

SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$DissectionGroupCard<-paste(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group,SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card)

pca_GabaGlu_byDissectionCard<-matrix(0, 52, 10)
row.names(pca_GabaGlu_byDissectionCard)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$DissectionGroupCard))

for(i in c(1:10)){
  pca_GabaGlu_byDissectionCard[,i]<-tapply(X=pca_GabaGlu$x[,i], INDEX=SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$DissectionGroupCard, FUN=mean)
}

#There are some very strong blocks - let's take a peek
pca_GabaGlu_byDissectionCard_CorMatrix<-cor(t(pca_GabaGlu_byDissectionCard))
dim(pca_GabaGlu_byDissectionCard_CorMatrix)
row.names(pca_GabaGlu_byDissectionCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$DissectionGroupCard))
colnames(pca_GabaGlu_byDissectionCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$DissectionGroupCard))

pdf("Heatmap_pca_GabaGlu_byDissectionCard_CorMatrix.pdf", width=10, height=10)
heatmap(pca_GabaGlu_byDissectionCard_CorMatrix)
dev.off()

pca_GabaGlu_byDissectionCard_CorMatrix

write.csv(pca_GabaGlu_byDissectionCard_CorMatrix, "pca_GabaGlu_byDissectionCard_NegDeltaCq_CorMatrix.csv")

pca_GabaGlu_byDissection<-matrix(0, 13, 10)
row.names(pca_GabaGlu_byDissection)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group))

for(i in c(1:10)){
  pca_GabaGlu_byDissection[,i]<-tapply(X=pca_GabaGlu$x[,i], INDEX=SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group, FUN=mean)
}

#There are some very strong blocks - let's take a peek
pca_GabaGlu_byDissection_CorMatrix<-cor(t(pca_GabaGlu_byDissection))
dim(pca_GabaGlu_byDissection_CorMatrix)
row.names(pca_GabaGlu_byDissection_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group))
colnames(pca_GabaGlu_byDissection_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group))

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
row.names(pca_GabaGlu_byCard)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))

for(i in c(1:10)){
  pca_GabaGlu_byCard[,i]<-tapply(X=pca_GabaGlu$x[,i], INDEX=SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, FUN=mean)
}


pca_GabaGlu_byCard_CorMatrix<-cor(t(pca_GabaGlu_byCard))
dim(pca_GabaGlu_byCard_CorMatrix)
row.names(pca_GabaGlu_byCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))
colnames(pca_GabaGlu_byCard_CorMatrix)<-names(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card))

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


SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$CardBlock2<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card

SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$CardBlock2[SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$CardBlock2%in%c("Card 24.eds", "Card 4.eds", "Card 2.eds", "Card 28.eds", "Card 30.eds", "Card 22.eds", "Card 26.eds", "Card 12.eds", "Card 18.eds", "Card 20.eds", "Card 34.eds", "Card 32.eds", "Card 36.eds", "Card 8.eds", "Card 14.eds", "Card 16.eds", "Card 10.eds", "Card 6.eds")]<-"EvenCards"

SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$CardBlock2[SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$CardBlock2%in%c("Card 27.eds", "Card 9.eds", "Card 13.eds", "Card 29.eds", "Card 31.eds", "Card 21.eds", "Card 1.eds", "Card 15.eds", "Card 23.eds", "Card 25.eds", "Card 7.eds", "Card 19.eds", "Card 5.eds", "Card 33.eds", "Card 35.eds", "Card 17.eds", "Card 11.eds", "Card 3.eds")]<-"OddCards"

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$CardBlock2)
# EvenCards  OddCards 
# 67        69 

#PC1 - which model fits best?
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

# Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + CardBlock+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 837.72 869.76 -407.86   815.72                             
# Model1 46 816.14 950.13 -362.07   724.14 91.581     35  5.891e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The model with card has better AIC but worse BIC
#According to ANOVA, card is still a better model fit.

anova(Model, Model3)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 837.72 869.76 -407.86   815.72                             
# Model3 12 807.00 841.95 -391.50   783.00 32.727      1   1.06e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Same with CardBlock2 - also better than no technical co-variate.
#And better in terms of both AIC and BIC.

anova(Model3, Model1)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Model3 12 807.00 841.95 -391.50   783.00                            
# Model1 46 816.14 950.13 -362.07   724.14 58.853     34   0.005137 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#CardBlock2 has better AIC and BIC... but ANOVA says that Card is still a better model???
#That is confusing.


#Nope.

#Let's take a peek at PC2 also: 

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Model  11 596.37 628.41 -287.18   574.37                            
# Model1 46 600.85 734.83 -254.43   508.85 65.515     35   0.001337 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#AIC & BIC is worse with Card 
#... but Anova seems to suggest using card is better than no technical co-variate. (?!)

anova(Model, Model3)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 596.37 628.41 -287.18   574.37                         
# Model3 12 596.53 631.48 -286.26   572.53 1.8401      1     0.1749

#card block 2 is not an improvement over no covariates.

anova(Model3, Model1)
# Data: Temp
# Models:
#   Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + CardBlock2 + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Model3 12 596.53 631.48 -286.26   572.53                            
# Model1 46 600.85 734.83 -254.43   508.85 63.675     34   0.001516 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Card has worse AIC and BIC than CardBlock2
#But the anova still suggests a significant improvement??


#One more check: PC3

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 648.41 680.45 -313.20   626.41                             
# Model1 46 572.37 706.35 -240.19   480.37 146.04     35  1.628e-15 ***
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
# Df    AIC    BIC logLik deviance Chisq Chi Df Pr(>Chisq)    
# Model  11 648.41 680.45 -313.2   626.41                            
# Model3 12 546.41 581.36 -261.2   522.41   104      1  < 2.2e-16 ***
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
# Model3 12 546.41 581.36 -261.20   522.41                         
# Model1 46 572.37 706.35 -240.19   480.37 42.036     34      0.162

#Card is not significantly better than CardBlock2

#So the most efficient model is just odd/even cards for PC3



#So it seems like a toss up regarding which choice of co-variates is best for the data as a whole.
#Also a toss-up: I've included our normal biological co-variates, but some don't seem to matter with this dataset nearly as much as the technical co-variates. Toss out?  And which technical co-variates should we include - there were some non-traditional ones that seemed related to the PCs in this cleaned up analysis (purity, 28s/18s)

#Let's look at the results:
Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                  8.4404  1   0.003670 ** 
# Diagnosis                    4.4888  2   0.105990    
# Age                         25.0907  1  5.469e-07 ***
# Gender                       0.1636  1   0.685849    
# pH                           6.8332  1   0.008948 ** 
# Hours.Final                 14.0187  1   0.000181 ***
# TZP_BioAnalyzer_RIN          1.4290  1   0.231933    
# TZP_Average.RNAConc..ng.uL.  1.6582  1   0.197851    
# CardBlock2                  39.3176  1  3.602e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           0.0000  1  0.9944193    
# Diagnosis                             3.4879  2  0.1748253    
# Age                                  27.2167  1  1.819e-07 ***
# Gender                                0.3754  1  0.5400590    
# pH                                    6.0932  1  0.0135705 *  
# Hours.Final                          14.8373  1  0.0001172 ***
# TZP_BioAnalyzer_RIN                   1.8681  1  0.1716888    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.4638  1  0.4958369    
# TZP_Average.RNAConc..ng.uL.           0.3609  1  0.5480165    
# TZP_Average.260.230                   2.7762  1  0.0956778 .  
# TZP_Average.260.280                   0.9177  1  0.3380850    
# CardBlock2                           39.2459  1  3.736e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 14.0047  1  0.0001824 ***
# Diagnosis                    2.2069  2  0.3317198    
# Age                          0.4167  1  0.5185913    
# Gender                       0.0005  1  0.9819305    
# pH                           8.7053  1  0.0031729 ** 
# Hours.Final                  5.0585  1  0.0245050 *  
# TZP_BioAnalyzer_RIN          6.6748  1  0.0097783 ** 
# TZP_Average.RNAConc..ng.uL.  0.4803  1  0.4882781    
# CardBlock2                   1.8700  1  0.1714751    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           6.3290  1  0.0118780 *  
# Diagnosis                             1.2300  2  0.5406345    
# Age                                   0.1485  1  0.6999851    
# Gender                                0.0385  1  0.8443989    
# pH                                   11.0691  1  0.0008778 ***
# Hours.Final                           5.8081  1  0.0159530 *  
# TZP_BioAnalyzer_RIN                   7.3705  1  0.0066301 ** 
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.9297  1  0.3349385    
# TZP_Average.RNAConc..ng.uL.           0.0882  1  0.7664872    
# TZP_Average.260.230                   0.1415  1  0.7067708    
# TZP_Average.260.280                   1.6551  1  0.1982694    
# CardBlock2                            1.8983  1  0.1682682    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   4.1699  1    0.04115 *  
# Diagnosis                     0.8071  2    0.66794    
# Age                          19.8695  1  8.291e-06 ***
# Gender                        0.1669  1    0.68287    
# pH                            6.1397  1    0.01322 *  
# Hours.Final                   1.8970  1    0.16842    
# TZP_BioAnalyzer_RIN           0.2441  1    0.62127    
# TZP_Average.RNAConc..ng.uL.   1.0001  1    0.31727    
# CardBlock2                  227.0003  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+ CardBlock2+ (1 | ID), data = Temp, REML=F)
car::Anova(Model3, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            3.9233  1  0.0476223 *  
# Diagnosis                              3.7015  2  0.1571174    
# Age                                   16.5965  1  4.624e-05 ***
# Gender                                 1.2271  1  0.2679676    
# pH                                     1.6482  1  0.1991990    
# Hours.Final                            1.2944  1  0.2552399    
# TZP_BioAnalyzer_RIN                    5.0005  1  0.0253400 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  13.1316  1  0.0002904 ***
# TZP_Average.RNAConc..ng.uL.            9.1816  1  0.0024446 ** 
# TZP_Average.260.230                    0.8022  1  0.3704242    
# TZP_Average.260.280                    7.4011  1  0.0065186 ** 
# CardBlock2                           228.3428  1  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#So going back to this question:
#I've included our normal biological co-variates. As expected, Age, pH, and PMI all seem to matter quite a bit for PC1-PC3 in this cleaned/normalized dataset (unlike for the noisier Cq). 
#Gender *never* matters (which is true in our microarray and RNA-Seq data too!)...but it is also the only variable that is significantly correlated with diagnosis and therefore seems important to include to prevent confound... and reviewers always ask about it. :(
#RIN still matters, although not as much as before.
#RNA concentration no longer seems to matter as much - I'm guessing the normalization using average housekeeping gene took care of these effects in the dataset (as one might hope they would...). Remove from model? Still matters a little for PC3, but maybe just due to a few outlier subjects.
#Also other more esoteric RNA measurements seem to potentially matter for PC3: TZP_BioAnalyzer_rRNA.Ratio..28s.18s., TZP_Average.260.280. Interestingly, RNA concentration only seems to matter when these other co-variates are included.
#Including all of those with Card is a pretty unwieldy model. With cardblock, maybe doable.


Temp<-data.frame(y=pca_GabaGlu$x[,1], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

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

Temp<-data.frame(y=pca_GabaGlu$x[,2], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

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

Temp<-data.frame(y=pca_GabaGlu$x[,3], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

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

#To be thorough, I came back and looped the basic statistical output for each variable in relationship to PC1-4:
colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)

GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats<-matrix(0, 20, 4)
row.names(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[,c(2,5:9,11,14,16,24:28,38:42,50)])
colnames(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats)<-c("PC1_ANOVA_Pval","PC2_ANOVA_Pval","PC3_ANOVA_Pval","PC4_ANOVA_Pval")

Temp<-data.frame(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[,c(2,5:9,11,14,16,24:28,38:42,50)], pca_GabaGlu$x[,c(1:4)], ID=SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[,1])

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
 
setwd("~/Documents/Microarray Gen/FrontalPole/Output/GabaGlu")
write.csv(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats, "GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats.csv")

#Potentially worth exploring:
#PC1 vs. Age, pH, PMI, Card/Cardblock2, Trends: TZP_BioAnalyzer_rRNA.Ratio..28s.18s.,TZP_Average.260.280
#PC2 vs. pH, PMI, RIN, Card, Cohort(?), Trends: PMI, Dissection Group, TTP_Average260/280
#PC3 vs. Age, pH, 28s/18s, Card/Cardblock2, Trends: Cohort(?) 


#What about with the three outliers removed?

GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers<-matrix(0, 20, 4)
row.names(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[,c(2,5:9,11,14,16,24:28,38:42,50)])
colnames(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers)<-c("PC1_ANOVA_Pval","PC2_ANOVA_Pval","PC3_ANOVA_Pval","PC4_ANOVA_Pval")

Temp2<-Temp[Temp$PC1>(-15),]

for(i in c(1:20)){
  
  Model<-lmer(Temp2$PC1~Temp2[,i]+ (1 | Temp2$ID), data = Temp2, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers[i,1]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp2$PC2~Temp2[,i]+ (1 | Temp2$ID), data = Temp2, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers[i,2]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp2$PC3~Temp2[,i]+ (1 | Temp2$ID), data = Temp2, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers[i,3]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp2$PC4~Temp2[,i]+ (1 | Temp2$ID), data = Temp2, REML=F)
  GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers[i,4]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
}

setwd("~/Documents/Microarray Gen/FrontalPole/Output/GabaGlu")
write.csv(GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers, "GabaGlu_NegDeltaCq_PCA_vsSubjVariables_Stats_NoOutliers.csv")
#basically the same results, although it occurred to me that this analysis might not really make sense since PCA is defined by the main sources of variation in the data - really, if we wanted to go farther with this question we would need to re-run the PCA itself.
#... Sigh. Well, let's do it.

#**************************