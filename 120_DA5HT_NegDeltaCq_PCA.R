#*******************************


#Calculating negative delta Cq:
str(DA5HT_Cq_AllSubjects_QCed3)
# num [1:46, 1:154] 24.2 25.1 24.6 27.8 32.8 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:46] "ADRB1" "ADRB2" "COMT" "DBH" ...
# ..$ : chr [1:154] "1" "2" "3" "4" ...
head(DA5HT_Cq_AllSubjects_QCed3)
# 1      2      3      4      5      6      7      8      1      2      3      4      5      6      7      8      9
# ADRB1 24.245 23.982 24.608 24.017 23.872 24.560 24.412 23.761 23.581 24.125 24.656 24.517 24.143 24.746 24.971 24.227 24.317
# ADRB2 25.087 25.066 24.983 25.385 24.772 25.816 25.637 24.968 24.934 25.044 24.682 25.126 24.667 25.531 25.443 24.944 24.799
# COMT  24.579 24.757 24.539 24.838 24.480 24.476 25.011 24.249 23.588 24.507 25.110 24.355 24.159 24.332 24.862 24.633 24.230
# DBH   27.823 28.153 29.204 28.992 27.872 27.979 28.558 28.265 28.040 28.739 29.888 29.194 28.477 28.355 29.177 28.859 28.903
# DDC   32.773 33.446 32.634 32.706 33.048 32.914 34.089 34.110 32.681 33.205 32.143 32.642 32.937 32.687 33.408 35.332 33.250
# DRD1  25.515 24.959 25.285 25.047 24.952 25.217 25.165 25.232 25.167 24.983 25.248 25.053 24.928 25.226 25.251 25.309 25.150

#Sanity check:
DA5HT_Cq_AllSubjects_QCed3[1,]
#24.245 23.982 24.608 24.017 ...
DA5HT_MeanHousekeeping_QCed
#22.35756 22.46625 22.67619 22.37669
DA5HT_Cq_AllSubjects_QCed3[1,]-DA5HT_MeanHousekeeping_QCed
#1.8874375  1.5157500  1.9318125  1.6403125

DA5HT_NegDeltaCq_AllSubjects_QCed<-apply(DA5HT_Cq_AllSubjects_QCed3, 1, function(y) ((y-DA5HT_MeanHousekeeping_QCed)*(-1)))
str(DA5HT_NegDeltaCq_AllSubjects_QCed)
# num [1:154, 1:46] -2.16 -1.78 -2.21 -1.92 -1.86 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:154] "1" "2" "3" "4" ...
# ..$ : chr [1:46] "ADRB1" "ADRB2" "COMT" "DBH" ...

head(DA5HT_NegDeltaCq_AllSubjects_QCed)
#The apply function transposed things. Huh. I wonder why.

sum(is.na(DA5HT_Cq_AllSubjects_QCed3))
#[1] 0

pdf("DA5HT_Boxplot_Cq_perSample_AllGenes_QCed.pdf", width=20, height=5)
boxplot(data.frame(DA5HT_Cq_AllSubjects_QCed3), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq per sample for all genes", xlab="Sample ID", ylab="Cq")
dev.off()

sum(is.na(DA5HT_NegDeltaCq_AllSubjects_QCed))
#[1] 0

pdf("DA5HT_Boxplot_NegDeltaCq_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(t(DA5HT_NegDeltaCq_AllSubjects_QCed)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq per sample for all genes", xlab="Sample ID", ylab="-DeltaCq")
dev.off()
#That seems to have reduced the variability between samples. Excellent.

#***************************


#Scaling -DeltaCq and Running PCA again:

#While scaling I'm going to also change the orientation back to match my previous code:
DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed<-t(scale(DA5HT_NegDeltaCq_AllSubjects_QCed, center=T, scale=T))
head(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed)
#Sanity check:
apply(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed, 1, function(y) mean(y, na.rm=T))
apply(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed, 1, function(y) sd(y, na.rm=T))
#Looks good.

DA5HT_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed<-cor(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")

boxplot(DA5HT_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed)

pdf("DA5HT_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed.pdf", height=14, width=14)
heatmap(DA5HT_SubjectSubjectCorMatrix_NegDeltaCq_Scaled_QCed)
dev.off()
#Again - most of the replicates seem to cluster together. There are still some correlation blocks.


#Dealing with the NA values and looking at overall variation in -DeltaCq values across samples:

#This was skipped, since there aren't any NA values anymore:
#Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
# DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed[is.na(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed)]<-0
# sum(is.na(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed))
# #[1] 0

pdf("DA5HT_Boxplot_NegDeltaCqZscore_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq Z-scores per sample for all genes", xlab="Sample ID", ylab="-DeltaCq Z-score", ylim=c(-6.5,6.5))
dev.off()
#Very pretty. Definitely decreased variability.

pdf("DA5HT_Boxplot_NegDeltaCqZscore_perSample_HousekeepingGenes.pdf", width=20, height=5)
boxplot(data.frame(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed[c(31:46),]), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of -DeltaCq Z-scores per sample for housekeeping genes", xlab="Sample ID", ylab="-DeltaCq Z-score")
dev.off()
#Looks good.


#Trying PCA:

pca_DA5HT<-prcomp(t(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed))
tmp<-pca_DA5HT$x[,1:10]
dim(tmp)
#[1] 154  10
rownames(tmp)<-colnames(DA5HT_NegDeltaCq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_DA5HT_NegaDeltaCq.csv")

tmp<-pca_DA5HT$rotation[,1:10]
write.csv(tmp, "pca_DA5HT_NegDeltaCq_Eigenvectors.csv")

png("PCA_ScreePlot_DA5HT_NegDeltaCq.png")
plot(summary(pca_DA5HT)$importance[2,]~(c(1:length(summary(pca_DA5HT)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 explains a smaller portion of variance now (>12%)

png("PCA_ScreePlot2_DA5HT_NegDeltaCq.png")
plot(summary(pca_DA5HT)$importance[3,]~(c(1:length(summary(pca_DA5HT)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()

png("PC1vsPC2_DA5HT_NegDeltaCq_byDiagnosis.png")
plot(pca_DA5HT$x[,1]~pca_DA5HT$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Diagnosis))
dev.off()
#PC1 may be slightly related to diagnosis. No major outliers.

pdf("PC1_DA5HT_NegDeltaCq_byMeanHousekeeping.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~DA5HT_MeanHousekeeping_QCed)
dev.off()
#There seems to be a negative correlation with housekeeping expression - weird
cor(pca_DA5HT$x[,1],DA5HT_MeanHousekeeping_QCed)
#[1] -0.1843164

colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
cor(cbind(pca_DA5HT$x[,c(1:4)], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[,c(7,9,11,16,24:28, 38:42)]))
#                                       PC1           PC2           PC3           PC4          Age          pH        Hours.Final
# PC1                               1.000000e+00  2.689936e-16 -3.760085e-16 -2.367172e-16 -0.142236367  0.31509062  0.33205252
# PC2                               2.689936e-16  1.000000e+00 -3.742645e-16 -1.577792e-16  0.223330819  0.14273820  0.06957542
# PC3                              -3.760085e-16 -3.742645e-16  1.000000e+00  1.228779e-17 -0.211070157 -0.04820233  0.15229010
# PC4                              -2.367172e-16 -1.577792e-16  1.228779e-17  1.000000e+00  0.028039562 -0.01094156 -0.08353505
# Age                              -1.422364e-01  2.233308e-01 -2.110702e-01  2.803956e-02  1.000000000 -0.10159127  0.06708153
# pH                                3.150906e-01  1.427382e-01 -4.820233e-02 -1.094156e-02 -0.101591268  1.00000000  0.09760066
# Hours.Final                       3.320525e-01  6.957542e-02  1.522901e-01 -8.353505e-02  0.067081532  0.09760066  1.00000000
# Block.Weight..g.                  1.869037e-02  8.953521e-02 -1.668603e-02 -3.559502e-02 -0.033416539  0.05436122 -0.16210409
# TZP_BioAnalyzer_RIN                   3.152752e-01  3.584380e-01 -7.588951e-02  1.209720e-01  0.022310735  0.13296453 -0.18736430
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  3.754138e-01  3.650797e-02 -6.206581e-02 -3.203780e-01 -0.230728914  0.26879948 -0.03893164
#TZP_Average.RNAConc..ng.uL.          -6.217074e-04  1.102332e-01 -4.683776e-01 -1.528333e-01  0.020862408 -0.04190603 -0.12380585
#TZP_Average.260.280                   4.983754e-02  8.488659e-05  3.417290e-01  1.261761e-01  0.073282922 -0.10335614  0.01693539
#TZP_Average.260.230                  -1.019875e-01  3.305522e-02 -2.798308e-01  2.429189e-01 -0.005405992  0.02271813  0.01514819
#                                     Block.Weight..g. TZP_BioAnalyzer_RIN TZP_BioAnalyzer_rRNA.Ratio..28s.18s.TZP_Average.RNAConc..ng.uL.
# PC1                                    0.01869037      0.31527525                       0.37541378           -0.0006217074
# PC2                                    0.08953521      0.35843800                       0.03650797            0.1102332129
# PC3                                   -0.01668603     -0.07588951                      -0.06206581           -0.4683775970
# PC4                                   -0.03559502      0.12097196                      -0.32037800           -0.1528333317
# Age                                   -0.03341654      0.02231073                      -0.23072891            0.0208624080
# pH                                     0.05436122      0.13296453                       0.26879948           -0.0419060282
# Hours.Final                           -0.16210409     -0.18736430                      -0.03893164           -0.1238058485
# Block.Weight..g.                       1.00000000      0.19026888                       0.41739157            0.1328631054
# TZP_BioAnalyzer_RIN                        0.19026888      1.00000000                       0.50783900            0.0234041952
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.       0.41739157      0.50783900                       1.00000000            0.1468758696
#TZP_Average.RNAConc..ng.uL.                0.13286311      0.02340420                       0.14687587            1.0000000000
#TZP_Average.260.280                       -0.07283000      0.07968736                      -0.07488756           -0.6539717641
#TZP_Average.260.230                        0.03922430     -0.08018082                      -0.19885240            0.2765182975
#                                    TZP_Average.260.280TZP_Average.260.230
# PC1                                 4.983754e-02    -0.101987483
# PC2                                 8.488659e-05     0.033055223
# PC3                                 3.417290e-01    -0.279830791
# PC4                                 1.261761e-01     0.242918851
# Age                                 7.328292e-02    -0.005405992
# pH                                 -1.033561e-01     0.022718125
# Hours.Final                         1.693539e-02     0.015148188
# Block.Weight..g.                   -7.283000e-02     0.039224302
# TZP_BioAnalyzer_RIN                     7.968736e-02    -0.080180816
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   -7.488756e-02    -0.198852400
#TZP_Average.RNAConc..ng.uL.            -6.539718e-01     0.276518297
#TZP_Average.260.280                     1.000000e+00    -0.365965190
#TZP_Average.260.230                    -3.659652e-01     1.000000000

Temp<-cor(cbind(pca_DA5HT$x[,c(1:4)], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[,c(7,9,11,16,24:28, 38:42)]))

write.csv(Temp, "DA5HT_CorMatrix_PCA_vs_NumericVar_NegDeltaCq.csv")

#To be thorough, I came back and looped the basic statistical output for each variable in relationship to PC1-4:
colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)

DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats<-matrix(0, 19, 4)
row.names(DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[,c(2,5:9,11,14,16,24:28,38:42)])
colnames(DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats)<-c("PC1_ANOVA_Pval","PC2_ANOVA_Pval","PC3_ANOVA_Pval","PC4_ANOVA_Pval")

Temp<-data.frame(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[,c(2,5:9,11,14,16,24:28,38:42)], pca_DA5HT$x[,c(1:4)], ID=SubjectInfo_OrderedForDA5HTCqMatrix_QCed2[,1])

for(i in c(1:19)){
  
  Model<-lmer(Temp$PC1~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats[i,1]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp$PC2~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats[i,2]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp$PC3~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats[i,3]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
  Model<-lmer(Temp$PC4~Temp[,i]+ (1 | Temp$ID), data = Temp, REML=F)
  DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats[i,4]<-car::Anova(Model, type="III")[3][2,1]
  rm(Model)
  
}

setwd("~/Documents/Microarray Gen/FrontalPole/Output/DA5HT")
write.csv(DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats, "DA5HT_NegDeltaCq_PCA_vsSubjVariables_Stats.csv")

#PC1 correlates with pH, PMI, and RIN, 28s/18s - so things related to tissue degradation - and Card.
#PC2 correlates with RIN, cohort, and Card - looks like RNA degradation.
#PC3 correlates with 28s/18s, RNAconc, 260/280, dissection group, and Card - looks like batch effects.
#PC4 correlates with 260/230 and 28s/18s, cohort, dissection group, and 260/230, 260/280 in TTP - I'm not sure what this is. (and we're getting down to pretty low levels of variation accounted for at this level)


#Diagnosis:
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)
# (Intercept) 1.6153  1     0.2038
# Diagnosis   3.3447  2     0.1878

#Technical variables
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)                          0.0813  1    0.77557  
# Diagnosis                            2.0807  2    0.35333  
# Block.Weight..g.                     1.1447  1    0.28466  
# TZP_BioAnalyzer_RIN                  1.3430  1    0.24651  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 5.8561  1    0.01552 *
# TZP_Average.RNAConc..ng.uL.          0.0360  1    0.84956  
# TZP_Average.260.230                  0.0553  1    0.81413  
# TZP_Average.260.280                  0.0002  1    0.99019  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Basically no relationship.... except maybe with TZP_BioAnalyzer_rRNA.Ratio..28s.18s.?

#Biological variables:
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept) 10.0785  1   0.001500 **
# Diagnosis    3.2424  2   0.197665   
# Age          0.3496  1   0.554362   
# Gender       0.7160  1   0.397462   
# pH           9.8257  1   0.001721 **
# Hours.Final  6.9188  1   0.008530 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Biological and technical variables
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           3.2901  1  0.0696985 .  
# Diagnosis                             0.9679  2  0.6163420    
# Age                                   0.4337  1  0.5101690    
# Gender                                0.1544  1  0.6943726    
# pH                                    6.3073  1  0.0120240 *  
# Hours.Final                          12.8675  1  0.0003343 ***
# Block.Weight..g.                      0.3647  1  0.5459290    
# TZP_BioAnalyzer_RIN                   4.3366  1  0.0373007 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  1.6433  1  0.1998798    
# TZP_Average.RNAConc..ng.uL.           0.6306  1  0.4271271    
# TZP_Average.260.230                   0.0007  1  0.9784929    
# TZP_Average.260.280                   0.3956  1  0.5293905    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Interesting - RNA metrics less after controlling for HK expression, but other biological variables matter more.

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                  13.6604  1  0.0002190 ***
# Diagnosis                     1.6578  2  0.4365247    
# Age                           0.6583  1  0.4171760    
# Gender                        0.0169  1  0.8965375    
# pH                            6.4378  1  0.0111716 *  
# Hours.Final                  13.7206  1  0.0002121 ***
# TZP_BioAnalyzer_RIN          12.0747  1  0.0005111 ***
# TZP_Average.RNAConc..ng.uL.   0.2296  1  0.6318250    
# Card                        152.5120 19  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Biological and more technical variables & card
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)

# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            2.2232  1  0.1359504    
# Diagnosis                              1.4967  2  0.4731405    
# Age                                    0.2246  1  0.6355913    
# Gender                                 0.0000  1  0.9966902    
# pH                                     4.6300  1  0.0314183 *  
# Hours.Final                           13.1743  1  0.0002838 ***
# TZP_BioAnalyzer_RIN                    5.3773  1  0.0204005 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   2.2413  1  0.1343719    
# TZP_Average.RNAConc..ng.uL.            0.1420  1  0.7063003    
# TZP_Average.260.230                    0.0812  1  0.7757177    
# TZP_Average.260.280                    0.1920  1  0.6612507    
# Card                                 154.2672 19  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Interesting - the relationship with purity seen in the GabaGlu dataset isn't seen here.

#With dissection group instead of card:
Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)   
# (Intercept)                          1.8696  1   0.171523   
# Diagnosis                            1.5332  2   0.464599   
# Age                                  0.3982  1   0.528041   
# Gender                               0.1512  1   0.697383   
# pH                                   4.5112  1   0.033673 * 
# Hours.Final                          9.6099  1   0.001935 **
# TZP_BioAnalyzer_RIN                  2.1859  1   0.139276   
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 1.8583  1   0.172825   
# TZP_Average.RNAConc..ng.uL.          0.1216  1   0.727353   
# TZP_Average.260.230                  0.2931  1   0.588224   
# TZP_Average.260.280                  0.2194  1   0.639505   
# Dissecton.Group                      6.7201 12   0.875546   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Again, Dissection group doesn't seem to matter


#PC1 varies with pH, PMI and Card, and RIN, 28s/18s. It doesn't seem to vary with RNAconc after normalizing by HK.


png("PC1_DA5HT_NegDeltaCq_byDiagnosis.png")
boxplot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Diagnosis, ylab="PC1")
dev.off()
#pretty unconvincing
#Notably, those outlier samples are outliers within the control distribution too

pdf("PC1_DA5HT_NegDeltaCq_byDissection.pdf", width=10, height=5)
boxplot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Dissecton.Group)
dev.off()

pdf("PC1_DA5HT_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Hours.Final)
dev.off()

pdf("PC1_DA5HT_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$pH)
dev.off()

pdf("PC1_DA5HT_NegDeltaCq_byAge.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Age)
dev.off()
#Maybe Age...

pdf("PC1_DA5HT_NegDeltaCq_byGender.pdf", width=5, height=5)
boxplot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Gender)
dev.off()
#The outliers are all male, but otherwise the distributions are not different.

#I'm going to check on some other metrics that didn't matter when looking at Cq:
pdf("PC1_DA5HT_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Block.Weight..g.)
dev.off()
#no relationship

pdf("PC1_DA5HT_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC1_DA5HT_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.230)
dev.off()

pdf("PC1_DA5HT_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()
#Maybe?

pdf("PC1_DA5HT_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
dev.off()

pdf("PC1_DA5HT_NegDeltaCq_byCard.pdf", width=20, height=5)
boxplot(pca_DA5HT$x[,1]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, las=3, par(cex.axis=0.75))
dev.off()
#The effect of card seems almost perfectly imitated in the replicate card - I'm not sure this is worth including, unlike in the GabaGlu dataset.


Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Card + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)   8.3036  1   0.003957 ** 
#   Card        154.6402 19  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary.lm(lm(y~Card, data = Temp))
# Residual standard error: 2.468 on 134 degrees of freedom
# Multiple R-squared:  0.07773,	Adjusted R-squared:  -0.05304 
# F-statistic: 0.5944 on 19 and 134 DF,  p-value: 0.9052


#PC2:


#Technical variables
Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           1.1345  1  0.2868183    
# Diagnosis                             3.4415  2  0.1789326    
# Block.Weight..g.                      2.7672  1  0.0962115 .  
# TZP_BioAnalyzer_RIN                  12.4253  1  0.0004236 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  4.3169  1  0.0377350 *  
# TZP_Average.RNAConc..ng.uL.           0.9757  1  0.3232691    
# TZP_Average.260.230                   1.1365  1  0.2863983    
# TZP_Average.260.280                   0.6433  1  0.4225124    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#RIN & 28s/18s again for PC2, maybe block weight (?)

#Biological variables:
Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept) 0.3729  1    0.54145  
# Diagnosis   4.4603  2    0.10751  
# Age         3.3206  1    0.06842 .
# Gender      0.7317  1    0.39233  
# pH          0.3372  1    0.56147  
# Hours.Final 0.6572  1    0.41754  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#maybe Age?

#Biological and technical variables
Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           2.0902  1  0.1482455    
# Diagnosis                             2.9533  2  0.2284046    
# Age                                   0.7981  1  0.3716549    
# Gender                                0.5029  1  0.4782356    
# pH                                    0.5220  1  0.4699943    
# Hours.Final                           0.1288  1  0.7196300    
# Block.Weight..g.                      3.0094  1  0.0827833 .  
# TZP_BioAnalyzer_RIN                  10.9103  1  0.0009563 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  4.1531  1  0.0415584 *  
# TZP_Average.RNAConc..ng.uL.           1.2341  1  0.2666121    
# TZP_Average.260.230                   0.9095  1  0.3402499    
# TZP_Average.260.280                   1.0724  1  0.3004133    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#RIN, 28s/18s, maybe block weight.

#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   2.0662  1    0.15059    
# Diagnosis                     1.8706  2    0.39246    
# Age                           1.9100  1    0.16697    
# Gender                        0.2358  1    0.62726    
# pH                            0.1518  1    0.69686    
# Hours.Final                   0.0005  1    0.98256    
# TZP_BioAnalyzer_RIN           4.8064  1    0.02835 *  
# TZP_Average.RNAConc..ng.uL.   2.9868  1    0.08395 .  
# Card                        167.6844 19    < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#PC2 is related to RIN, pH, PMI, CARD

#Biological and more technical variables & card
Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            1.6116  1    0.20426    
# Diagnosis                              1.3090  2    0.51970    
# Age                                    1.1707  1    0.27926    
# Gender                                 0.1166  1    0.73272    
# pH                                     0.6007  1    0.43832    
# Hours.Final                            0.0220  1    0.88209    
# TZP_BioAnalyzer_RIN                    5.8227  1    0.01582 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   1.0899  1    0.29650    
# TZP_Average.RNAConc..ng.uL.            4.0367  1    0.04452 *  
# TZP_Average.260.230                    0.2074  1    0.64884    
# TZP_Average.260.280                    0.6333  1    0.42615    
# Card                                 166.6208 19    < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Maybe RNA conc in addition to RIN?  Card seems important.

#With dissection group instead of card:
Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           1.8278  1    0.17638    
# Diagnosis                             2.4115  2    0.29946    
# Age                                   0.8231  1    0.36429    
# Gender                                2.0504  1    0.15217    
# pH                                    0.1167  1    0.73264    
# Hours.Final                           0.0517  1    0.82008    
# TZP_BioAnalyzer_RIN                  19.0411  1  1.279e-05 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.7391  1    0.38996    
# TZP_Average.RNAConc..ng.uL.           1.4094  1    0.23516    
# TZP_Average.260.230                   0.1034  1    0.74778    
# TZP_Average.260.280                   0.5920  1    0.44165    
# Dissecton.Group                      25.9205 12    0.01102 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Mabye dissection group, but it is related to Card so probably not. 


png("PC2_DA5HT_NegDeltaCq_byDiagnosis.png")
boxplot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Diagnosis, ylab="PC2")
dev.off()
#Nope

pdf("PC2_DA5HT_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$pH)
dev.off()

pdf("PC2_DA5HT_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Hours.Final)
dev.off()

pdf("PC2_DA5HT_NegDeltaCq_byRIN.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_RIN)
dev.off()
#Very convincing.

#Double-checking some of the other metrics that didn't matter for Cq:

pdf("PC2_DA5HT_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Block.Weight..g.)
dev.off()
#Maybe?

pdf("PC2_DA5HT_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.280)
dev.off()
#no relationship

pdf("PC2_DA5HT_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.230)
dev.off()

pdf("PC2_DA5HT_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()
#U-shaped

pdf("PC2_DA5HT_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
dev.off()
#U-shaped

pdf("PC2_DA5HT_NegDeltaCq_byCard.pdf", width=20, height=5)
boxplot(pca_DA5HT$x[,2]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, las=3, par(cex.axis=0.75))
dev.off()
#The effect of Card mostly seems parallel in replicate cards again.

#PC3:

#Technical variables
Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)                          0.1959  1    0.65805  
# Diagnosis                            1.7004  2    0.42733  
# Block.Weight..g.                     0.5142  1    0.47331  
# TZP_BioAnalyzer_RIN                  0.1047  1    0.74627  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 2.1457  1    0.14297  
# TZP_Average.RNAConc..ng.uL.          5.4952  1    0.01907 *
# TZP_Average.260.230                  0.1075  1    0.74298  
# TZP_Average.260.280                  0.4657  1    0.49498  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Maybe RNAConc

#Biological variables:
Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept) 0.5591  1    0.45464  
# Diagnosis   4.1202  2    0.12744  
# Age         2.7900  1    0.09485 .
# Gender      4.8688  1    0.02735 *
# pH          0.6753  1    0.41120  
# Hours.Final 4.8425  1    0.02777 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Gender...really?  and PMI, maybe Age.


#Biological and technical variables
Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Block.Weight..g.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)  
# (Intercept)                          0.0272  1    0.86901  
# Diagnosis                            2.1567  2    0.34016  
# Age                                  5.8264  1    0.01579 *
# Gender                               1.4247  1    0.23263  
# pH                                   0.1301  1    0.71828  
# Hours.Final                          3.3375  1    0.06772 .
# Block.Weight..g.                     0.4287  1    0.51264  
# TZP_BioAnalyzer_RIN                  0.5819  1    0.44558  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. 3.3093  1    0.06889 .
# TZP_Average.RNAConc..ng.uL.          4.7324  1    0.02960 *
# TZP_Average.260.230                  0.5769  1    0.44752  
# TZP_Average.260.280                  0.3507  1    0.55369  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Age, RNAconc, maybe PMI, 28s/18s
#I'm guessing the effect of "gender" was actually standing in for RNAconc (which differes a little by gender, for whatever reason)


#Biological and the most important technical variables & Card:
Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                   1.2766  1    0.25853    
# Diagnosis                     2.3517  2    0.30855    
# Age                           4.0865  1    0.04323 *  
# Gender                        3.1591  1    0.07551 .  
# pH                            0.6846  1    0.40802    
# Hours.Final                   2.7124  1    0.09957 .  
# TZP_BioAnalyzer_RIN           0.0644  1    0.79963    
# TZP_Average.RNAConc..ng.uL.   4.9768  1    0.02569 *  
# Card                        112.7376 19  2.489e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Age, RNAConc, Card, maybe PMI, maybe Gender.

#Biological and more technical variables & card
Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                            0.2014  1    0.65357    
# Diagnosis                              2.3513  2    0.30862    
# Age                                    4.5901  1    0.03216 *  
# Gender                                 3.0845  1    0.07904 .  
# pH                                     0.4397  1    0.50726    
# Hours.Final                            2.8931  1    0.08896 .  
# TZP_BioAnalyzer_RIN                    0.0068  1    0.93421    
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   0.4267  1    0.51360    
# TZP_Average.RNAConc..ng.uL.            3.1662  1    0.07518 .  
# TZP_Average.260.230                    0.3648  1    0.54583    
# TZP_Average.260.280                    0.0053  1    0.94180    
# Card                                 110.0282 19  7.877e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Age, Card, Maybe Gender, PMI, RNAConc

#With dissection group instead of card:
Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#Dissection group maybe, but probably mediated by Card


png("PC3_DA5HT_NegDeltaCq_byDiagnosis.png")
boxplot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Diagnosis, ylab="PC2")
dev.off()
#Nope

pdf("PC3_DA5HT_NegDeltaCq_byAge.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Age)
dev.off()
#Maybe?

pdf("PC3_DA5HT_NegDeltaCq_byGender.pdf", width=5, height=5)
boxplot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Gender)
dev.off()
#Maybe? the sample size for F is pretty small...

pdf("PC3_DA5HT_NegDeltaCq_byPH.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$pH)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_byPMI.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Hours.Final)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_byRIN.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_RIN)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_byRNAConc.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
dev.off()
#Reasonably convincing

pdf("PC3_DA5HT_NegDeltaCq_byBlockWeight.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Block.Weight..g.)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_by260280.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.280)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_by260230.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_Average.260.230)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_by28s18s.pdf", width=5, height=5)
plot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
dev.off()

pdf("PC3_DA5HT_NegDeltaCq_byCard.pdf", width=20, height=5)
boxplot(pca_DA5HT$x[,3]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, las=3, par(cex.axis=0.75))
dev.off()
#Replicate cards seem pretty similar again.

#This is code that was included in parallel with GabaGlu, but maybe doesn't make a lot of sense except to demonstrate that this dataset doesn't show the same pattern of Card-level effects as the GabaGlu dataset:

#Is there some way to collapse card? (so that it doesn't take up so many df in the model???)

pca_DA5HT_byCard<-matrix(0, 20, 10)
row.names(pca_DA5HT_byCard)<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card))

for(i in c(1:10)){
  pca_DA5HT_byCard[,i]<-tapply(X=pca_DA5HT$x[,i], INDEX=SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card, FUN=mean)
}


pca_DA5HT_byCard_CorMatrix<-cor(t(pca_DA5HT_byCard))
dim(pca_DA5HT_byCard_CorMatrix)
row.names(pca_DA5HT_byCard_CorMatrix)<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card))
colnames(pca_DA5HT_byCard_CorMatrix)<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed2$Card))

pdf("Heatmap_pca_DA5HT_byCard_CorMatrix.pdf", width=10, height=10)
heatmap(pca_DA5HT_byCard_CorMatrix)
dev.off()
#Unlike the DA5HT experiment, there are no large obvious blocks - the blocks mostly seem to be defined by replicate cards, although a few replicates don't match up quite as well as they probably should (e.g., 43 & 44, 49 & 50).

pca_DA5HT_byCard_CorMatrix
#The correlation between replicate cards is typically >0.50

#Example plot:
row.names(pca_DA5HT_byCard)
plot(pca_DA5HT_byCard[1,]~pca_DA5HT_byCard[2,])
#the correlation is almost entirely driven by what I'm assuming are PC1 and PC2 or PC2.

write.csv(pca_DA5HT_byCard_CorMatrix, "pca_DA5HT_byCard_NegDeltaCq_CorMatrix.csv")


#*************************


#There isn't any final outlier removal for the DA5HT dataset

DA5HT_NegDeltaCq_AllSubjects_QCed2<-DA5HT_NegDeltaCq_AllSubjects_QCed
str(DA5HT_NegDeltaCq_AllSubjects_QCed2)
# num [1:154, 1:46] -2.16 -1.78 -2.21 -1.92 -1.86 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:154] "1" "2" "3" "4" ...
# ..$ : chr [1:46] "ADRB1" "ADRB2" "COMT" "DBH" ...

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed2
str(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
# 'data.frame':	154 obs. of  47 variables:
#   $ ID                                  : chr  "1" "2" "3" "4" ...
# $ Card                                : chr  "Card 39.eds" "Card 39.eds" "Card 39.eds" "Card 39.eds" ...
# $ Subject.Number                      : int  5000 5313 2169 4087 5066 2950 4383 4619 5000 5313 ...
# $ Barcode                             : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Cohort                              : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis                           : Factor w/ 3 levels "Control","BP",..: 1 3 1 2 1 3 2 3 1 3 ...
# $ Age                                 : int  77 42 18 51 47 45 40 34 77 42 ...
# $ Gender                              : Factor w/ 2 levels "M","F": 1 1 1 2 1 1 1 1 1 1 ...
# $ pH                                  : num  6.79 6.77 6.97 6.81 6.6 7.05 6.77 6.7 6.79 6.77 ...
# $ AFS                                 : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final                         : num  16.5 23.7 22 16.8 21 20 26.3 17.2 16.5 23.7 ...
# $ Slab.Format                         : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number                         : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group                     : Factor w/ 13 levels "1","2","3","4",..: 1 1 1 1 1 1 2 2 1 1 ...
# $ TTP.ID                              : chr  "TTP010001" "TTP010002" "TTP010003" "TTP010004" ...
# $ Block.Weight..g.                    : num  0.72 0.66 0.76 0.56 0.7 0.89 0.85 1.02 0.72 0.66 ...
# $ TZP.ID                              : chr  "TZP010001" "TZP010002" "TZP010003" "TZP010004" ...
# $ TZP_MBNI_Nanodrop_RNAConc..ng.uL.   : num  436 343 394 417 399 ...
# $ TZP_MBNI_Nanodrop_260.280           : num  1.96 2.01 1.99 1.97 2 1.99 2 2 1.96 2.01 ...
# $ TZP_MBNI_Nanodrop_260.230           : num  1.93 2.18 2.02 2.1 2.19 2.23 2.3 2.3 1.93 2.18 ...
# $ TZP_DSC_Nanodrop_RNAConc..ng.uL.    : num  562 354 451 449 435 ...
# $ TZP_DSC_Nanodrop_260.280            : num  2.14 1.92 2.03 2.03 2.05 2.05 2.05 2.06 2.14 1.92 ...
# $ TZP_DSC_Nanodrop_260.230            : num  1.9 2.06 2.01 2.13 2.22 2.27 2.29 2.28 1.9 2.06 ...
# $ TZP_BioAnalyzer_RIN                 : num  7.7 8.1 7.9 8 8.2 8.1 6.6 7.7 7.7 8.1 ...
# $ TZP_BioAnalyzer_rRNA.Ratio..28s.18s.: num  2.4 2.3 3.1 2 2.4 2.7 1.44 2.4 2.4 2.3 ...
# $ TZP_Average.RNAConc..ng.uL.         : num  499 349 422 433 417 ...
# $ TZP_Average.260.280                 : num  2.05 1.97 2.01 2 2.03 2.02 2.03 2.03 2.05 1.97 ...
# $ TZP_Average.260.230                 : num  1.92 2.12 2.02 2.12 2.21 2.25 2.3 2.29 1.92 2.12 ...
# $ TTP.ID.1                            : chr  "TTP010001" "TTP010002" "TTP010003" "TTP010004" ...
# $ Initial.Volume..uL.                 : int  720 660 760 560 700 890 850 1020 720 660 ...
# $ Final.Volume..uL.                   : int  620 560 660 460 600 790 750 920 620 560 ...
# $ TTP_MBNI_Nanodrop_RNAConc..ng.uL.   : num  528 398 452 523 463 ...
# $ TTP_MBNI_Nanodrop_260.280           : num  2 1.82 1.77 1.96 1.78 1.8 1.78 1.79 2 1.82 ...
# $ TTP_MBNI_Nanodrop_260.230           : num  1.87 2.02 1.99 1.99 1.9 2.04 2.24 1.94 1.87 2.02 ...
# $ TTP_DSC_Nanodrop_RNAConc..ng.uL.    : num  631 429 555 644 551 ...
# $ TTP_DSC_Nanodrop_260.280            : num  2.02 1.92 2.01 2.02 2.03 1.91 1.97 1.98 2.02 1.92 ...
# $ TTP_DSC_Nanodrop_260.230            : num  1.83 2.04 2.05 1.96 1.9 2.09 2.29 1.95 1.83 2.04 ...
# $ TTP_BioAnalyzer_RIN                 : num  7.4 7.8 7.7 7.7 8 7.9 7 7.2 7.4 7.8 ...
# $ TTP_BioAnalyzer_rRNA.Ratio..28s.18s.: num  2.2 2.2 2.8 2.1 2.4 2.5 1.9 2.3 2.2 2.2 ...
# $ TTP_Average.RNAConc..ng.uL.         : num  579 414 503 583 507 ...
# $ TTP_Average.260.280                 : num  2.01 1.87 1.89 1.99 1.91 ...
# $ TTP_Average.260.230                 : num  1.85 2.03 2.02 1.98 1.9 ...
# $ MissingPH                           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ MissingExtractionData               : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ LowQualityRNA                       : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ SampleNumber                        : logi  NA NA NA NA NA NA ...
# $ QC_AmplificationIssue               : chr  "N" "N" "N" "N" ...

DA5HT_MeanHousekeeping_QCed2<-DA5HT_MeanHousekeeping_QCed
length(DA5HT_MeanHousekeeping_QCed2)
#[1] 154

#After final "outlier removal", double check the relationship between variables (e.g., diagnosis and gender)

length(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$ID)
#[1] 154
table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$ID)
# 1 10 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 37 38 39  4 40 41 42 44 45 46 47 48 49  5 51 52 53 54 55 56 
# 4  2  2  2  2  2  2  4  2  2  2  2  4  2  2  2  2  2  2  2  2  2  4  2  2  2  2  2  2  2  2  4  2  2  2  4  2  2  2  2  2  2  2  2  2  2 
# 57 58 59  6 60 61 62 63 64 65 66 67 68 69  7 70 71 72 73 74 75  8  9 
# 2  4  2  2  2  4  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 

#No subjects have been removed beyond the original basic QC, so the relationship between diagnosis and subject variables should be the same as identified earlier.


table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Card)
# Card 39.eds Card 40.eds Card 41.eds Card 42.eds Card 43.eds Card 44.eds Card 45.eds Card 46.eds Card 47.eds Card 48.eds
# Control           3           3           3           3           3           3           2           2           4           4
# BP                2           2           3           3           2           2           2           2           1           1
# Schiz             3           3           1           1           3           3           4           4           2           2
# 
# Card 49.eds Card 50.eds Card 51.eds Card 52.eds Card 53.eds Card 54.eds Card 55.eds Card 56.eds Card 57.eds Card 58.eds
# Control           3           3           2           2           4           4           2           2           4           4
# BP                3           3           3           3           2           2           3           3           2           2
# Schiz             1           1           3           3           2           2           3           3           2           2

write.csv(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Card), "Table_DA5HT_DiagnosisVsCard.csv")

library(nlme)

#The relationship with housekeeping gene expression should be the same as evaluated earlier also (no outlier removal)



pdf("Boxplot_DA5HT_RIN_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_BioAnalyzer_RIN~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_BioAnalyzer_RIN~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

Model<-lme(TZP_BioAnalyzer_RIN~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, na.action = na.omit, method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForDA5HTCqMatrix_QCed3 
# AIC      BIC   logLik
# -5462.954 -5447.77 2736.477
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)    Residual
# StdDev:   0.5629399 9.29243e-16
# 
# Fixed effects: TZP_BioAnalyzer_RIN ~ Diagnosis 
# Value Std.Error DF  t-value p-value
# (Intercept)     7.817016 0.1106491 85 70.64693  0.0000
# DiagnosisBP     0.036049 0.1654861 66  0.21784  0.8282
# DiagnosisSchiz -0.177916 0.1633899 66 -1.08891  0.2802
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.669       
# DiagnosisSchiz -0.677  0.453
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.9116172 -0.9558086  0.0000000  0.9558086  2.8674257 
# 
# Number of Observations: 154
# Number of Groups: 69 


pdf("Boxplot_DA5HT_RNAconc_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

Model<-lme(TZP_Average.RNAConc..ng.uL.~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, na.action = na.omit, method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForDA5HTCqMatrix_QCed3 
# AIC       BIC   logLik
# -4117.041 -4101.856 2063.521
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:    63.63683 5.492186e-14
# 
# Fixed effects: TZP_Average.RNAConc..ng.uL. ~ Diagnosis 
# Value Std.Error DF  t-value p-value
# (Intercept)    394.0458  12.26479 85 32.12822  0.0000
# DiagnosisBP      4.6715  18.32986 66  0.25486  0.7996
# DiagnosisSchiz -14.9536  18.09700 66 -0.82630  0.4116
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.669       
# DiagnosisSchiz -0.678  0.453
# 
# Standardized Within-Group Residuals:
#   Min        Q1       Med        Q3       Max 
# -3.104961 -1.034987  0.000000  1.034987  3.104961 
# 
# Number of Observations: 154
# Number of Groups: 69


# Model<-lme(TZP_Average.260.280~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
# summary(Model)
# Error in solve.default(pdMatrix(a, factor = TRUE)) : 
#   Lapack routine dgesv: system is exactly singular: U[1,1] = 0
#I at first thought this was being caused by ID being accidentally coded as an integer, but I went back and fixed that, re-ran everything, and discovered it didn't change any of the results (?!)
#Ah -ha - the problem is the lmeControl(opt='optim') argument. Huh.

Model<-lme(TZP_Average.260.280~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, na.action = na.omit,  method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForDA5HTCqMatrix_QCed3 
# AIC       BIC   logLik
# -6078.809 -6063.624 3044.404
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:  0.02121032 3.553915e-16
# 
# Fixed effects: TZP_Average.260.280 ~ Diagnosis 
# Value   Std.Error DF  t-value p-value
# (Intercept)     2.0293408 0.004200492 85 483.1198  0.0000
# DiagnosisBP    -0.0083136 0.006284029 66  -1.3230  0.1904
# DiagnosisSchiz -0.0075074 0.006204519 66  -1.2100  0.2306
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.668       
# DiagnosisSchiz -0.677  0.453
# 
# Standardized Within-Group Residuals:
#   Min        Q1       Med        Q3       Max 
# -1.249578 -1.249578  0.000000  0.000000  1.249578 
# 
# Number of Observations: 154
# Number of Groups: 69 


pdf("Boxplot_DA5HT_260280_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_Average.260.280~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Quality (260/280)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_Average.260.280~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


Model<-lme(TZP_Average.260.230~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, na.action = na.omit,  method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForDA5HTCqMatrix_QCed3 
# AIC       BIC   logLik
# -5815.886 -5800.701 2912.943
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:   0.1754663 3.002446e-16
# 
# Fixed effects: TZP_Average.260.230 ~ Diagnosis 
# Value  Std.Error DF  t-value p-value
# (Intercept)     2.0610358 0.03450696 85 59.72812  0.0000
# DiagnosisBP    -0.0066054 0.05160944 66 -0.12799  0.8985
# DiagnosisSchiz  0.0236215 0.05095574 66  0.46357  0.6445
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.669       
# DiagnosisSchiz -0.677  0.453
# 
# Standardized Within-Group Residuals:
#   Min        Q1       Med        Q3       Max 
# -1.479091  0.000000  0.000000  0.000000  1.479091 
# 
# Number of Observations: 154
# Number of Groups: 69 

pdf("Boxplot_DA5HT_260230_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_Average.260.230~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Quality (260/230)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_Average.260.230~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

Model<-lme(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, na.action = na.omit, method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForDA5HTCqMatrix_QCed3 
# AIC   BIC   logLik
# -5738.185 -5723 2874.092
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:   0.4024914 2.417044e-16
# 
# Fixed effects: TZP_BioAnalyzer_rRNA.Ratio..28s.18s. ~ Diagnosis 
# Value  Std.Error DF  t-value p-value
# (Intercept)     2.3994084 0.07546945 85 31.79311  0.0000
# DiagnosisBP     0.0555091 0.11268140 66  0.49262  0.6239
# DiagnosisSchiz -0.0335971 0.11124447 66 -0.30201  0.7636
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.670       
# DiagnosisSchiz -0.678  0.454
# 
# Standardized Within-Group Residuals:
#   Min        Q1       Med        Q3       Max 
# -1.837323  0.000000  1.837323  1.837323  3.674647 
# 
# Number of Observations: 154
# Number of Groups: 69 

pdf("Boxplot_DA5HT_rRNA28s18s_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="TZP RNA Quality (28s/18s rRNA ratio)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


#I inverted this to make the results more intuitive:
pdf("Boxplot_DA5HT_meanHK_ByDiagnosis.pdf", width=4, height=6)
boxplot((-1*DA5HT_MeanHousekeeping_QCed2)~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart((-1*DA5HT_MeanHousekeeping_QCed2)~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


pdf("Boxplot_DA5HT_MeanHKexpression_vs_Card.pdf", width=10, height=5)
boxplot(y~Card, data=Temp, ylab="Mean Housekeeping Gene Expression", las=2)
dev.off()



#Trying stepwise to see if that helps increase my confidence regarding which variables to keep/toss from the model:


library(lmerTest)
#This is backward elimination for mixed-effects models.

Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + (1 | ID)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + Card + (1 | ID)
#hmm... interesting. A little different, but not much - contains 28s/18s when card is in the model. 

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + Card + (1 | ID)


Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + TZP_BioAnalyzer_RIN + Dissecton.Group + (1 | ID)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + TZP_BioAnalyzer_RIN + Card + (1 | ID)
#Interesting.

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + TZP_BioAnalyzer_RIN + Dissecton.Group + Card + (1 | ID)


Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Gender + TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)


Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Age + Gender + Hours.Final + TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.RNAConc..ng.uL.+TZP_Average.260.230 +TZP_Average.260.280+Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)

step(Model3, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.1, 
     reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, keep = c("Diagnosis", "Gender"))

# Model found:
#   y ~ Diagnosis + Age + Gender + pH + Hours.Final + Dissecton.Group + Card + (1 | ID)
# Interesting - in the DA5HT dataset we can hypothetically include both Dissection Group and Card (!!!), although it is a huge amount of df again.

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Dissecton.Group, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Card)
#Variability due to dissection group is mostly captured by card - each card typically only contains samples from 2 dissection groups.


#*****************************

#I'm still nervous about all of the df used here. Feels like overfitting.
#Again, evaluating whether card, dissection group, and 28s/18s is worth having in the model:

Temp<-data.frame(y=pca_DA5HT$x[,1], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Card+ (1 | ID), data = Temp, REML=F)

Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Dissecton.Group+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 514.26 547.66 -246.13   492.26                             
# Model1 30 466.32 557.43 -203.16   406.32 85.935     19  1.727e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#AIC is better for the model with Card, BIC is better for the model without Card.
#According to the ANOVA, adding Card is a significant improvement in model fit.

#What about 28s/18s?

anova(Model, Model2)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model2: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model2:     TZP_Average.RNAConc..ng.uL. + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + 
#   Model2:     (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 514.26 547.66 -246.13   492.26                         
# Model2 12 514.79 551.23 -245.39   490.79 1.4698      1     0.2254

#Not an improvement.

anova(Model, Model3)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 514.26 547.66 -246.13   492.26                         
# Model3 23 532.92 602.77 -243.46   486.92 5.3413     12     0.9456


anova(Model, Model4)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + Card + (1 | 
#                                                                         Model4:     ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 514.26 547.66 -246.13   492.26                             
# Model4 42 485.51 613.06 -200.75   401.51 90.751     31  8.994e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Better AIC but worse BIC, sig ANOVA

anova(Model1, Model4)
# Data: Temp
# Models:
#   Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + Card + (1 | 
#                                                                         Model4:     ID)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# Model1 30 466.32 557.43 -203.16   406.32                        
# Model4 42 485.51 613.06 -200.75   401.51 4.816     12     0.9638

#Dissection group is not an improvement - fits stepwise results.


#Let's take a peek at PC2 also: 

Temp<-data.frame(y=pca_DA5HT$x[,2], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)

Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Card+ (1 | ID), data = Temp, REML=F)

Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Dissecton.Group+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)

# Data: Temp
# Models:
# Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
# Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
# Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 510.14 543.55 -244.07   488.14                             
# Model1 30 453.10 544.21 -196.55   393.10 95.044     19  4.185e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#AIC is better for the model with Card, BIC is better for the model without Card.
#According to the ANOVA, adding Card is a significant improvement in model fit.

anova(Model, Model2)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model2: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model2:     TZP_Average.RNAConc..ng.uL. + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + 
#   Model2:     (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 510.14 543.55 -244.07   488.14                         
# Model2 12 510.80 547.25 -243.40   486.80 1.3433      1     0.2465

#No improvement with 28s/18s

anova(Model, Model3)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Model  11 510.14 543.55 -244.07   488.14                           
# Model3 23 510.80 580.65 -232.40   464.80 23.349     12     0.0249 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Slight improvement with dissection group, but could actually be card, and AIC and BIC are both actually worse with dissection.

anova(Model, Model4)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + Card + (1 | 
#                                                                         Model4:     ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 510.14 543.55 -244.07   488.14                             
# Model4 42 450.51 578.07 -183.26   366.51 121.63     31  1.115e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Better AIC, worse BIC, sig ANOVA.

anova(Model1, Model4)

# Data: Temp
# Models:
#   Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + Card + (1 | 
#                                                                         Model4:     ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Model1 30 453.10 544.21 -196.55   393.10                            
# Model4 42 450.51 578.07 -183.26   366.51 26.586     12   0.008859 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Better AIC, worse BIC, sig ANOVA compared to model with just card. Huh. So dissection group is an improvement maybe, even after controlling for card.


#One more check: PC3

Temp<-data.frame(y=pca_DA5HT$x[,3], SubjectInfo_OrderedForDA5HTCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + (1 | ID), data = Temp, REML=F)

Model1<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Card+ (1 | ID), data = Temp, REML=F)

Model2<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+ (1 | ID), data = Temp, REML=F)

Model3<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Dissecton.Group+ (1 | ID), data = Temp, REML=F)

Model4<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL. + Dissecton.Group+Card+ (1 | ID), data = Temp, REML=F)

anova(Model, Model1)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Model  11 499.08 532.48 -238.54   477.08                             
# Model1 30 466.69 557.80 -203.35   406.69 70.382     19  7.946e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(Model, Model2)
# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model2: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model2:     TZP_Average.RNAConc..ng.uL. + TZP_BioAnalyzer_rRNA.Ratio..28s.18s. + 
#   Model2:     (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Model  11 499.08 532.48 -238.54   477.08                         
# Model2 12 498.65 535.09 -237.32   474.65 2.4268      1     0.1193

anova(Model, Model3)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model3: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model3:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + (1 | ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Model  11 499.08 532.48 -238.54   477.08                           
# Model3 23 502.86 572.71 -228.43   456.86 20.216     12    0.06312 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(Model, Model4)

# Data: Temp
# Models:
#   Model: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model:     TZP_Average.RNAConc..ng.uL. + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + Card + (1 | 
#                                                                         Model4:     ID)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# Model  11 499.08 532.48 -238.54   477.08                            
# Model4 42 468.58 596.13 -192.29   384.58  92.5     31  4.908e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Better AIC, worse BIC, sig ANOVA

anova(Model1, Model4)

# Data: Temp
# Models:
#   Model1: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model1:     TZP_Average.RNAConc..ng.uL. + Card + (1 | ID)
# Model4: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN + 
#   Model4:     TZP_Average.RNAConc..ng.uL. + Dissecton.Group + Card + (1 | 
#                                                                         Model4:     ID)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Model1 30 466.69 557.80 -203.35   406.69                           
# Model4 42 468.58 596.13 -192.29   384.58 22.118     12    0.03622 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Worse AIC, worse BIC, slightly sig ANOVA when adding dissection group on top of Card. Hmmm...


#Card may or may not be a good co-variate (depending on metric), dissection group is a good addition some of the time (PC2).




#*********************************************
