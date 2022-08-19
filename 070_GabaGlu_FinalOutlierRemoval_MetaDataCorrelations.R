#**************************

#Final Outlier removal

str(GabaGlu_NegDeltaCq_AllSubjects_QCed)
# num [1:136, 1:93] 0.934 1.447 1.575 1.674 1.626 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:136] "1" "2" "3" "4" ...
# ..$ : chr [1:93] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...

GabaGlu_NegDeltaCq_AllSubjects_QCed2<-GabaGlu_NegDeltaCq_AllSubjects_QCed[pca_GabaGlu$x[,1]>(-15),]
str(GabaGlu_NegDeltaCq_AllSubjects_QCed2)
# num [1:133, 1:93] 1.45 1.57 1.67 1.63 1.68 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:133] "2" "3" "4" "17" ...
# ..$ : chr [1:93] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...

sum(is.na(GabaGlu_NegDeltaCq_AllSubjects_QCed2))
#[1] 16
length(is.na(GabaGlu_NegDeltaCq_AllSubjects_QCed2))
#[1] 12369
16/12369
#0.001293556
#0.13% of measurements are NA following sample-level and probe-level QC

str(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
# 'data.frame':	136 obs. of  50 variables:

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[pca_GabaGlu$x[,1]>(-15),]
str(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
# 'data.frame':	133 obs. of  50 variables:

# length(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$ID))
# [1] 69

str(GabaGlu_Cq_AllSubjects_QCed3)
# num [1:93, 1:136] 20.5 25.7 23.7 27.3 22.8 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:93] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# ..$ : chr [1:136] "1" "2" "3" "4" ...

GabaGlu_Cq_AllSubjects_QCed4<-GabaGlu_Cq_AllSubjects_QCed3[,pca_GabaGlu$x[,1]>(-15)]
str(GabaGlu_Cq_AllSubjects_QCed4)
# num [1:93, 1:133] 20.7 25 23.1 26.8 22.8 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:93] "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# ..$ : chr [1:133] "2" "3" "4" "17" ...

GabaGlu_MeanHousekeeping_QCed2<-GabaGlu_MeanHousekeeping_QCed[pca_GabaGlu$x[,1]>(-15)]
length(GabaGlu_MeanHousekeeping_QCed2)
#[1] 133

dim(GabaGlu_Cq_AllSubjects_QCed3_Imputed[,pca_GabaGlu$x[,1]>(-15)])
GabaGlu_Cq_AllSubjects_QCed3_Imputed<-GabaGlu_Cq_AllSubjects_QCed3_Imputed[,pca_GabaGlu$x[,1]>(-15)]
dim(GabaGlu_Cq_AllSubjects_QCed3_Imputed)
#[1]  93 133



###########################

#After final outlier removal, double check the relationship between variables (e.g., diagnosis and gender)

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$ID)
# 1 10 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 37 38 39  4 40 41 42 44 45 46 47 48 49  5 51 52 53 54 55 56 
# 1  2  2  2  2  1  2  2  2  2  2  2  2  2  2  2  2  1  2  2  2  2  2  2  2  2  1  1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 
# 57 58 59  6 60 61 62 63 64 65 66 67 68 69  7 70 71 72 73 74 75  8  9 
# 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 

length(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$ID))
#[1]  69

#We haven't removed any additional subjects beyond the first Basic QC culling, so the relationships between subject level variables and diagnosis should be the same as reported previously.

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$CardBlock2)
#           EvenCards OddCards
# Control        26       23
# BP             20       21
# Schiz          21       22
#Pretty balanced - as it should be.

write.csv(table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card), "Table_GabaGlu_DiagnosisVsCard.csv")


Model<-lme(TZP_BioAnalyzer_RIN~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForGabaGluCqMatrix_QCed3 
# AIC       BIC  logLik
# -4070.3 -4055.848 2040.15
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:    0.562509 9.865146e-16
# 
# Fixed effects: TZP_BioAnalyzer_RIN ~ Diagnosis 
# Value Std.Error DF  t-value p-value
# (Intercept)     7.822136 0.1107999 66 70.59697  0.0000
# DiagnosisBP     0.042207 0.1657105 66  0.25470  0.7997
# DiagnosisSchiz -0.181792 0.1636133 66 -1.11111  0.2706
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.669       
# DiagnosisSchiz -0.677  0.453
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.8006392 -0.9003196  0.0000000  0.9003196  2.7009588 
# 
# Number of Observations: 133
# Number of Groups: 69 

pdf("Boxplot_GabaGlu_RIN_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_RIN~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Integrity Number (RIN)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_RIN~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


Model<-lme(TZP_Average.RNAConc..ng.uL.~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForGabaGluCqMatrix_QCed3 
# AIC       BIC  logLik
# -2893.041 -2878.589 1451.52
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:    63.44478 5.968386e-14
# 
# Fixed effects: TZP_Average.RNAConc..ng.uL. ~ Diagnosis 
# Value Std.Error DF  t-value p-value
# (Intercept)    394.6535  12.28618 66 32.12175  0.0000
# DiagnosisBP      6.2822  18.36221 66  0.34212  0.7333
# DiagnosisSchiz -15.3961  18.12969 66 -0.84922  0.3988
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.669       
# DiagnosisSchiz -0.678  0.453
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.8572258  0.0000000  0.0000000  0.9524086  2.8572258 
# 
# Number of Observations: 133
# Number of Groups: 69 

pdf("Boxplot_GabaGlu_RNAconc_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Concentration (ng/uL)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.RNAConc..ng.uL.~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

Model<-lme(TZP_Average.RNAConc..ng.uL.~Block.Weight..g.+Diagnosis+Gender+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
#                                         Value Std.Error DF   t-value p-value
# (Intercept)                            5.76157  48.85719 64  0.117927  0.9065
# Block.Weight..g.                     147.62508  38.41338 59  3.843064  0.0003
# DiagnosisBP                            9.50395  19.06398 59  0.498529  0.6200
# DiagnosisSchiz                       -19.44941  18.18246 59 -1.069679  0.2891
# GenderF                              -66.98875  26.00615 59 -2.575881  0.0125
# Age                                    0.54556   0.52048 59  1.048201  0.2988
# pH                                     7.62109   6.53038 59  1.167022  0.2479
# Hours.Final                            0.28630   0.95469 59  0.299883  0.7653
# TZP_BioAnalyzer_RIN                   20.03800   5.76405 59  3.476377  0.0010
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  21.14877  17.01854 59  1.242690  0.2189

#RNA concentration tracks block weight and RIN, maybe gender.  
#Block weight could represent something biological in addition to dissection variation.

Model<-lme((TZP_Average.RNAConc..ng.uL./Block.Weight..g.)~Diagnosis+Gender+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
#After controlling for block weight, RNAconc is mildly related to RIN, PMI, pH, and gender.

# Value Std.Error DF    t-value p-value
# (Intercept)                            10.69910  84.14234 64  0.1271547  0.8992
# DiagnosisBP                            14.08143  43.56789 60  0.3232065  0.7477
# DiagnosisSchiz                        -40.43979  41.09728 60 -0.9840017  0.3291
# GenderF                              -125.96799  56.06460 60 -2.2468364  0.0283
# Age                                     1.17677   1.09488 60  1.0747874  0.2868
# pH                                     22.51251  10.83858 60  2.0770718  0.0421
# Hours.Final                             3.47886   2.05437 60  1.6933962  0.0956
# TZP_BioAnalyzer_RIN                    22.72283   9.87362 60  2.3013667  0.0249
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   40.36859  30.75755 60  1.3124777  0.1944

Model<-lme(Block.Weight..g.~Diagnosis+Gender+Age+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
#                                           Value  Std.Error DF   t-value p-value
# (Intercept)                          -0.02682017 0.11208633 64 -0.239281  0.8117
# DiagnosisBP                           0.01116423 0.04534073 60  0.246230  0.8063
# DiagnosisSchiz                        0.07570104 0.04270969 60  1.772456  0.0814
# GenderF                               0.20450406 0.05889729 60  3.472215  0.0010
# Age                                   0.00225226 0.00121290 60  1.856927  0.0682
# pH                                    0.03014427 0.01471569 60  2.048444  0.0449
# Hours.Final                           0.00040759 0.00225697 60  0.180592  0.8573
# TZP_BioAnalyzer_RIN                  -0.00178814 0.01323400 60 -0.135117  0.8930
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  0.15464296 0.03713945 60  4.163846  0.0001

#Block weight tracks integrity, pH, maybe gender, diagnosis, and age.


Model<-lme(TZP_Average.260.280~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForGabaGluCqMatrix_QCed3 
# AIC       BIC   logLik
# -4641.313 -4626.861 2325.657
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:  0.02121268 3.903433e-16
# 
# Fixed effects: TZP_Average.260.280 ~ Diagnosis 
# Value   Std.Error DF  t-value p-value
# (Intercept)     2.0293876 0.004207606 66 482.3141  0.0000
# DiagnosisBP    -0.0083395 0.006294673 66  -1.3248  0.1898
# DiagnosisSchiz -0.0075114 0.006215029 66  -1.2086  0.2311
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.668       
# DiagnosisSchiz -0.677  0.453
# 
# Standardized Within-Group Residuals:
#   Min        Q1       Med        Q3       Max 
# -1.137689 -1.137689  0.000000  0.000000  1.137689 
# 
# Number of Observations: 133
# Number of Groups: 69 


pdf("Boxplot_GabaGlu_260280_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.280~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Quality (260/280)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.280~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

Model<-lme(TZP_Average.260.230~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# 
# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForGabaGluCqMatrix_QCed3 
# AIC       BIC   logLik
# -4368.308 -4353.856 2189.154
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:   0.1751941 3.382e-16
# 
# Fixed effects: TZP_Average.260.230 ~ Diagnosis 
# Value  Std.Error DF  t-value p-value
# (Intercept)     2.0621955 0.03455098 66 59.68559  0.0000
# DiagnosisBP    -0.0040031 0.05167652 66 -0.07747  0.9385
# DiagnosisSchiz  0.0262493 0.05102254 66  0.51447  0.6086
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.669       
# DiagnosisSchiz -0.677  0.453
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.3130966 -0.6565483  0.0000000  0.0000000  1.3130966 
# 
# Number of Observations: 133
# Number of Groups: 69 

pdf("Boxplot_GabaGlu_260230_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.230~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Quality (260/230)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_Average.260.230~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

Model<-lme(TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~Diagnosis, random=~1|ID, data=SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
# Linear mixed-effects model fit by maximum likelihood
# Data: SubjectInfo_OrderedForGabaGluCqMatrix_QCed3 
# AIC       BIC   logLik
# -4282.686 -4268.234 2146.343
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)     Residual
# StdDev:   0.3964174 2.737442e-16
# 
# Fixed effects: TZP_BioAnalyzer_rRNA.Ratio..28s.18s. ~ Diagnosis 
# Value  Std.Error DF  t-value p-value
# (Intercept)     2.4247112 0.07526888 66 32.21399  0.0000
# DiagnosisBP     0.0507655 0.11240540 66  0.45163  0.6530
# DiagnosisSchiz -0.0309356 0.11098118 66 -0.27875  0.7813
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.670       
# DiagnosisSchiz -0.678  0.454
# 
# Standardized Within-Group Residuals:
#   Min        Q1       Med        Q3       Max 
# -1.622278  0.000000  0.000000  1.622278  3.244556 
# 
# Number of Observations: 133
# Number of Groups: 69 

pdf("Boxplot_GabaGlu_rRNA28s18s_ByDiagnosis.pdf", width=4, height=6)
boxplot(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="RNA Quality (28s/18s rRNA ratio)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$TZP_BioAnalyzer_rRNA.Ratio..28s.18s~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

#I inverted the Cq values for this figure so that it would more intuitively match the NegDeltaCq results:
pdf("Boxplot_GabaGlu_HK_ByDiagnosis.pdf", width=4, height=6)
boxplot((-1*GabaGlu_MeanHousekeeping_QCed2)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, xlab="Diagnosis", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart((-1*GabaGlu_MeanHousekeeping_QCed2)~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


#I inverted the HK measurements to make the results easier to interpret:
Temp<-data.frame(y=(-1*GabaGlu_MeanHousekeeping_QCed2), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC    logLik
# 146.4161 160.8679 -68.20807
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:  0.09642602 0.3927498
# 
# Fixed effects: y ~ Diagnosis 
# Value  Std.Error DF   t-value p-value
# (Intercept)    -22.517645 0.05997003 66 -375.4816  0.0000
# DiagnosisBP     -0.032604 0.08890183 66   -0.3667  0.7150
# DiagnosisSchiz  -0.160809 0.08776940 66   -1.8322  0.0714
# Correlation: 
#   (Intr) DgnsBP
# DiagnosisBP    -0.675       
# DiagnosisSchiz -0.683  0.461
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.0135143 -0.6112058  0.1941969  0.6951280  2.3035781 
# 
# Number of Observations: 133
# Number of Groups: 69  

#Relationship between housekeeping gene expression and diagnosis persists following outlier removal, although it is only a trend now. Sort of:

car::Anova(Model, type=3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 1.4424e+05  1     <2e-16 ***
#   Diagnosis   3.7308e+00  2     0.1548    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=(-1*GabaGlu_MeanHousekeeping_QCed2), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Block.Weight..g.+Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
#As before, no relationship with block weight. 

#I inverted the HK measurements to make the results easier to interpret:
Temp<-data.frame(y=(-1*GabaGlu_MeanHousekeeping_QCed2), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC   logLik
# 64.40814 197.3642 13.79593
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:   0.1148629 0.1900481
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_BioAnalyzer_RIN +      TZP_Average.RNAConc..ng.uL. + Card 
# Value Std.Error DF    t-value p-value
# (Intercept)                 -26.657824 1.5564616 60 -17.127197  0.0000
# DiagnosisBP                  -0.051885 0.0760858 60  -0.681933  0.4979
# DiagnosisSchiz               -0.092072 0.0697948 60  -1.319176  0.1921
# Age                           0.005535 0.0024432 60   2.265407  0.0271
# GenderF                      -0.024269 0.1055736 60  -0.229876  0.8190
# pH                            0.144729 0.2152186 60   0.672474  0.5039
# Hours.Final                   0.001332 0.0046757 60   0.284878  0.7767
# TZP_BioAnalyzer_RIN           0.320178 0.0603656 60   5.303978  0.0000
# TZP_Average.RNAConc..ng.uL.   0.001518 0.0005076 60   2.990863  0.0040

car::Anova(Model, type=3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 438.3633  1  < 2.2e-16 ***
# Diagnosis                     2.6325  2   0.268144    
# Age                           7.6693  1   0.005617 ** 
# Gender                        0.0790  1   0.778701    
# pH                            0.6758  1   0.411040    
# Hours.Final                   0.1213  1   0.727653    
# TZP_BioAnalyzer_RIN          42.0402  1  8.942e-11 ***
# TZP_Average.RNAConc..ng.uL.  13.3676  1   0.000256 ***
# Card                        239.0453 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Well that leaves things completely ambiguous... The outlier removal makes a difference. Not sure how to interpret that - it means that the results from the final analysis are unlikely to be driven by diagnosis-related differences in housekeeping genes, but it also means that it is possible that diagnosis-related differences in housekeeping genes still exist and/or that we are somehow systematically excluding control subjects based on extreme HK expression, which is a feature that could be biological in nature (maybe???)

Temp<-data.frame(y=(-1*GabaGlu_MeanHousekeeping_QCed2), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
car::Anova(Model, type=3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 506.5553  1  < 2.2e-16 ***
#   Diagnosis                     2.9011  2    0.23444    
# Age                           6.1714  1    0.01298 *  
#   Gender                        0.7771  1    0.37802    
# pH                            1.0605  1    0.30310    
# Hours.Final                   0.0836  1    0.77251    
# TZP_BioAnalyzer_RIN          32.7337  1  1.057e-08 ***
#   TZP_Average.RNAConc..ng.uL.  18.9761  1  1.324e-05 ***
#   CardBlock2                   98.8612  1  < 2.2e-16 ***

#No RNAConc?
Temp<-data.frame(y=(-1*GabaGlu_MeanHousekeeping_QCed2), SubjectInfo_OrderedForGabaGluCqMatrix_QCed3)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)
#Similar results

pdf("Plot_GabaGlu_MeanHK_vs_RIN_finalsample.pdf", width=5, height=5)
plot(y~TZP_BioAnalyzer_RIN, data=Temp, ylab="Mean Housekeeping Gene Expression (-Cq)", xlab="RNA Integrity Number (RIN)")
BestFitLine<-lm(y~TZP_BioAnalyzer_RIN, data=Temp)
abline(BestFitLine)
dev.off()

pdf("Plot_GabaGlu_MeanHK_vs_RNAConc_finalsample.pdf", width=5, height=5)
plot(y~TZP_Average.RNAConc..ng.uL., data=Temp, ylab="Mean Housekeeping Gene Expression (-Cq)", xlab="RNA Concentration (ng/uL)")
BestFitLine<-lm(y~TZP_Average.RNAConc..ng.uL., data=Temp)
abline(BestFitLine)
dev.off()

pdf("Plot_GabaGlu_MeanHK_vs_Age_finalsample.pdf", width=5, height=5)
plot(y~Age, data=Temp, ylab="Mean Housekeeping Gene Expression (-Cq)", xlab="Age (yrs)")
BestFitLine<-lm(y~Age, data=Temp)
abline(BestFitLine)
dev.off()

#Getting a sense of how much of the variation in HK in the final sample is technical:
summary.lm(lm(y~TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL., data=Temp))
# Residual standard error: 0.3572 on 130 degrees of freedom
# Multiple R-squared:  0.2597,	Adjusted R-squared:  0.2483 
# F-statistic:  22.8 on 2 and 130 DF,  p-value: 3.257e-09

summary.lm(lm(y~TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230, data=Temp))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.03610 -0.23553  0.02879  0.23461  0.85399 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          -2.121e+01  3.727e+00  -5.692 8.24e-08 ***
# TZP_BioAnalyzer_RIN                   3.386e-01  6.154e-02   5.501 2.00e-07 ***
# TZP_Average.RNAConc..ng.uL.           1.577e-03  6.105e-04   2.583   0.0109 *  
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s. -1.616e-01  9.483e-02  -1.705   0.0907 .  
# TZP_Average.260.280                  -2.304e+00  1.741e+00  -1.323   0.1882    
# TZP_Average.260.230                   2.039e-01  1.811e-01   1.126   0.2622    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3512 on 127 degrees of freedom
# Multiple R-squared:  0.3007,	Adjusted R-squared:  0.2731 
# F-statistic: 10.92 on 5 and 127 DF,  p-value: 9.346e-09

#Adding additional RNA-metric variables only accounts for another 5% of variation... but that is more than all of the biological co-variates combined... 

#I've never heard about people controlling for these other measurements
#Info about how RNA integrity and RNA purity affects qPCR:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3089491/
#http://gene-quantification.org/fleige-pfaffl-rin-2006.pdf
#https://www.researchgate.net/post/What_is_the_effect_of_the_260_280_value_on_PCR

summary.lm(lm(y~TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2, data=Temp))
# Residual standard error: 0.2848 on 129 degrees of freedom
# Multiple R-squared:  0.5331,	Adjusted R-squared:  0.5222 
# F-statistic: 49.09 on 3 and 129 DF,  p-value: < 2.2e-16

summary.lm(lm(y~TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+CardBlock2, data=Temp))
# Residual standard error: 0.2769 on 126 degrees of freedom
# Multiple R-squared:  0.5688,	Adjusted R-squared:  0.5483 
# F-statistic:  27.7 on 6 and 126 DF,  p-value: < 2.2e-16

summary.lm(lm(y~TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, data=Temp))
# Residual standard error: 0.2792 on 95 degrees of freedom
# Multiple R-squared:  0.6695,	Adjusted R-squared:  0.5408 
# F-statistic: 5.202 on 37 and 95 DF,  p-value: 4.799e-11

summary.lm(lm(y~TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+Card, data=Temp))
# Residual standard error: 0.2737 on 92 degrees of freedom
# Multiple R-squared:  0.6924,	Adjusted R-squared:  0.5586 
# F-statistic: 5.177 on 40 and 92 DF,  p-value: 3.894e-11

#So technical co-variates alone account for 56-70% of the variation in HK expression.

#Biological confounds:
summary.lm(lm(y~Age+Gender+pH+Hours.Final, data=Temp))

# Residual standard error: 0.4087 on 128 degrees of freedom
# Multiple R-squared:  0.04579,	Adjusted R-squared:  0.01597 
# F-statistic: 1.535 on 4 and 128 DF,  p-value: 0.1958

summary.lm(lm(y~Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2, data=Temp))
# Residual standard error: 0.2793 on 125 degrees of freedom
# Multiple R-squared:  0.5647,	Adjusted R-squared:  0.5403 
# F-statistic: 23.16 on 7 and 125 DF,  p-value: < 2.2e-16

summary.lm(lm(y~Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+CardBlock2, data=Temp))
# Residual standard error: 0.2705 on 122 degrees of freedom
# Multiple R-squared:  0.6014,	Adjusted R-squared:  0.5687 
# F-statistic: 18.41 on 10 and 122 DF,  p-value: < 2.2e-16

summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+CardBlock2, data=Temp))
# Residual standard error: 0.2782 on 123 degrees of freedom
# Multiple R-squared:  0.5751,	Adjusted R-squared:  0.544 
# F-statistic:  18.5 on 9 and 123 DF,  p-value: < 2.2e-16

summary.lm(lm(y~Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, data=Temp))
# Residual standard error: 0.2715 on 91 degrees of freedom
# Multiple R-squared:  0.7005,	Adjusted R-squared:  0.5656 
# F-statistic: 5.191 on 41 and 91 DF,  p-value: 3.319e-11

summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card, data=Temp))
# Residual standard error: 0.271 on 89 degrees of freedom
# Multiple R-squared:  0.7083,	Adjusted R-squared:  0.5673 
# F-statistic: 5.025 on 43 and 89 DF,  p-value: 6.769e-11

#So Diagnosis+Biological confounds only change the R-squared by 3-4%, and only ~1% of that is Diagnosis.


#Card makes a big difference!

pdf("Boxplot_GabaGlu_MeanHKexpression_vs_Card_finalsample.pdf", width=20, height=5)
boxplot(y~Card, data=Temp, ylab="Mean Housekeeping Gene Expression", las=2)
dev.off()

pdf("Boxplot_GabaGlu_MeanHKexpression_vs_CardBlock2_finalsample.pdf", width=4, height=5)
boxplot(y~CardBlock2, data=Temp, ylab="Mean Housekeeping Gene Expression", las=2)
dev.off()


table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Cohort)
#           Cohort 11 Cohort 12 Cohort 13 Cohort 7 Cohort 8 Dep Cohort 1 Dep Cohort 5 Dep Cohort 6 Schiz Cohort 2
# Control         9        17        10        0        0            4            2            0              7
# BP              2         8         0       11        4            4            8            4              0
# Schiz           6         4         4        9        8            0            0            2             10

#unfortunately, some of these cohorts don't include controls, so I don't think I can use it easily as a co-variate. :(


table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Diagnosis, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3$Card)

# Card 1.eds Card 10.eds Card 11.eds Card 12.eds Card 13.eds Card 14.eds Card 15.eds Card 16.eds Card 17.eds Card 18.eds
# Control          1           0           2           2           1           1           2           2           1           1
# BP               1           2           1           1           2           2           0           0           1           1
# Schiz            1           2           1           1           1           0           2           2           2           2
# 
# Card 19.eds Card 2.eds Card 20.eds Card 21.eds Card 22.eds Card 23.eds Card 24.eds Card 25.eds Card 26.eds Card 27.eds
# Control           1          2           2           2           2           1           1           2           2           0
# BP                1          1           0           0           0           2           2           1           1           2
# Schiz             1          1           1           1           1           1           1           0           0           2
# 
# Card 28.eds Card 29.eds Card 3.eds Card 30.eds Card 31.eds Card 32.eds Card 33.eds Card 34.eds Card 35.eds Card 36.eds Card 4.eds
# Control           0           2          1           2           2           2           1           1           2           2          1
# BP                2           1          1           1           1           1           2           2           1           1          1
# Schiz             2           1          2           1           1           1           1           1           1           1          2
# 
# Card 5.eds Card 6.eds Card 7.eds Card 8.eds Card 9.eds
# Control          1          1          1          2          0
# BP               1          1          1          1          2
# Schiz            1          1          1          1          2

#*********************************************

