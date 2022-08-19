#**************************************************************

#Before normalizing by housekeeping gene:
#Double-checking whether there is a relationship between housekeeping genes and diagnosis: 

dim(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#[1] 136  48
dim(GabaGlu_Cq_AllSubjects_QCed2)
#[1]  93 136

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$ID)
# 1 10 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 37 38 39  4 40 41 42 44 45 46 47 48 49  5 51 52 53 54 55 56 
# 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  1  2  2  2  2  2  2  2  2  2  1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 
# 57 58 59  6 60 61 62 63 64 65 66 67 68 69  7 70 71 72 73 74 75  8  9 
# 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 


#Running analyses requires an MLM because there are 1-2 samples per ID.  

#This can be done in 2 ways:

#I originally was using the nlme method.


library(nlme)
summary(lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim')))
car::Anova(lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML", contrasts=list(Diagnosis=contr.sum, Gender=contr.sum)), type="III")



for(i in c(83:93)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspects_nlme.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspects_nlme.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#...but then discovered that the nlme method might limit my ability to add other batch-related variables if they were defined as random effects (independent/crossed random effects):
#https://errickson.net/stats-notes/vizrandomeffects.html
#This seems to suggest that it is difficult to specify crossed random effects in nlme and that I should switch to lme4 instead. Sigh:
#https://biostatmatt.com/archives/2718
#...so I switched to lme4-style analysis (below), which work similarly except that there are no p-values. I remember this being a problem last time I used lme4 too, and that the work around for it was weird.
#https://rdrr.io/cran/lme4/man/pvalues.html 
#This led me to choose to calculate p-values using a car package anova summary (type III), which uses a Wald chisquare test to compare the model with a reduced model. 
#Also, one thing that came up while reading was that apparently REML (which is the standard for evaluating multilevel models) shouldn't be used in combination with anova, but anova() function may automatically correct it (make it ML). I tested that out and it was true:
#https://www.researchgate.net/post/Problems_performing_ANOVA_and_Mixed_Linear_model_on_R
#https://www.researchgate.net/post/REML_FALSE_versus_REML_TRUE_lme4_package_in_R-any_thoughts

#but then later discovered that all of this might not actually matter in the end, because it seemed like the addition of any batch effect that was multicollinear with diagnosis in some way as a random variable ended up producing results with the shared variance assigned to the random effect batch variable (type 1 style), so I ended up only including dissection batch as a fixed effect and then later exploring other ways of controlling for card by averaging across replicate cards following controlling for housekeeping gene expression.

library(lme4)


#Collapsed - this should probably just be moved to the end and/or eliminated (since we have a nice sensitivity analysis later that includes HK Cq measurements)

#The usual suspects:

for(i in c(83:93)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspects.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspects.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}

#Sanity check: This output looks almost exactly the same as the output from nlme. (*phew*)
#Some of the housekeeping genes show relationships with age, pH
#Three genes show relationships with diagnosis: TBP (p=0.0192) and TFRC (p=0.0124), IPO8 (p=0.0277), B2M and GUSB  have non significant trends

for(i in c(83:93)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndDissection.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndDissection.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Dissection group is sig for about half the genes

for(i in c(83:93)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndCard.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndCard.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Card seems to be really important as a technical co-variate

for(i in c(83:93)){
  Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
  Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndCardRINConc.txt")
  stats_output <- c(
    print(row.names(GabaGlu_Cq_AllSubjects_QCed2)[i]),
    capture.output(car::Anova(Model, type="III")),
    print("******************************************************************************")
  )
  cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndCardRINConc.txt", sep="\n", append=TRUE)
  close(OutputtingStats)
  rm(stats_output)
  rm(Temp)
  rm(Model)
}
#Noteworthy: RNA concentration, RIN, and Card are significant co-variates for almost every gene.
#Age is sometimes a significant co-variate, pH is for one gene.
#There is a relationship with Diagnosis for GUSB (p=0.016828), HMBS (p=0.039774), and TFRC (p=0.0003418), a trend for IPO8 (p=0.06016) and TBP (p=0.061195)



# for(i in c(83:93)){
#   Temp<-data.frame(y=GabaGlu_Cq_AllSubjects_QCed2[i,], SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#   Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group+Card + (1 | ID), data = Temp, REML=F)
#   OutputtingStats<-file("GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndDissectionCard.txt")
#   stats_output <- c(
#     print(row.names(GabaGlu_Cq_AllSubjects_QCed2)[i]),
#     capture.output(car::Anova(Model, type="III")),
#     print("******************************************************************************")
#   )
#   cat(stats_output, file="GabaGlu_MLM_HousekeepingGenesByTheUsualSuspectsAndDissectionCard.txt", sep="\n", append=TRUE)
#   close(OutputtingStats)
#   rm(stats_output)
#   rm(Temp)
#   rm(Model)
# }
#This model is rank deficient (36 coefficients just for card + 18 coefficients for all other fixed effects) - so it is not recommended that I include both dissection group and card. :(

table(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group)
#Ah-they are often partially redundant. Each card tends to contain only samples from 1-2 dissection groups.
#So card mostly covers variation due to dissection group, and variation due to dissection group is partially exprssed in variation by card.

table(paste(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group))

# Card 1.eds 1  Card 10.eds 3  Card 10.eds 4  Card 11.eds 4  Card 12.eds 4  Card 13.eds 5 
# 4              2              2              4              4              4 
# Card 14.eds 5  Card 15.eds 5  Card 15.eds 6  Card 16.eds 5  Card 16.eds 6  Card 17.eds 6 
# 3              2              2              2              2              4 
# Card 18.eds 6  Card 19.eds 7   Card 2.eds 1  Card 20.eds 7  Card 21.eds 7  Card 21.eds 8 
# 4              4              4              3              2              2 
# Card 22.eds 7  Card 22.eds 8  Card 23.eds 8  Card 24.eds 8  Card 25.eds 9  Card 26.eds 9 
# 2              2              4              4              4              4 
# Card 27.eds 10  Card 27.eds 9 Card 28.eds 10  Card 28.eds 9 Card 29.eds 10   Card 3.eds 1 
# 2              2              2              2              4              2 
# Card 3.eds 2 Card 30.eds 10 Card 31.eds 11 Card 32.eds 11 Card 33.eds 11 Card 33.eds 12 
# 2              4              4              4              2              2 
# Card 34.eds 11 Card 34.eds 12 Card 35.eds 12 Card 36.eds 12   Card 4.eds 1   Card 4.eds 2 
# 2              2              4              4              2              2 
# Card 5.eds 12   Card 5.eds 2  Card 6.eds 12   Card 6.eds 2   Card 7.eds 3   Card 8.eds 3 
# 1              3              1              2              4              4 
# Card 9.eds 3   Card 9.eds 4 
# 2              2 

length(table(paste(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group)))
#[1] 52
#So just a little bit more than 13 (# of dissection groups) + 36 (# of cards)=49
#There is definitely no way to include both - the subsetted sample sizes are just too puny.
#Card matters strongly for every single gene, so I'm inclined to use the model with card as a co-variate and drop dissection group.
#Alternatively, we could try to see if the results in some of the dissection groups or cards are correlated and collapse them in some way.
#Or do some sort of model selection procedure + test/train?
#Good description of the problem of rank decificiency:
#https://stats.stackexchange.com/questions/35071/what-is-rank-deficiency-and-how-to-deal-with-it/35077#35077


#########################


#Calculating mean housekeeping gene expression


#What about either the mean or trimmed mean?

GabaGlu_TrimmedMeanHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed2[c(83:93),], 2, function(y) mean(y, trim=0.20, na.rm=T))
hist(GabaGlu_TrimmedMeanHousekeeping, breaks=20)

GabaGlu_TrimmedMeanNonHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed2[-c(83:93),], 2, function(y) mean(y, trim=0.20, na.rm=T))
hist(GabaGlu_TrimmedMeanNonHousekeeping, breaks=20)

GabaGlu_MeanHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed2[c(83:93),], 2, function(y) mean(y, na.rm=T))
hist(GabaGlu_MeanHousekeeping, breaks=20)

GabaGlu_MeanNonHousekeeping<-apply(GabaGlu_Cq_AllSubjects_QCed2[-c(83:93),], 2, function(y) mean(y, na.rm=T))
hist(GabaGlu_MeanNonHousekeeping, breaks=20)

plot(GabaGlu_MeanNonHousekeeping~GabaGlu_MeanHousekeeping)
cor(GabaGlu_MeanNonHousekeeping,GabaGlu_MeanHousekeeping)
#[1] 0.840277
#That's weird - I wonder if it is being driven by the influence of NAs on the mean housekeeping score
plot(GabaGlu_TrimmedMeanNonHousekeeping~GabaGlu_TrimmedMeanHousekeeping)
#Yep - the relationship between trimmed means for housekeeping and non-housekeeping is pretty linear.
cor(GabaGlu_TrimmedMeanNonHousekeeping, GabaGlu_TrimmedMeanHousekeeping)
#[1] 0.8192341
#...but still not as strong as it was for Z-score, which had mean housekeeping calculated following imputation.


#Since there are only 11 genes being used to calculate the average housekeeping score that will be used for normalization, missing data from any particular gene (especially if it has much higher or lower expression than average) could throw off the normalization.
#Note: this originally came on my radar because an earlier version of the analysis still included 18S (which has a lot of missing data)

#How many housekeeping genes are missing data?

sum(is.na(GabaGlu_Cq_AllSubjects_QCed2[c(83:93),]))
#[1] 10
sum(is.na(GabaGlu_Cq_AllSubjects_QCed2[-c(83:93),]))
#[1]  9
#Interesting - a disproportionate number of the NAs are housekeeping genes. Is it a particular gene that is problematic?
apply(GabaGlu_Cq_AllSubjects_QCed2[c(83:93),], 1, function(y) sum(is.na(y)))
# HPRT1  GUSB  ACTB   B2M  HMBS  IPO8  PGK1 RPLP0   TBP  TFRC   UBC 
# 0     2     0     0     4     0     2     0     0     0     2 

#Maybe we could impute for just the housekeeping genes? I would rather not impute for the genes that we are actually interested in...



# GabaGlu_Cq_AllSubjects_QCed3<-GabaGlu_Cq_AllSubjects_QCed2[GabaGlu_NumberNA_PerGene<20,]
# dim(GabaGlu_Cq_AllSubjects_QCed3)
# #[1]  93 141


#Replacing NA values with imputed values based on the mean for that gene:
#I will still use the non-imputed dataset for the actual genes of interest analyses, but I will use the mean housekeeping gene value calculated from the imputed dataset to calculate -deltaCq
GabaGlu_Cq_AllSubjects_QCed3_Imputed<-GabaGlu_Cq_AllSubjects_QCed2
head(GabaGlu_Cq_AllSubjects_QCed3_Imputed)
for(i in c(1:93)){
  GabaGlu_Cq_AllSubjects_QCed3_Imputed[i, is.na(GabaGlu_Cq_AllSubjects_QCed3_Imputed[i,])]<-mean(GabaGlu_Cq_AllSubjects_QCed3_Imputed[i,], na.rm=T)
}

#Sanity check:
sum(is.na(GabaGlu_Cq_AllSubjects_QCed3_Imputed))
#[1] 0
sum((apply(GabaGlu_Cq_AllSubjects_QCed2, 1, function(y) mean(y, na.rm=T)))==(apply(GabaGlu_Cq_AllSubjects_QCed3_Imputed, 1, function(y) mean(y, na.rm=T))))
#[1] 93
#So the mean for each gene hasn't shifted.

GabaGlu_TrimmedMeanHousekeeping_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:93),], 2, function(y) mean(y, trim=0.20, na.rm=T))
hist(GabaGlu_TrimmedMeanHousekeeping_QCed, breaks=20)

GabaGlu_TrimmedMeanNonHousekeeping_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed3_Imputed[-c(83:93),], 2, function(y) mean(y, trim=0.20, na.rm=T))
hist(GabaGlu_TrimmedMeanNonHousekeeping_QCed, breaks=20)

GabaGlu_MeanHousekeeping_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed3_Imputed[c(83:93),], 2, function(y) mean(y, na.rm=T))
hist(GabaGlu_MeanHousekeeping_QCed, breaks=20)

GabaGlu_MeanNonHousekeeping_QCed<-apply(GabaGlu_Cq_AllSubjects_QCed3_Imputed[-c(83:93),], 2, function(y) mean(y, na.rm=T))
hist(GabaGlu_MeanNonHousekeeping_QCed, breaks=20)

#Sanity check:
hist(GabaGlu_MeanHousekeeping_QCed)
#Centered around 22.25
hist(GabaGlu_MeanNonHousekeeping_QCed)
#Centered around 23.25
#So typically housekeeping genes are more highly expressed than non-housekeeping genes.


pdf("Plot_GabaGlu_MeanNonHK_vs_MeanHK.pdf", width=5, height=5)
plot(GabaGlu_MeanNonHousekeeping_QCed~GabaGlu_MeanHousekeeping_QCed)
abline(a=0, b=1, col="red")
dev.off()

cor(GabaGlu_MeanNonHousekeeping_QCed,GabaGlu_MeanHousekeeping_QCed)
#[1] 0.8872947
plot(GabaGlu_TrimmedMeanNonHousekeeping_QCed~GabaGlu_TrimmedMeanHousekeeping_QCed)
cor(GabaGlu_TrimmedMeanNonHousekeeping_QCed, GabaGlu_TrimmedMeanHousekeeping_QCed)
#[1] 0.8696665
#Conclusion: Both of these are better. The regular mean appears to perform slightly worse than trimmed mean.
#Also, the correlation mostly seems to fall apart for the 5 samples with the lowest Cq (highest expression):
GabaGlu_MeanHousekeeping_QCed[GabaGlu_MeanHousekeeping_QCed<21.75]
# 1       37       49       54       16 
# 21.44739 21.03345 21.71168 21.68564 21.59409
#The 5 particularly high housekeeping gene measurements come from 5 subjects, 4 of which are controls. Snooping through the metadata, there doesn't seem to be anything particularly striking about these subjects/RNA. The samples are from 5 different cards.
#The sample from 37, 1, and 16 are particularly far off the Housekeeping vs. NonHousekeeping line.
SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[GabaGlu_MeanHousekeeping_QCed<21.75,]
# ID        Card Subject.Number  Barcode         Cohort Diagnosis Age Gender   pH AFS Hours.Final Slab.Format Slab.Number
# 1      1  Card 1.eds           5000 B007375A      Cohort 11   Control  77      M 6.79   0        16.5        SLAB           1
# 3841  37 Card 19.eds           5041 B008721A      Cohort 12   Control  64      M 6.71   0        20.2        SLAB           1
# 6529  49 Card 25.eds           5287 B008460A      Cohort 13   Control  48      M 6.73   0        28.1        SLAB           1
# 7393  54 Card 27.eds           5272 B009543A      Cohort 13     Schiz  72      M 7.01   0         8.5        SLAB           1
# 12961 16  Card 7.eds           3265 B004251A Schiz Cohort 2   Control  53      M 6.64   0        14.5        SLAB           1

colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)


#####################

#Let's take a peek at the correlation between housekeeping gene expression and other variables:

Temp<-cor(cbind(GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2[, c(7,9,11,16,24:28,38:42)]))
#Similar to PC1, looks like the strongest correlation is with RIN (R=-0.42), followed by RNAconc (R=-0.328)

write.csv(Temp, "GabaGlu_CorMatrix_HK_vs_NumericVar_Cq_QCed2.csv")


#To make plots more intuitive, I inverted the Cq measurements (so that increase = greater expression)

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_RIN_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_RIN, xlab="TZP RNA Integrity Number (RIN)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_RIN)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#HK expression increases with RIN

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_RNAConc_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL., xlab="TZP RNA Concentration (ng/uL)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.RNAConc..ng.uL.)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#HK expression increases with RNA conc - this is interesting because the same *quantity* of RNA was included for each sample in the reaction.

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_260280_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.280, xlab="TZP RNA Purity (260/280)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.280)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#very slight negative relationship, but maybe useful for correcting RNA conc?

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_260230_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.230, xlab="TZP RNA Purity (260/230)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_Average.260.230)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#very slight positive relationship, but maybe useful for correcting RNA conc?

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_28s18s_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s., xlab="TZP RNA Integrity (28s/18s rRNA)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#Very slight positive relationship.


pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_PMI_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Hours.Final, xlab="Post-Mortem Interval (PMI)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Hours.Final)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#A slight decrease with higher PMI - probably mediated by RIN

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_pH_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$pH, xlab="Brain pH", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$pH)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#No relationship with pH

pdf("GabaGlu_Scatterplot_MeanHousekeepingCq_vs_Age_QCed2.pdf", width=6, height=6)
plot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Age, xlab="Age (yrs)", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, cex=1.7)
ForBestFitLine<-lm(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Age)
abline(ForBestFitLine, col="black", lwd=3)
dev.off()
#A slight increase in HK with age.

pdf("GabaGlu_Boxplot_MeanHousekeepingCq_vs_DissectionGroup_QCed2.pdf", width=12, height=6)
boxplot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group, xlab="Dissection/Extraction Group", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey")
stripchart(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Dissecton.Group, xlab="Dissection/Extraction Group", ylab="Mean Housekeeping Gene Expression (-Cq)", vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()


pdf("Boxplot_MeanHousekeepingCq_vs_Card_QCed2.pdf", width=25, height=6)
par(mar=c(10,10,5,5))
boxplot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col="grey", las=2)
stripchart(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Card, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black', las=2)
dev.off()

pdf("GabaGlu_Boxplot_MeanHousekeepingCq_vs_Diagnosis_QCed2.pdf", width=4, height=6)
boxplot(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis, xlab="Diagnosis", ylab="Mean Housekeeping Gene Expression (-Cq)", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,3,4), outline=FALSE)
stripchart(-GabaGlu_MeanHousekeeping_QCed~SubjectInfo_OrderedForGabaGluCqMatrix_QCed2$Diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
dev.off()

#So Bipolar and Schizophrenia are associated with slightly lower housekeeping gene expression (higher Cq)


#Alright, back to looking at the relationship with diagnosis and the usual suspects:

#Just diagnosis:

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+ (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 1.3626e+05  1    < 2e-16 ***
#   Diagnosis   6.3893e+00  2    0.04098 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Diagnosis and the usual suspects
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept) 197.9756  1    < 2e-16 ***
# Diagnosis     7.8419  2    0.01982 *  
# Age           6.3869  1    0.01150 *  
# Gender        0.0477  1    0.82718    
# pH            1.0327  1    0.30952    
# Hours.Final   5.3619  1    0.02058 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Diagnosis and the usual suspects+dissection/extraction
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+Dissecton.Group + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)     189.8426  1    < 2e-16 ***
# Diagnosis         6.2458  2    0.04403 *  
# Age               6.4608  1    0.01103 *  
# Gender            0.0972  1    0.75519    
# pH                0.6575  1    0.41743    
# Hours.Final       2.5156  1    0.11273    
# Dissecton.Group  14.7408 12    0.25592    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Diagnosis and RNA technical variables:

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                          34.9042  1  3.463e-09 ***
# Diagnosis                             4.6921  2   0.095745 .  
# Block.Weight..g.                      0.0407  1   0.840107    
# TZP_Average.RNAConc..ng.uL.           8.3414  1   0.003875 ** 
# TZP_BioAnalyzer_RIN                  27.6420  1  1.460e-07 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.  1.9895  1   0.158395    
# TZP_Average.260.280                   0.4558  1   0.499609    
# TZP_Average.260.230                   1.2608  1   0.261503    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                          2865.1034  1  < 2.2e-16 ***
#   Diagnosis                               4.4538  2     0.1079    
# TZP_Average.RNAConc..ng.uL.            15.7050  1  7.403e-05 ***
#   TZP_BioAnalyzer_RIN                    28.6930  1  8.481e-08 ***
#   TZP_BioAnalyzer_rRNA.Ratio..28s.18s.    3.6401  1     0.0564 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Mostly RNAConc & RIN (maybe 28s/18s)

#How much of the variation is just RNAconc and RIN?
#Note that this quick analysis doesn't take into account the replicate structure of the dataset
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
summary.lm(lm(y~TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN, data = Temp))
# Residual standard error: 0.3898 on 133 degrees of freedom
# Multiple R-squared:  0.2502,	Adjusted R-squared:  0.2389 
# F-statistic: 22.19 on 2 and 133 DF,  p-value: 4.838e-09

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
summary.lm(lm(y~TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s., data = Temp))
# Residual standard error: 0.3858 on 132 degrees of freedom
# Multiple R-squared:  0.2708,	Adjusted R-squared:  0.2542 
# F-statistic: 16.34 on 3 and 132 DF,  p-value: 4.336e-09
#28s/18s doesn't add much to the model.

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN, data = Temp))
# Residual standard error: 0.3833 on 127 degrees of freedom
# Multiple R-squared:  0.3076,	Adjusted R-squared:  0.264 
# F-statistic: 7.054 on 8 and 127 DF,  p-value: 1.102e-07
#Biological variables don't add much either - a little bit though.
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN, data = Temp))
# Residual standard error: 0.3806 on 130 degrees of freedom
# Multiple R-squared:  0.3013,	Adjusted R-squared:  0.2744 
# F-statistic: 11.21 on 5 and 130 DF,  p-value: 5.37e-09
#And mostly that is age and diagnosis...

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
summary.lm(lm(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+Card, data = Temp))
# Residual standard error: 0.2819 on 92 degrees of freedom
# Multiple R-squared:  0.7287,	Adjusted R-squared:  0.6019 
# F-statistic: 5.746 on 43 and 92 DF,  p-value: 1.152e-12
#There is a huge jump in R-Squared after adding Card (admittedly, it also uses up a ton of df but the jump is also seen in adjusted R-Squared)

#Technical variables and card:
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Block.Weight..g.+TZP_Average.RNAConc..ng.uL.+TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+TZP_Average.260.280+TZP_Average.260.230+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                           48.2518  1  3.748e-12 ***
# Diagnosis                              7.7987  2   0.020255 *  
# Block.Weight..g.                       0.4659  1   0.494886    
# TZP_Average.RNAConc..ng.uL.            7.2758  1   0.006989 ** 
# TZP_BioAnalyzer_RIN                   39.2166  1  3.793e-10 ***
# TZP_BioAnalyzer_rRNA.Ratio..28s.18s.   2.4277  1   0.119210    
# TZP_Average.260.280                    0.3686  1   0.543793    
# TZP_Average.260.230                    0.9288  1   0.335170    
# Card                                 238.0843 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Most important technical variables + the usual suspects
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 281.8803  1  < 2.2e-16 ***
# Diagnosis                     3.2143  2  0.2004618    
# Age                           5.3962  1  0.0201807 *  
# Gender                        0.8455  1  0.3578419    
# pH                            0.1207  1  0.7283143    
# Hours.Final                   0.3620  1  0.5473808    
# TZP_Average.RNAConc..ng.uL.  12.2456  1  0.0004663 ***
# TZP_BioAnalyzer_RIN          20.7274  1  5.295e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Diagnosis and the usual suspects+most important technical variables+dissection/extraction
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+Dissecton.Group + (1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 268.3078  1  < 2.2e-16 ***
# Diagnosis                     2.0255  2   0.363219    
# Age                           2.6371  1   0.104393    
# Gender                        0.5842  1   0.444690    
# pH                            0.0252  1   0.873917    
# Hours.Final                   0.0019  1   0.964949    
# TZP_Average.RNAConc..ng.uL.   8.6286  1   0.003309 ** 
# TZP_BioAnalyzer_RIN          24.8137  1  6.315e-07 ***
# Dissecton.Group              13.3831 12   0.341820    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Most important technical variables + the usual suspects+Card
Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                 415.2403  1  < 2.2e-16 ***
# Diagnosis                     5.0925  2  0.0783761 .  
# Age                          12.4710  1  0.0004133 ***
# Gender                        0.1661  1  0.6835771    
# pH                            0.0591  1  0.8079577    
# Hours.Final                   0.0314  1  0.8593294    
# TZP_Average.RNAConc..ng.uL.  14.8520  1  0.0001163 ***
# TZP_BioAnalyzer_RIN          42.0730  1  8.793e-11 ***
# Card                        245.9820 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lme(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+Card, random=~1|ID, data=Temp, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: Temp 
# AIC      BIC   logLik
# 78.12503 212.1072 6.937484
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept)  Residual
# StdDev:  0.09855401 0.2098339
# 
# Fixed effects: y ~ Diagnosis + Age + Gender + pH + Hours.Final + TZP_Average.RNAConc..ng.uL. +      TZP_BioAnalyzer_RIN + Card 
#                                 Value Std.Error DF   t-value p-value
# (Intercept)                 25.826781 1.5408874 60 16.760978  0.0000
# DiagnosisBP                  0.085900 0.0755761 60  1.136609  0.2602
# DiagnosisSchiz               0.125049 0.0693702 60  1.802637  0.0765
# Age                         -0.006938 0.0023885 60 -2.904695  0.0051
# GenderF                      0.035433 0.1056904 60  0.335253  0.7386
# pH                          -0.042555 0.2128488 60 -0.199930  0.8422
# Hours.Final                  0.000677 0.0046444 60  0.145776  0.8846
# TZP_Average.RNAConc..ng.uL. -0.001562 0.0004928 60 -3.169869  0.0024
# TZP_BioAnalyzer_RIN         -0.320876 0.0601432 60 -5.335205  0.0000
#I didn't copy and paste all of the card values..
#Also, note that the betas are still inverted (Cq).

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+ TZP_BioAnalyzer_RIN+TZP_BioAnalyzer_rRNA.Ratio..28s.18s.+Block.Weight..g.+TZP_Average.260.280+TZP_Average.260.230+Card+(1 | ID), data = Temp, REML=F)
car::Anova(Model, type="III")
#Again, the other technical variables don't add much and just use up df
#260/280 and 260/230  should hypothetically correct for bias in RNAconc, but don't seem to be doing much.
#I guess that makes sense, since 260/280 and 260/230 are strongly multicollinear


Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#TZP_Average.RNAConc..ng.uL._Uncontaminated2<-lm(Temp$TZP_Average.RNAConc..ng.uL.~Temp$TZP_Average.260.280+Temp$TZP_Average.260.230)$residuals+mean(Temp$TZP_Average.RNAConc..ng.uL.)
#TZP_Average.RNAConc..ng.uL._Uncontaminated2<-lm(Temp$TZP_Average.RNAConc..ng.uL.~Temp$TZP_Average.260.280)$residuals+mean(Temp$TZP_Average.RNAConc..ng.uL.)
TZP_Average.RNAConc..ng.uL._Uncontaminated2<-Temp$TZP_Average.RNAConc..ng.uL./Temp$TZP_Average.260.280
Temp2<-data.frame(Temp,TZP_Average.RNAConc..ng.uL._Uncontaminated2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL._Uncontaminated2+TZP_BioAnalyzer_RIN+Card + (1 | ID), data = Temp2, REML=F)
car::Anova(Model, type="III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                                 418.7911  1  < 2.2e-16 ***
# Diagnosis                                     5.3950  2  0.0673752 .  
# Age                                          12.8759  1  0.0003328 ***
# Gender                                        0.1415  1  0.7068199    
# pH                                            0.0359  1  0.8496763    
# Hours.Final                                   0.0412  1  0.8391314    
# TZP_Average.RNAConc..ng.uL._Uncontaminated2  15.5924  1  7.857e-05 ***
# TZP_BioAnalyzer_RIN                          42.3715  1  7.548e-11 ***
# Card                                        246.3118 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#That does seem like a minor improvement.
#Although apparently it doesn't work if I use the lm version (stronger correction) instead of the ratio

#Hmmm... another possibility: 
#I believe the protocol is actually set up to use exactly 1 ug of RNA
#Is it possible that the issue isn't the concentration, but the miscalculation of concentration (based on contamination) which was then used to calculate the 1 ug for the protocol?

Temp<-data.frame(y=GabaGlu_MeanHousekeeping_QCed, SubjectInfo_OrderedForGabaGluCqMatrix_QCed2)
#TZP_Average.RNAConc..ng.uL._Uncontaminated2<-lm(Temp$TZP_Average.RNAConc..ng.uL.~Temp$TZP_Average.260.280+Temp$TZP_Average.260.230)$residuals+mean(Temp$TZP_Average.RNAConc..ng.uL.)
TZP_Average.RNAConc..ng.uL._Uncontaminated2<-Temp$TZP_Average.RNAConc..ng.uL./Temp$TZP_Average.260.280
Temp2<-data.frame(Temp,TZP_Average.RNAConc..ng.uL._Uncontaminated2)
Model<-lmer(y~Diagnosis+Age+Gender+pH+Hours.Final+TZP_Average.RNAConc..ng.uL.+TZP_Average.RNAConc..ng.uL._Uncontaminated2+TZP_BioAnalyzer_RIN+Card + (1 | ID), data = Temp2, REML=F)
car::Anova(Model, type="III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: y
# Chisq Df Pr(>Chisq)    
# (Intercept)                                 370.5918  1  < 2.2e-16 ***
# Diagnosis                                     7.0952  2  0.0287931 *  
# Age                                          14.8468  1  0.0001166 ***
# Gender                                        0.0611  1  0.8048121    
# pH                                            0.0490  1  0.8248703    
# Hours.Final                                   0.1974  1  0.6568548    
# TZP_Average.RNAConc..ng.uL.                   2.0527  1  0.1519359    
# TZP_Average.RNAConc..ng.uL._Uncontaminated2   2.6814  1  0.1015291    
# TZP_BioAnalyzer_RIN                          44.1580  1  3.029e-11 ***
# Card                                        248.5798 35  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#I think I need to come back and think about this some more. I'm not sure that last analysis even made sense...


#Conclusions: 

#Diagnosis shows a slight relationship with average housekeeping gene expression in most models. That certainly complicates things.

#Card seems to really matter as a co-variate, dissection group less so. RIN, RNA concentration, and Age seem like important co-variates too.


#########################################################
