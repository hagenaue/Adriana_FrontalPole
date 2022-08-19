#############################################


##### At this point I focused exclusively on the Gaba Glu cards for quite a while:

#QC 4. How much do the replicate samples correlate?

#There are 96 measurements per sample on each card.  All cards should follow the same order of targets.
#Sanity Check:
temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="1",]
cbind(temp$GeneSymbol[c(1:96)],temp$GeneSymbol[c(97:192)])
#Yep.
plot(temp$Cq[c(1:96)], temp$Cq[c(97:192)])
#Ooh - nice correlation.
#Although it has one super high leverage point-18S has very low Cq in both samples (highly expressed) 

#Percentage difference between the cards for each target:
(abs(temp$Cq[c(1:96)]-temp$Cq[c(97:192)])/apply(cbind(temp$Cq[c(1:96)],temp$Cq[c(97:192)]), 1, mean, na.rm=T))*100

#Let's loop it:

length(names(table(Concatenated_GabaGlu_wSubjectInfo$ID)))
#[1] 72

ReplicateCor_ForEachSubjectID_GabaGlu<-matrix(0, 72, 1)
row.names(ReplicateCor_ForEachSubjectID_GabaGlu)<-names(table(Concatenated_GabaGlu_wSubjectInfo$ID))

Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<-rep(0, length(Concatenated_GabaGlu_wSubjectInfo[,1]))

for (i in c(1:72)){
  
  print(names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i])
  
  temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i],] 
  
  ReplicateCor_ForEachSubjectID_GabaGlu[i]<-cor(x=temp$Cq[c(1:96)], y=temp$Cq[c(97:192)], use="pairwise.complete.obs")
  
  PercentageDiff_BetweenReplicates<-(abs(temp$Cq[c(1:96)]-temp$Cq[c(97:192)])/apply(cbind(temp$Cq[c(1:96)],temp$Cq[c(97:192)]), 1, mean, na.rm=T))*100
  print(length(PercentageDiff_BetweenReplicates))
  print(length(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i]] ))
  Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_GabaGlu_wSubjectInfo$ID==names(table(Concatenated_GabaGlu_wSubjectInfo$ID))[i]] <-rep(PercentageDiff_BetweenReplicates, 2)
  
  rm(temp)
  rm(PercentageDiff_BetweenReplicates)
  
}

#Original notes - this problem is now fixed:
#Ah ha - I found a problem. "Card 22.eds" isn't in my PCR data folder. Interesting. That means two things:
##1) My decoder is wrong regarding which subject IDs are found on which cards (a minor point, since I didn't end up using it, just curious).
##2) I'm missing a file.  I double-checked the server and it isn't there either - so not a download issue on my end. E-mailed Adriana.
##3) She re-uploaded the the missing file, and I re-ran the analysis on 08-25-2020 



pdf("GabaGlu_Histogram_PercentageDiffBetweenReplicates.pdf", width=4, height=4)
hist(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates, xlab="Per Gene: % Diff Between Replicates", main="", breaks=100)
dev.off()

length(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates)
#[1] 13824

sum(is.na(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates)==F)
#[1] 13256

sum(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<20, na.rm=T)
#[1] 13208

13208/13256
#[1] 0.996379
#More than 99% of the measurements match between replicates within 20%. :)
1-0.996379
#0.003621

table(Concatenated_GabaGlu_wSubjectInfo$ID, Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<20)
# FALSE TRUE
# 1      6  182
# 10     0  190
# 12     2  186
# 13     4  184
# 14     2  188
# 15     2  188
# 16     0  190
# 17     0  188
# 18     0  190
# 19     0  188
# 2      0  188
# 20     0  190
# 21     0  188
# 22     0  188
# 23     0  188
# 24     0  190
# 25     0  188
# 26     0  188
# 27     2    0
# 28     0  188
# 29     0  188
# 3      0  188
# 30     0  190
# 31     0  190
# 32     0  190
# 33     0  188
# 34     0  188
# 37     0  190
# 38     4    0
# 39     0  190
# 4      0  188
# 40     2  186
# 41     2  188
# 42     2  186
# 43     0  190
# 44     2  188
# 45     0  188
# 46     0  190
# 47     0  188
# 48     0  188
# 49     6  184
# 5      0  190
# 50     0  192
# 51     2  188
# 52     2  186
# 53     0  190
# 54     0  192
# 55     0  190
# 56     0  190
# 57     0  192
# 58     0  190
# 59     0  192
# 6      0  188
# 60     0  190
# 61     2  188
# 62     0  188
# 63     0  190
# 64     0  190
# 65     2  190
# 66     0  192
# 67     0  188
# 68     0  188
# 69     2  188
# 7      0  188
# 70     0  188
# 71     0  190
# 72     0  190
# 73     0  190
# 74     0  190
# 75     0  188
# 8      0  188
# 9      2  188

hist(table(Concatenated_GabaGlu_wSubjectInfo$ID, Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<20)[,1])

sum(table(Concatenated_GabaGlu_wSubjectInfo$ID, Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<20)[,1]>0)
sum(table(Concatenated_GabaGlu_wSubjectInfo$ID, Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<20)[,1]>2)

#Very few bad measurements per subject (except for #27 & #38, which we already know have a bad sample)

write.csv(table(Concatenated_GabaGlu_wSubjectInfo$ID, Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates<20), "GabaGlu_Table_IDvsGenesWithPercentageDiffMoreThan20Percent.csv")

write.csv(ReplicateCor_ForEachSubjectID_GabaGlu, "ReplicateCor_ForEachSubjectID_GabaGlu.csv")

pdf("GabaGlu_Boxplot_CorrBetweenReplicates.pdf", width=3, height=4)
boxplot(ReplicateCor_ForEachSubjectID_GabaGlu[,1], ylab="Correlation between replicate samples (R)")
dev.off()

summary((ReplicateCor_ForEachSubjectID_GabaGlu[,1]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.8652  0.9922  0.9949  0.9897  0.9962  1.0000       1 

# #Pretty tight. Whiskers & IQ all fall above 0.98 (maybe more like 0.985?). There are 10 "outlier" datapoints outside of that range. Let's see who they are:
ReplicateCor_ForEachSubjectID_GabaGlu[ReplicateCor_ForEachSubjectID_GabaGlu[,1]<0.984,]
# 1        13        15        16        25      <NA>        37        49        61        69         7 
# 0.8651958 0.9583229 0.9795163 0.9715362 0.9652474        NA 0.9737637 0.9398720 0.9773817 0.9835464 0.9802135 

#Note: 27 is  NA for correlation because one of its replicates (Card 14.eds) is bad and includes exactly one measurement (18S - the most highly expressed gene) and everything else is NA. Ouch. 
##...and 38 says cor=1, because it only has two measurements in one of its replicates (Card 17.eds).

#Other notes:
#Sample ID #1 has a particularly low replicate sample correlation (0.86)
#49 and 13 also have notably low replicate sample correlations.

#I'd like to quickly snoop at their correlation plots to see if those results are driven by outlier measurements:
temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="1",]
plot(temp$Cq[c(1:96)], temp$Cq[c(97:192)])
#There are 3 measurements that are pretty far off the trend line, and another 6 that are also pretty noticeable.
#This sample also had the highest number of genes with percentage differences between replicates >20% (6/2=3)

temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="49",]
plot(temp$Cq[c(1:96)], temp$Cq[c(97:192)])
#This sample also has 2-3 outlier measurements, but hte rest of the measurements are pretty tight on the line.
#This sample also had the highest number of genes with percentage differences between replicates >20% (6/2=3)

temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="13",]
plot(temp$Cq[c(1:96)], temp$Cq[c(97:192)])
#2 outlier measurements, and everything else correlates beautifully.

#Based on these plots and the small number of replicate measurements with greater than 20% percentage difference, I'm inclined to just replace those few questionable measurements with NA rather than mark these samples as being particularly low quality.
sum(is.na(Concatenated_GabaGlu_wSubjectInfo$Cq))
#[1] 322
sum(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates>20, na.rm=T)
#[1] 48

plot(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates, Concatenated_GabaGlu_wSubjectInfo$Cq) 
#Need to remove the Card 14 ID 27 samples and Card 20 ID 38 samples. Also axis may need log-scaling
plot(log(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates,2), Concatenated_GabaGlu_wSubjectInfo$Cq) 

Concatenated_GabaGlu_wSubjectInfo$Cq[(Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")]
Concatenated_GabaGlu_wSubjectInfo$Cq[(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds")]
Concatenated_GabaGlu_wSubjectInfo$Cq[(Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds")]

plot(PercentageDiff_BetweenReplicates~Cq, data=Concatenated_GabaGlu_wSubjectInfo[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE,], cex=0.3, col=rgb(c(0,1),c(0,1),c(0,1),alpha=0.5)) 
#It would be nice to combine this with the data from the DA5HT Card to see what replicates...
#The 18S measurements are really weird. Otherwise, there doesn't seem to be much relationship between Cq and accuracy in the other measurements.

Concatenated_GabaGlu_wSubjectInfo$Cq_NoBad<-Concatenated_GabaGlu_wSubjectInfo$Cq
Concatenated_GabaGlu_wSubjectInfo$Cq_NoBad[Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates>20]<-NA
Concatenated_GabaGlu_wSubjectInfo$Cq_NoBad[Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates>20]
#Looks good.

write.csv(Concatenated_GabaGlu_wSubjectInfo, "Concatenated_GabaGlu_wSubjectInfo_AfterQC.csv")

#Let's see how much that improves things for these samples with lower replicate correlations:
temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="1",]
cor(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)], use="pairwise.complete.obs")
#[1] 0.9710771
plot(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)])
#1 still has some measurements that don't correlate fantastically, but that was a huge improvement.

temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="49",]
cor(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)], use="pairwise.complete.obs")
#[1] 0.9901122
plot(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)])
#49 looks beautiful now.
temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="13",]
cor(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)], use="pairwise.complete.obs")
#[1] 0.9963091
plot(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)])
#13 also looks beautiful now.

temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="25",]
cor(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)], use="pairwise.complete.obs")
#[1] 0.9652474
plot(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)])


temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="16",]
cor(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)], use="pairwise.complete.obs")
#[1] 0.9715362
plot(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)])

temp<-Concatenated_GabaGlu_wSubjectInfo[Concatenated_GabaGlu_wSubjectInfo$ID=="37",]
cor(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)], use="pairwise.complete.obs")
#[1] 0.9737637
plot(temp$Cq_NoBad[c(1:96)], temp$Cq_NoBad[c(97:192)])

#It would be helpful to interpret this data in the context of the correlations across subjects (vs. within subjects). Our data frame is currently not well set up for that 

#If I could rearrange it, I could reuse a lot of my code from microarray/RNA-Seq analyses.

dim(Concatenated_GabaGlu_wSubjectInfo)
#[1] 13824    73

13824/96
#[1] 144
144/2
#[1] 72

Concatenated_GabaGlu_wSubjectInfo$ID[c(1:96)]
Concatenated_GabaGlu_wSubjectInfo$ID[c(97:(96+96))]

table(Concatenated_GabaGlu_wSubjectInfo$ID[seq(1, 13824, 96)])
#looks good: two samples per subject

Concatenated_GabaGlu_wSubjectInfo$ID[seq(1, 13824, 96)]

GabaGlu_Cq_AllSubjects<-matrix(0,96,144)
colnames(GabaGlu_Cq_AllSubjects)<-Concatenated_GabaGlu_wSubjectInfo$ID[seq(1, 13824, 96)]
colnames(GabaGlu_Cq_AllSubjects)
#If we use the subject IDs, the columnames sometimes repeat - I thought R would be unhappy about it, but I don't see a .2 added to the names. Interesting.

FirstGeneForSample_iterative<-seq(1, 13824, 96)
LastGeneForSample_iterative<-seq(96, 13824, 96)

i<-1
#Sanity check:
Concatenated_GabaGlu_wSubjectInfo$ID[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]
row.names(GabaGlu_Cq_AllSubjects)<-Concatenated_GabaGlu_wSubjectInfo$GeneSymbol[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]

#This has been updated from the original version of the analysis, which included the Cq measurements that were >20% different from their replicate.
for(i in c(1:144)){
  GabaGlu_Cq_AllSubjects[,i]<-Concatenated_GabaGlu_wSubjectInfo$Cq_NoBad[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]
}

head(GabaGlu_Cq_AllSubjects)
# 1      2      3      4     17     18     19     20     21     22     23     24     21     22     23     24     25     26
# ABAT    20.513 20.655 20.957 20.553 21.166 21.669 20.939 21.598 21.161 20.989 21.125 20.875 20.683 22.012 21.678 20.426 21.716 21.090
# ADCY7   25.734 24.956 25.178 25.612 26.201 26.061 25.861 26.895 25.576 24.902 25.002 25.145 25.469 25.750 25.353 25.210 26.383 25.021
# ADORA1  23.737 23.150 23.728 23.443 23.945 24.597 24.181 24.652 24.189 23.408 23.796 23.617 24.012 24.702 24.492 23.849 24.586 23.241
# ADORA2A 27.270 26.826 26.960 26.980 27.683 27.646 27.191 28.035 27.333 26.829 26.842 26.795 27.148 28.093 27.573 26.906 27.727 26.831
# ALDH5A1 22.833 22.772 23.126 22.469 23.415 23.650 23.139 23.605 23.060 23.085 23.426 22.841 22.605 23.769 23.818 22.326 23.894 23.223
# SST     22.400 22.200 21.470 22.466 23.223 23.710 24.028 25.061 22.427 21.856 22.119 21.505 22.285 22.924 22.867 21.866 23.407 23.036

#yesss....
#colnames 27(v2) and 38 (v2) are mostly NA data.
colnames(GabaGlu_Cq_AllSubjects)[23]
colnames(GabaGlu_Cq_AllSubjects)[50]

colnames(Concatenated_GabaGlu_wSubjectInfo)

SubjectInfo_OrderedForGabaGluCqMatrix<-Concatenated_GabaGlu_wSubjectInfo[FirstGeneForSample_iterative, c(2, 24:69)]
head(SubjectInfo_OrderedForGabaGluCqMatrix)
str(SubjectInfo_OrderedForGabaGluCqMatrix)
# 'data.frame':	144 obs. of  47 variables:
#   $ ID                                  : chr  "1" "2" "3" "4" ...
# $ Card                                : chr  "Card 1.eds" "Card 1.eds" "Card 1.eds" "Card 1.eds" ...
# $ Subject.Number                      : int  5000 5313 2169 4087 3772 4656 4419 4114 5001 3618 ...
# $ Barcode                             : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Cohort                              : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis                           : Factor w/ 3 levels "Control","BP",..: 1 3 1 2 2 3 2 3 1 2 ...
# $ Age                                 : int  77 42 18 51 70 43 51 58 73 36 ...
# $ Gender                              : Factor w/ 2 levels "M","F": 1 1 1 2 1 1 2 2 1 2 ...
# $ pH                                  : num  6.79 6.77 6.97 6.81 6.96 6.6 6.63 6.59 6.95 7 ...
# $ AFS                                 : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final                         : num  16.5 23.7 22 16.8 27.5 14.7 9.2 27.3 30.2 25.5 ...
# $ Slab.Format                         : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number                         : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group                     : Factor w/ 13 levels "1","2","3","4",..: 1 1 1 1 3 3 4 4 4 4 ...
# $ TTP.ID                              : chr  "TTP010001" "TTP010002" "TTP010003" "TTP010004" ...
# $ Block.Weight..g.                    : num  0.72 0.66 0.76 0.56 0.82 0.62 1.13 1.03 0.72 0.9 ...
# $ TZP.ID                              : chr  "TZP010001" "TZP010002" "TZP010003" "TZP010004" ...
# $ TZP_MBNI_Nanodrop_RNAConc..ng.uL.   : num  436 343 394 417 367 ...
# $ TZP_MBNI_Nanodrop_260.280           : num  1.96 2.01 1.99 1.97 2 1.96 1.99 2.05 1.96 2.06 ...
# $ TZP_MBNI_Nanodrop_260.230           : num  1.93 2.18 2.02 2.1 1.8 1.84 2.28 2.14 2.23 1.92 ...
# $ TZP_DSC_Nanodrop_RNAConc..ng.uL.    : num  562 354 451 449 431 ...
# $ TZP_DSC_Nanodrop_260.280            : num  2.14 1.92 2.03 2.03 2.04 2.05 2.03 2.08 2.05 2.09 ...
# $ TZP_DSC_Nanodrop_260.230            : num  1.9 2.06 2.01 2.13 1.78 1.82 2.31 2.15 2.32 1.91 ...
# $ TZP_BioAnalyzer_RIN                 : num  7.7 8.1 7.9 8 8.4 6.6 8.3 8.5 8.4 8.4 ...
# $ TZP_BioAnalyzer_rRNA.Ratio..28s.18s.: num  2.4 2.3 3.1 2 3.3 2.4 2.3 2.4 2.3 3 ...
# $ TZP_Average.RNAConc..ng.uL.         : num  499 349 422 433 399 ...
# $ TZP_Average.260.280                 : num  2.05 1.97 2.01 2 2.02 2.01 2.01 2.07 2.01 2.08 ...
# $ TZP_Average.260.230                 : num  1.92 2.12 2.02 2.12 1.79 1.83 2.3 2.15 2.28 1.92 ...
# $ TTP.ID.1                            : chr  "TTP010001" "TTP010002" "TTP010003" "TTP010004" ...
# $ Initial.Volume..uL.                 : int  720 660 760 560 820 620 1130 1030 720 900 ...
# $ Final.Volume..uL.                   : int  620 560 660 460 720 520 1030 930 620 800 ...
# $ TTP_MBNI_Nanodrop_RNAConc..ng.uL.   : num  528 398 452 523 483 ...
# $ TTP_MBNI_Nanodrop_260.280           : num  2 1.82 1.77 1.96 1.75 2.11 1.77 1.89 1.95 1.89 ...
# $ TTP_MBNI_Nanodrop_260.230           : num  1.87 2.02 1.99 1.99 1.75 2.2 2.2 1.84 2.23 2.01 ...
# $ TTP_DSC_Nanodrop_RNAConc..ng.uL.    : num  631 429 555 644 590 ...
# $ TTP_DSC_Nanodrop_260.280            : num  2.02 1.92 2.01 2.02 2.04 2 2 1.94 2.01 1.93 ...
# $ TTP_DSC_Nanodrop_260.230            : num  1.83 2.04 2.05 1.96 1.7 1.71 2.28 1.76 2.27 1.96 ...
# $ TTP_BioAnalyzer_RIN                 : num  7.4 7.8 7.7 7.7 8 6.7 7.8 8.2 8 8 ...
# $ TTP_BioAnalyzer_rRNA.Ratio..28s.18s.: num  2.2 2.2 2.8 2.1 3.1 3.1 2.1 2.3 2.1 2.5 ...
# $ TTP_Average.RNAConc..ng.uL.         : num  579 414 503 583 537 ...
# $ TTP_Average.260.280                 : num  2.01 1.87 1.89 1.99 1.9 ...
# $ TTP_Average.260.230                 : num  1.85 2.03 2.02 1.98 1.73 ...
# $ MissingPH                           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ MissingExtractionData               : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ LowQualityRNA                       : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ SampleNumber                        : int  1 2 3 4 4 1 2 3 4 1 ...
# $ QC_AmplificationIssue               : chr  "N" "N" "N" "N" ...

#These subjects failed earlier (basic) QC (low RIN, highly contaminated):
SubjectInfo_OrderedForGabaGluCqMatrix$ID[SubjectInfo_OrderedForGabaGluCqMatrix$LowQualityRNA==T]
#[1] "43" "43" "50" "50" "12" "12"

#These are the subjects that had amplification issues in the GABAGLU dataset (and also lots of NA measurements):

#     Subject.Number ID Folder SampleNumber QC_AmplificationIssue
# 52           4301 27     NA           NA                     Y
# 69           4359 38     NA           NA                     Y

#Specifically, samples 27(v2) and 38 (v2) are bad.
#They are these samples:
GabaGlu_Cq_AllSubjects[,c(23,50)]

SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC<-FALSE
SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC[c(23,50)]<-TRUE
SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC

sum((SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC|SubjectInfo_OrderedForGabaGluCqMatrix$LowQualityRNA)==T)
#[1] 8
#8 samples need to be removed due to either basic subject-level or sample-level QC
sum((SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC|SubjectInfo_OrderedForGabaGluCqMatrix$LowQualityRNA)==F)
#[1] 136
8/144
#[1] 0.05555556
# a little over 5% of samples removed. :(

GabaGlu_Cq_AllSubjects_QCed<-GabaGlu_Cq_AllSubjects[,(SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC|SubjectInfo_OrderedForGabaGluCqMatrix$LowQualityRNA)==F]
SubjectInfo_OrderedForGabaGluCqMatrix_QCed<-SubjectInfo_OrderedForGabaGluCqMatrix[(SubjectInfo_OrderedForGabaGluCqMatrix$FailedSampleQC|SubjectInfo_OrderedForGabaGluCqMatrix$LowQualityRNA)==F, ]
#Sanity check:
colnames(GabaGlu_Cq_AllSubjects_QCed)
SubjectInfo_OrderedForGabaGluCqMatrix_QCed$ID
#Yep, bad samples now gone & ID orders match in the two dataframes.

#What about bad probes?

GabaGlu_NumberNA_PerGene<-apply(GabaGlu_Cq_AllSubjects_QCed, 1, function(y) sum(is.na(y)))
hist(GabaGlu_NumberNA_PerGene, breaks=30)
#Looks like there are three genes with really high numbers of NAs

#Which ones are they?
GabaGlu_NumberNA_PerGene

# ABAT   ADCY7  ADORA1 ADORA2A ALDH5A1     SST   PVALB    BDNF CACNA1A CACNA1B     18S  CDK5R1   CALB1    DLG4  GABBR1  GABBR2  GABRA1 
# 0       0       0       0       0       0       0       0       0       0      24       1       0       0       0       0       0 
# GABRA2  GABRA4  GABRA5  GABRA6  GABRB1  GABRB3   GABRD   GABRE  GABRG1  GABRG2  GABRG3   GABRQ  GABRR1  GABRR2    GAD1     GLS    GLUL 
# 0       0       0       0       0       0       0       0       0       0       0       0      68       0       0       0       0 
# GNAI1    GNAQ    GPHN   GRIA1   GRIA2   GRIA3   GRIA4   GRIK1   GRIK2   GRIK4   GRIK5   GRIN1  GRIN2A  GRIN2B  GRIN2C    GRM1    GRM2 
# 0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0 
# GRM3    GRM4    GRM5    GRM6    GRM7    GRM8  HOMER1  HOMER2    GFAP   ITPR1   MAPK1     NSF   P2RX7   PHGDH    AQP4    GJA1   S100B 
# 0       0       0      60       0       0       0       0       0       0       0       0       0       0       0       0       0 
# SHANK2 SLC17A6 SLC17A7 SLC17A8  SLC1A1  SLC1A2  SLC1A3  SLC1A6 SLC32A1 SLC38A1  SLC6A1 SLC6A11 SLC6A12 SLC6A13 SLC7A11    SNCA     SRR 
# 0       0       0       0       0       0       0       0       0       0       0       0       0       2       0       0       6 
# HPRT1    GUSB    ACTB     B2M    HMBS    IPO8    PGK1   RPLP0     TBP    TFRC     UBC 
# 0       2       0       0       4       0       2       0       0       0       2 


#Looks like the culprits are 18S, GABRR1, and GRM6

dim(GabaGlu_Cq_AllSubjects_QCed)
#[1]  96 136

#18S
24/136
#[1] 0.1764706
#So more than 17% of the 18S measurements had either no measurement or replicates with >20% difference.

#GABRR1 
68/136
#[1] 0.5
#50% of the samples had GABRR1 measurements with either no measurement or >20% difference between replicates.

#GRM6:
60/136
#[1] 0.4411765
#Over 44% of the samples had GRM6 measurements with either no measurement or >20% difference between replicates.


#Do the bad genes tend to have particularly high or low expression?

GabaGlu_AverageCq_PerGene<-apply(GabaGlu_Cq_AllSubjects_QCed, 1, function(y) mean(y, na.rm=T))

pdf("GabaGlu_Scatterplot_AverageCq_vs_NumberNAPerGene.pdf", width=5, heigh=5)
plot(GabaGlu_AverageCq_PerGene~GabaGlu_NumberNA_PerGene, ylab="Average Cq per probe", xlab= "# of bad measurements per probe")
dev.off()
#Yeah, the 3 genes with ridiculously high numbers of NA have either really high Cq or really low Cq

#That's pretty useful information. Let's output that along with average Cq as a summary table:
write.csv(cbind(GabaGlu_AverageCq_PerGene,GabaGlu_NumberNA_PerGene), "GabaGlu_AverageCq_NumberNA_PerGene.csv")

dim(GabaGlu_Cq_AllSubjects_QCed)
#[1]  96 136

GabaGlu_Cq_AllSubjects_QCed2<-GabaGlu_Cq_AllSubjects_QCed[GabaGlu_NumberNA_PerGene<20,]
dim(GabaGlu_Cq_AllSubjects_QCed2)
#[1]  93 136


#Determining which indices are the housekeeping genes now:
row.names(GabaGlu_Cq_AllSubjects_QCed)[c(11, 86:96)]
#[1] "18S"   "HPRT1" "GUSB"  "ACTB"  "B2M"   "HMBS"  "IPO8"  "PGK1"  "RPLP0" "TBP"   "TFRC"  "UBC"   
row.names(GabaGlu_Cq_AllSubjects_QCed2)[83:93]
#[1] "HPRT1" "GUSB"  "ACTB"  "B2M"   "HMBS"  "IPO8"  "PGK1"  "RPLP0" "TBP"   "TFRC"  "UBC"  


#Also out of curiousity: which samples still have NA?

#after removing the bad genes:
apply(GabaGlu_Cq_AllSubjects_QCed2, 2, function(y) sum(is.na(y)))

# 1  2  3  4 17 18 19 20 21 22 23 24 21 22 23 24 25 26 27 28 25 26 28 29 30 31 32 29 30 31 32 33 34 74 75 33 34 74 75 37 38 39 40  1  2  3 
# 3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  3  0  0 
# 4 37 39 40 41 42 44 41 42 44 45 46 47 48 45 46 47 48 49 51 52 49 51 52 53 54 55 56 53 54 55 56 57 58 59 60  5  6  7  8 57 58 59 60 61 62 
# 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  3  0  0  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 
# 63 64 61 62 63 64 65 66 67 68 65 66 67 68 69 70 71 72 69 70 71 72  5  6  7  8  9 10 73  9 10 73 13 14 15 16 13 14 15 16 17 18 19 20 
# 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  2  0  0  0  0  0  0  0 

#There are only a few samples with NAs, max=3

#After removing the bad genes:
hist(apply(GabaGlu_Cq_AllSubjects_QCed2, 2, function(y) sum(is.na(y))))

sum(apply(GabaGlu_Cq_AllSubjects_QCed2, 2, function(y) sum(is.na(y)))>0)
#[1] 9
sum(apply(GabaGlu_Cq_AllSubjects_QCed2, 2, function(y) sum(is.na(y)))>1)
#[1] 6
sum(apply(GabaGlu_Cq_AllSubjects_QCed2, 2, function(y) sum(is.na(y)))>2)
#[1] 4
#There are only 4 samples with more than 2 NAs now.

sum(is.na(GabaGlu_Cq_AllSubjects_QCed2))
#[1] 19
length(is.na(GabaGlu_Cq_AllSubjects_QCed2))
#[1] 12648

19/12648
#[1] 0.001502214

############

#Examining the sample-sample correlation matrix:

GabaGlu_SubjectSubjectCorMatrix<-cor(GabaGlu_Cq_AllSubjects_QCed2, use="pairwise.complete.obs")

pdf("GabaGlu_SubjectSubjectCorMatrix.pdf", height=14, width=14)
heatmap(GabaGlu_SubjectSubjectCorMatrix)
dev.off()
#Interesting - some of the replicates most closely resemble each other, but certainly not all (maybe 30%?). There are correlation blocks (batch effects??). So the data is definitely noisy.

plot(GabaGlu_Cq_AllSubjects_QCed2[,1]~GabaGlu_Cq_AllSubjects_QCed2[,2])
plot(GabaGlu_Cq_AllSubjects_QCed2[,4]~GabaGlu_Cq_AllSubjects_QCed2[,5])

# I should redo the correlation analysis after z-scoring the data first so that the relative values across subjects are emphasized instead of higher leverage points.

GabaGlu_Cq_AllSubjects_Scaled_QCed<-t(scale(t(GabaGlu_Cq_AllSubjects_QCed2), center=T, scale=T))
head(GabaGlu_Cq_AllSubjects_Scaled_QCed)
#Sanity check:
apply(GabaGlu_Cq_AllSubjects_Scaled_QCed, 1, function(y) mean(y, na.rm=T))
apply(GabaGlu_Cq_AllSubjects_Scaled_QCed, 1, function(y) sd(y, na.rm=T))
#Looks good.

GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed<-cor(GabaGlu_Cq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")
GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed

pdf("GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed.pdf", height=14, width=14)
heatmap(GabaGlu_SubjectSubjectCorMatrix_Scaled_QCed)
dev.off()
#Again - some of the replicates most closely resemble each other, but certainly not all (maybe 50% now?). There are correlation blocks (batch effects??).


#Dealing with the remaining NA values and looking at overall variation in Cq values across samples:

#Just to confirm that there still aren't many NA values:
sum(is.na(GabaGlu_Cq_AllSubjects_Scaled_QCed))
#[1] 19
length(GabaGlu_Cq_AllSubjects_Scaled_QCed)
#[1] 12648
#Percent NA:
(sum(is.na(GabaGlu_Cq_AllSubjects_Scaled_QCed))/length(GabaGlu_Cq_AllSubjects_Scaled_QCed))*100
#[1] 0.1502214
#A very small percentage are NA now that we have removed the bad samples and bad genes.

#Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
GabaGlu_Cq_AllSubjects_Scaled_QCed[is.na(GabaGlu_Cq_AllSubjects_Scaled_QCed)]<-0

pdf("GabaGlu_Boxplot_CqZscore_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score")
dev.off()
#Huge differences in relative Cq across samples. Seems to partially correlate with card. 
#So we definitely need to do something to normalize for overall signal levels on the cards (which is why PCR normally normalizes by housekeeping genes)
#... but that presumably won't help the lack of correlation between replicate samples, since it only corrects for overall Cq differences across samples.

#For comparison with later -DeltaCq plot:
pdf("Boxplot_CqZscore_perSample_AllGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score", ylim=c(-10,10))
dev.off()

#Let's look at the "housekeeping genes" specifically:

row.names(GabaGlu_Cq_AllSubjects_Scaled_QCed)[83:93]
#[1] "HPRT1" "GUSB"  "ACTB"  "B2M"   "HMBS"  "IPO8"  "PGK1"  "RPLP0" "TBP"   "TFRC"  "UBC"  

pdf("Boxplot_CqZscore_perSample_HousekeepingGenes.pdf", width=20, height=5)
boxplot(data.frame(GabaGlu_Cq_AllSubjects_Scaled_QCed)[83:93,], cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for housekeeping genes", xlab="Sample ID", ylab="Cq Z-score")
dev.off()
#Yes, the housekeeping genes are tightly variable across samples. They definitely seem to vary with card, but also vary by sample in a manner that does not necessarily correlate across replicates.
#One of the replicate samples from subjects 1 and 37 look particulary extreme/variable. This was also true for sample 12 in an earlier analysis (before we discovered it should be thrown out for basic RNA quality reasons).
#Notably, subject #1 also had a lower replicate sample-sample correlation than the other. 
#This was also true in an earlier analysis for one of the replicate samples for subject #12 also had bizarrely high Cq values for non-housekeeping genes. 
#We may want to consider throwing these samples out, since everything will be normalized against these funky-looking housekeeping gene values.

GabaGlu_TrimmedMeanHousekeepingZScore<-apply(GabaGlu_Cq_AllSubjects_Scaled_QCed[83:93,], 2, function(y) mean(y, trim=0.2, na.rm=T))
hist(GabaGlu_TrimmedMeanHousekeepingZScore, breaks=20)
#Note: In a previous analysis, 12.1 was the sample with particularly high housekeeping gene measurements. Now everything is roughly between -2 and 2.

GabaGlu_MeanHousekeepingZScore<-apply(GabaGlu_Cq_AllSubjects_Scaled_QCed[83:93,], 2, function(y) mean(y, na.rm=T))
hist(GabaGlu_MeanHousekeepingZScore, breaks=20)
#Note: In a previous analysis, 12.1 had really high mean housekeeping gene measurements. This subject also had a really low RIN (<3), but the other sample for 12 didn't look as weird (???). Subject ID 50 also didn't look as weird despite having a low RIN.

GabaGlu_MeanHousekeepingZScore[GabaGlu_MeanHousekeepingZScore<(-2)]
# 37 
# -2.361561
#So the mean for one of the 37 replicates is also a little weird (not just its distribution), although not as extreme as 12 had been.

GabaGlu_TrimmedMeanNonHousekeepingZScore<-apply(GabaGlu_Cq_AllSubjects_Scaled_QCed[-c(83:93),], 2, function(y) mean(y, trim=0.2,  na.rm=T))
GabaGlu_MeanNonHousekeepingZScore<-apply(GabaGlu_Cq_AllSubjects_Scaled_QCed[-c(83:93),], 2, function(y) mean(y, na.rm=T))

plot(GabaGlu_TrimmedMeanNonHousekeepingZScore~GabaGlu_TrimmedMeanHousekeepingZScore)
cor(GabaGlu_TrimmedMeanNonHousekeepingZScore,GabaGlu_TrimmedMeanHousekeepingZScore)
#[1] 0.8922094


plot(GabaGlu_MeanNonHousekeepingZScore~GabaGlu_MeanHousekeepingZScore)

cor(GabaGlu_MeanNonHousekeepingZScore,GabaGlu_MeanHousekeepingZScore)
#[1] 0.8855489

#Housekeeping gene and non-housekeeping gene measurements are strongly correlated in all but a few samples.


#Trying PCA:

pca_GabaGlu<-prcomp(t(GabaGlu_Cq_AllSubjects_Scaled_QCed))
tmp<-pca_GabaGlu$x[,1:10]
dim(tmp)
#[1] 136  10
rownames(tmp)<-colnames(GabaGlu_Cq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_GabaGlu.csv")

tmp<-pca_GabaGlu$rotation[,1:10]
write.csv(tmp, "pca_GabaGlu_Eigenvectors.csv")

png("PCA_ScreePlot_GabaGlu.png")
plot(summary(pca_GabaGlu)$importance[2,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 is a huge amount of the variance (>60%)

png("PCA_ScreePlot2_GabaGlu.png")
plot(summary(pca_GabaGlu)$importance[3,]~(c(1:length(summary(pca_GabaGlu)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()

png("PC1vsPC2_GabaGlu_byDiagnosis.png")
plot(pca_GabaGlu$x[,1]~pca_GabaGlu$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForGabaGluCqMatrix_QCed$Diagnosis))
dev.off()
#Not diagnosis.
#Also not driven by a single major outlier
#In a previous analysis, sample 12.1 was pretty far out on PC1:
hist(pca_GabaGlu$x[,1])
#skewed distribution


pdf("PC1_GabaGlu_byTrimmedMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_TrimmedMeanHousekeepingZScore)
dev.off()
cor(pca_GabaGlu$x[,1],GabaGlu_TrimmedMeanHousekeepingZScore)
#[1] 0.9115217
#PC1 *tightly* correlates with mean houskeeping Zscore. So samples with greater Cq must correlate better with samples with greater Cq.
#This is a little surprising, because the PC scores are based on the correlation across samples (which shouldn't be affected by overall Cq, unless overall Cq is having some sort of non-linear effect on the signal for some of the genes)
#Therefore, it is possible that the measurements for a sample are qualitatively different when there is greater Cq.
#... so I'm not sure if normalizing by mean Housekeeping gene will help or not. 

pdf("PC1_GabaGlu_byMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_GabaGlu$x[,1]~GabaGlu_MeanHousekeepingZScore)
dev.off()
cor(pca_GabaGlu$x[,1],GabaGlu_MeanHousekeepingZScore)
#[1] 0.8987354


#What about some of the other numeric variables, esp. those related to RNA quality?

colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed)

cor(cbind(pca_GabaGlu$x[,1],GabaGlu_TrimmedMeanHousekeepingZScore, GabaGlu_MeanHousekeepingZScore, SubjectInfo_OrderedForGabaGluCqMatrix_QCed[,c(7, 9, 11, 16, 24:28, 38:42)]))

#Nope, It is by far most closely related to Housekeeping gene expression. The correlation with other RNA quality metrics (e.g., RNA concentration, RIN) is much weaker (R=-0.40 & -0.46, respectively). 
#It must represent at least partially something at the level of the PCR card instead of the extraction.

Temp<-cor(cbind(pca_GabaGlu$x[,c(1:4)],GabaGlu_TrimmedMeanHousekeepingZScore, GabaGlu_MeanHousekeepingZScore, SubjectInfo_OrderedForGabaGluCqMatrix_QCed[,c(7, 9, 11, 16, 24:28, 38:42)]))

write.csv(Temp, "GabaGlu_CorMatrix_PCA_vs_NumericVar_Cq.csv")


#**************************************************************
