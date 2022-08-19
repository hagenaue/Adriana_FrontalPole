

##### At this point I went back to the analysis of the DA5HT Cards:


#QC 4. How much do the replicate samples correlate?


#DA5HT:
#There are 48 measurements per sample on each card.  All cards should follow the same order of targets.
#Sanity Check:
temp<-Concatenated_DA5HT_wSubjectInfo[Concatenated_DA5HT_wSubjectInfo$ID=="1",]
cbind(temp$GeneSymbol[c(1:48)],temp$GeneSymbol[c(49:96)])
#Yep.
plot(temp$Cq[c(1:48)], temp$Cq[c(49:96)])
#Ooh - nice correlation.

#Percentage difference between the cards for each target:
(abs(temp$Cq[c(1:48)]-temp$Cq[c(49:96)])/apply(cbind(temp$Cq[c(1:48)],temp$Cq[c(49:96)]), 1, mean, na.rm=T))*100

#Let's loop it:

length(names(table(Concatenated_DA5HT_wSubjectInfo$ID)))
#[1] 72

#This code is going to be a little bit weird, because 8 subjects (subject numbers 1, 18, 22, 31, 41, 46, 58, 61) for the DA5HT cards actually have quadruplicate measurements (4 cards), with cards 57 and 58 entirely composed of samples that had been run on other cards.
#For the sake of simplicity, perhaps we could try just looking at the data from cards up through 56? to begin with, and then create an expanded version to look at the later cards for just the subjects that have them?

#Or add an if/else statement? And perhaps just look at how the data from cards 57 and 58 compare to each other?? Otherwise adding in the additional replicates to our pipeline code could get complicated fast.

cbind(temp$GeneSymbol[c(97:144)],temp$GeneSymbol[c(145:192)])
#Same order, correct indices.
plot(temp$Cq[c(97:144)], temp$Cq[c(145:192)])
#Ooh...pretty.

ReplicateCor_ForEachSubjectID_DA5HT<-matrix(0, 80, 1)
row.names(ReplicateCor_ForEachSubjectID_DA5HT)<-rep("X", 80)
row.names(ReplicateCor_ForEachSubjectID_DA5HT)[c(1:72)]<-names(table(Concatenated_DA5HT_wSubjectInfo$ID))
row.names(ReplicateCor_ForEachSubjectID_DA5HT)[c(73:80)]<-c("1_extra", "18_extra", "22_extra", "31_extra", "41_extra", "46_extra", "58_extra", "61_extra")
row.names(ReplicateCor_ForEachSubjectID_DA5HT)

ExtraIndices<-data.frame(SubjectID=c("1", "18", "22", "31", "41", "46", "58", "61"), ExtraIndices=c(73:80), stringsAsFactors = F)
str(ExtraIndices)

Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates<-rep(0, length(Concatenated_DA5HT_wSubjectInfo[,1]))


for (i in c(1:72)){
  
  print(names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i])
  
  temp<-Concatenated_DA5HT_wSubjectInfo[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i],] 
  
  if(dim(temp)[1]==96){
    print(dim(temp))
  }else{
    print(dim(temp))
  }
}


for (i in c(1:72)){
  
  print(names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i])
  
  temp<-Concatenated_DA5HT_wSubjectInfo[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i],] 
  
  if(dim(temp)[1]==96){
    
    ReplicateCor_ForEachSubjectID_DA5HT[i,1]<-cor(x=temp$Cq[c(1:48)], y=temp$Cq[c(49:96)], use="pairwise.complete.obs")
    
    PercentageDiff_BetweenReplicates<-(abs(temp$Cq[c(1:48)]-temp$Cq[c(49:96)])/apply(cbind(temp$Cq[c(1:48)],temp$Cq[c(49:96)]), 1, mean, na.rm=T))*100
    
    Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i]]<-rep(PercentageDiff_BetweenReplicates, 2)
    
    #For tracking loop progress:
    print(length(PercentageDiff_BetweenReplicates))
    print(length(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i]] ))
    
    rm(temp)
    rm(PercentageDiff_BetweenReplicates)
    
  }else{
    
    #Recording the correlation between the first two replicate samples:
    ReplicateCor_ForEachSubjectID_DA5HT[i,1]<-cor(x=temp$Cq[c(1:48)], y=temp$Cq[c(49:96)], use="pairwise.complete.obs")
    
    #Recording the correlation between the extra two replicate samples:
    CurrentExtraIndex<-ExtraIndices[ExtraIndices[,1]==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i],2] 
    
    ReplicateCor_ForEachSubjectID_DA5HT[CurrentExtraIndex,1]<-cor(x=temp$Cq[c(145:192)], y=temp$Cq[c(97:144)], use="pairwise.complete.obs")
    
    #Recording the percentage difference between the first two replicate samples:
    PercentageDiff_BetweenReplicates<-(abs(temp$Cq[c(1:48)]-temp$Cq[c(49:96)])/apply(cbind(temp$Cq[c(1:48)],temp$Cq[c(49:96)]), 1, mean, na.rm=T))*100
    
    Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i]][c(1:96)]<-rep(PercentageDiff_BetweenReplicates, 2) 
    
    rm(PercentageDiff_BetweenReplicates)
    
    #Recording the percentage difference between the extra two replicate samples:
    PercentageDiff_BetweenReplicates<-(abs(temp$Cq[c(97:144)]-temp$Cq[c(145:192)])/apply(cbind(temp$Cq[c(97:144)],temp$Cq[c(145:192)]), 1, mean, na.rm=T))*100
    
    Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i]][c(97:192)] <-rep(PercentageDiff_BetweenReplicates, 2) 
    
    #For tracking loop progress:
    print(length(PercentageDiff_BetweenReplicates))
    print(length(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates[Concatenated_DA5HT_wSubjectInfo$ID==names(table(Concatenated_DA5HT_wSubjectInfo$ID))[i]] ))
    
    rm(PercentageDiff_BetweenReplicates)
    rm(temp)
  }
}


pdf("DA5HT_Histogram_PercentageDiffBetweenReplicates.pdf", width=4, height=4)
hist(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates, xlab="Per Gene: % Diff Between Replicates", main="", breaks=100)
dev.off()
#Very pretty - really high degree of convergence. The vast majority of measurements agree within 5% of their replicate.

length(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates)
#[1] 7680

sum(is.na(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates)==F)
#[1] 7586

sum(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates<20, na.rm=T)
#[1] 7586

#Wow - 100% of the measurements match between replicates within 20%! :)


setwd("~/Documents/Microarray Gen/FrontalPole/Output/FinalOutput_AcrossDatasets")

pdf("BothDatasets_PercentDiffBetweenReplicates_vs_Cq.pdf", width=5, height=5)
plot(PercentageDiff_BetweenReplicates~Cq, data=Concatenated_GabaGlu_wSubjectInfo[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE,], cex=0.3, col=rgb(red=1,green=0,blue=0,alpha=0.1), ylab="% Difference Between Replicates", xlab="Cq") 
points(PercentageDiff_BetweenReplicates~Cq, data=Concatenated_DA5HT_wSubjectInfo, cex=0.3, col=rgb(red=0,green=0,blue=1,alpha=0.1))
dev.off()

pdf("BothDatasets_PercentDiffBetweenReplicates_vs_Cq_ZoomedIn.pdf", width=5, height=5)
plot(PercentageDiff_BetweenReplicates~Cq, data=Concatenated_GabaGlu_wSubjectInfo[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE,], cex=0.3, col=rgb(red=1,green=0,blue=0,alpha=0.1), ylab="% Difference Between Replicates", xlab="Cq", ylim=c(0,50)) 
points(PercentageDiff_BetweenReplicates~Cq, data=Concatenated_DA5HT_wSubjectInfo, cex=0.3, col=rgb(red=0,green=0,blue=1,alpha=0.1))
dev.off()

#hmm...how about a more digestible summary graphic:
head(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates)
head(Concatenated_DA5HT_wSubjectInfo$Cq)
trunc(Concatenated_DA5HT_wSubjectInfo$Cq)

DA5HT_CqVsMedianPercentDiffBetweenReplicates<-tapply(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates, trunc(Concatenated_DA5HT_wSubjectInfo$Cq), function(y) median(y, na.rm=TRUE))

DA5HT_CqVs95_PercentDiffBetweenReplicates<-tapply(Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates, trunc(Concatenated_DA5HT_wSubjectInfo$Cq), function(y) quantile(y, probs=0.95, na.rm=TRUE))

GabaGlu_CqVsMedianPercentDiffBetweenReplicates<-tapply(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE],trunc(Concatenated_GabaGlu_wSubjectInfo$Cq[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE]), function(y) median(y, na.rm=TRUE))

GabaGlu_CqVs95_PercentDiffBetweenReplicates<-tapply(Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE],trunc(Concatenated_GabaGlu_wSubjectInfo$Cq[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE]), function(y) quantile(y, probs=0.95, na.rm=TRUE))

str(DA5HT_CqVsMedianPercentDiffBetweenReplicates)
DA5HT_Cq_Centered<-(as.numeric(names(DA5HT_CqVsMedianPercentDiffBetweenReplicates))+0.5)
GabaGlu_Cq_Centered<-(as.numeric(names(GabaGlu_CqVsMedianPercentDiffBetweenReplicates))+0.5)

pdf("BothDatasets_PercentDiffBetweenReplicates_vs_Cq_Summary.pdf", width=5, height=5)
plot(GabaGlu_CqVsMedianPercentDiffBetweenReplicates~GabaGlu_Cq_Centered, lwd=2, col="red", ylab="% Difference Between Replicates", xlab="Cq",pch=20)
lines(GabaGlu_CqVsMedianPercentDiffBetweenReplicates~GabaGlu_Cq_Centered, lwd=2, col="red")
lines(GabaGlu_CqVs95_PercentDiffBetweenReplicates~GabaGlu_Cq_Centered, lwd=1, col="red")
points(DA5HT_CqVsMedianPercentDiffBetweenReplicates~DA5HT_Cq_Centered, lwd=2, col="blue", pch=20)
lines(DA5HT_CqVsMedianPercentDiffBetweenReplicates~DA5HT_Cq_Centered, lwd=2, col="blue")
lines(DA5HT_CqVs95_PercentDiffBetweenReplicates~DA5HT_Cq_Centered, lwd=1, col="blue")
dev.off()

#Removing the values based on the single measurements at 9 and 10:
pdf("BothDatasets_PercentDiffBetweenReplicates_vs_Cq_Summary_2.pdf", width=5, height=5)
plot(GabaGlu_CqVsMedianPercentDiffBetweenReplicates[-c(4:5)]~GabaGlu_Cq_Centered[-c(4:5)], lwd=2, col="red", ylab="% Difference Between Replicates", xlab="Cq", pch=20, ylim=c(0,35))
lines(GabaGlu_CqVsMedianPercentDiffBetweenReplicates[-c(4:5)]~GabaGlu_Cq_Centered[-c(4:5)], lwd=2, col="red")
lines(GabaGlu_CqVs95_PercentDiffBetweenReplicates[-c(4:5)]~GabaGlu_Cq_Centered[-c(4:5)], lwd=1, col="red")
points(DA5HT_CqVsMedianPercentDiffBetweenReplicates~DA5HT_Cq_Centered, lwd=2, col="blue", pch=20)
lines(DA5HT_CqVsMedianPercentDiffBetweenReplicates~DA5HT_Cq_Centered, lwd=2, col="blue")
lines(DA5HT_CqVs95_PercentDiffBetweenReplicates~DA5HT_Cq_Centered, lwd=1, col="blue")
dev.off()

#Ooh... what if we combined them so that the values at the high end were a little more stable...

TempGabaGlu_PercentageDiff_BetweenReplicates<-Concatenated_GabaGlu_wSubjectInfo$PercentageDiff_BetweenReplicates[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE]

TempGabaGlu_Cq<-trunc(Concatenated_GabaGlu_wSubjectInfo$Cq[((Concatenated_GabaGlu_wSubjectInfo$ID==27 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 14.eds")|(Concatenated_GabaGlu_wSubjectInfo$ID==38 & Concatenated_GabaGlu_wSubjectInfo$Card=="Card 20.eds"))==FALSE])

BothDatasets_CqVsMedianPercentDiffBetweenReplicates<-tapply(c(TempGabaGlu_PercentageDiff_BetweenReplicates,Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates), c(TempGabaGlu_Cq, trunc(Concatenated_DA5HT_wSubjectInfo$Cq)), function(y) median(y, na.rm=TRUE))

BothDatasets_CqVs95_PercentDiffBetweenReplicates<-tapply(c(TempGabaGlu_PercentageDiff_BetweenReplicates,Concatenated_DA5HT_wSubjectInfo$PercentageDiff_BetweenReplicates), c(TempGabaGlu_Cq, trunc(Concatenated_DA5HT_wSubjectInfo$Cq)), function(y) quantile(y, probs=0.95, na.rm=TRUE))

#The values at 9 & 10 are really high probably because there aren't that many measurements from then. Let's see:

table(c(TempGabaGlu_Cq, trunc(Concatenated_DA5HT_wSubjectInfo$Cq)))
# 6    7    8    9   10   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36 
# 54   62   24    1    1    8  201  531  681 1592 2701 2730 2640 2628 1890 1444  906  535  472  552  363  380  299  209  144   56 
# 37   38 
# 12   10 

#Yeah... a median for 1 measurement is not particularly informative. Let's remove that.

#That might be useful information to include with the plot, along with the number of NAs in each Cq bin.



BothDatasets_Cqcentered<-as.numeric(names(BothDatasets_CqVsMedianPercentDiffBetweenReplicates))+0.5

pdf("BothDatasets_PercentDiffBetweenReplicates_vs_Cq_Summary_3.pdf", width=5, height=5)
plot(BothDatasets_CqVsMedianPercentDiffBetweenReplicates[-c(4:5)]~BothDatasets_Cqcentered[-c(4:5)], lwd=2, col="black", ylab="% Difference Between Replicates", xlab="Cq", pch=20, ylim=c(0,35))
lines(BothDatasets_CqVsMedianPercentDiffBetweenReplicates[-c(4:5)]~BothDatasets_Cqcentered[-c(4:5)], lwd=2, col="black")
lines(BothDatasets_CqVs95_PercentDiffBetweenReplicates[-c(4:5)]~BothDatasets_Cqcentered[-c(4:5)], lwd=1, col="black")
dev.off()


write.csv(ReplicateCor_ForEachSubjectID_DA5HT, "ReplicateCor_ForEachSubjectID_DA5HT.csv")

pdf("DA5HT_Boxplot_CorrBetweenReplicates.pdf", width=3, height=4)
boxplot(ReplicateCor_ForEachSubjectID_DA5HT[,1], ylab="Correlation betweeen replicate samples (R)")
dev.off()

summary((ReplicateCor_ForEachSubjectID_DA5HT[,1]))

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9844  0.9940  0.9957  0.9951  0.9972  0.9988 

# #Extremely tight. Whiskers & IQ all fall above 0.99. There are 3 "outlier" datapoints outside of that range. Let's see who they are:
ReplicateCor_ForEachSubjectID_DA5HT[ReplicateCor_ForEachSubjectID_DA5HT[,1]<0.99,]
# 21        34        50 
# 0.9844096 0.9860820 0.9881312
#Notably, these are not any of the samples that Adriana was concerned about, nor is their sample-sample correlation low enough to really arouse suspicioun.

#I'd like to quickly snoop at their correlation plots to see if those results are driven by outlier measurements:
temp<-Concatenated_DA5HT_wSubjectInfo[Concatenated_DA5HT_wSubjectInfo$ID=="21",]
plot(temp$Cq[c(1:48)], temp$Cq[c(49:96)])
#There is 1 measurement (with high Cq) that is pretty far off the trend line

temp<-Concatenated_DA5HT_wSubjectInfo[Concatenated_DA5HT_wSubjectInfo$ID=="34",]
plot(temp$Cq[c(1:48)], temp$Cq[c(49:96)])
#There is 1 measurement that is pretty far off the trend line

temp<-Concatenated_DA5HT_wSubjectInfo[Concatenated_DA5HT_wSubjectInfo$ID=="50",]
plot(temp$Cq[c(1:48)], temp$Cq[c(49:96)])
#There is 1 measurement that is pretty far off the trend line

#Based on these plots and the lack of replicate measurements with greater than 20% percentage difference, I'm inclined to just use the data as is without further cleaning.

Concatenated_DA5HT_wSubjectInfo$Cq_NoBad<-Concatenated_DA5HT_wSubjectInfo$Cq

#Looks good.

write.csv(Concatenated_DA5HT_wSubjectInfo, "Concatenated_DA5HT_wSubjectInfo.csv")


#It would be helpful to interpret this data in the context of the correlations across subjects (vs. within subjects). Our data frame is currently not well set up for that 

#If I could rearrange it, I could reuse a lot of my code from microarray/RNA-Seq analyses.

dim(Concatenated_DA5HT_wSubjectInfo)
#[1] 7680   73

7680/48
#[1] 160
160/2
#[1] 80

Concatenated_DA5HT_wSubjectInfo$ID[c(1:48)]
Concatenated_DA5HT_wSubjectInfo$ID[c(49:(48+48))]

table(Concatenated_DA5HT_wSubjectInfo$ID[seq(1, 7680, 48)])
#looks good: two samples per subject except for the 8 subjects that have 4 samples.

# 1 10 12 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31 32 33 34 37 38 39  4 40 41 42 43 44 45 46 47 48 49  5 50 51 
# 4  2  2  2  2  2  2  2  4  2  2  2  2  4  2  2  2  2  2  2  2  2  2  4  2  2  2  2  2  2  2  2  4  2  2  2  2  4  2  2  2  2  2  2 
# 52 53 54 55 56 57 58 59  6 60 61 62 63 64 65 66 67 68 69  7 70 71 72 73 74 75  8  9 
# 2  2  2  2  2  2  4  2  2  2  4  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 

DA5HT_Cq_AllSubjects<-matrix(0,48,160)
colnames(DA5HT_Cq_AllSubjects)<-Concatenated_DA5HT_wSubjectInfo$ID[seq(1, 7680, 48)]
colnames(DA5HT_Cq_AllSubjects)
#If we use the subject IDs, the columnames sometimes repeat - I thought R would be unhappy about it, but I don't see a .2 added to the names. Interesting.

FirstGeneForSample_iterative<-seq(1, 7680, 48)
LastGeneForSample_iterative<-seq(48, 7680, 48)

i<-1
#Sanity check:
Concatenated_DA5HT_wSubjectInfo$ID[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]
row.names(DA5HT_Cq_AllSubjects)<-Concatenated_DA5HT_wSubjectInfo$GeneSymbol[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]

for(i in c(1:160)){
  DA5HT_Cq_AllSubjects[,i]<-Concatenated_DA5HT_wSubjectInfo$Cq_NoBad[FirstGeneForSample_iterative[i]:LastGeneForSample_iterative[i]]
}

head(DA5HT_Cq_AllSubjects)
# 1      2      3      4      5      6      7      8      1      2      3      4      5      6      7      8      9
# ADRB1 24.245 23.982 24.608 24.017 23.872 24.560 24.412 23.761 23.581 24.125 24.656 24.517 24.143 24.746 24.971 24.227 24.317
# ADRB2 25.087 25.066 24.983 25.385 24.772 25.816 25.637 24.968 24.934 25.044 24.682 25.126 24.667 25.531 25.443 24.944 24.799
# COMT  24.579 24.757 24.539 24.838 24.480 24.476 25.011 24.249 23.588 24.507 25.110 24.355 24.159 24.332 24.862 24.633 24.230
# DBH   27.823 28.153 29.204 28.992 27.872 27.979 28.558 28.265 28.040 28.739 29.888 29.194 28.477 28.355 29.177 28.859 28.903
# DDC   32.773 33.446 32.634 32.706 33.048 32.914 34.089 34.110 32.681 33.205 32.143 32.642 32.937 32.687 33.408 35.332 33.250
# DRD1  25.515 24.959 25.285 25.047 24.952 25.217 25.165 25.232 25.167 24.983 25.248 25.053 24.928 25.226 25.251 25.309 25.150

colnames(Concatenated_DA5HT_wSubjectInfo)

SubjectInfo_OrderedForDA5HTCqMatrix<-Concatenated_DA5HT_wSubjectInfo[FirstGeneForSample_iterative, c(2, 24:69)]
head(SubjectInfo_OrderedForDA5HTCqMatrix)
str(SubjectInfo_OrderedForDA5HTCqMatrix)

# 'data.frame':	160 obs. of  47 variables:
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

#Checking to see whether dissection group and card are strongly related again.
table(SubjectInfo_OrderedForDA5HTCqMatrix$Card, SubjectInfo_OrderedForDA5HTCqMatrix$Dissecton.Group)
#Yes - each card contains samples from 1-3 dissection groups (with the exception of cards 57 & 58)


#As far as I can tell, no outlier samples need to be removed based on basic QC, but the samples with super low RIN and bad purity should be taken out to make basic standards:
colnames(SubjectInfo_OrderedForDA5HTCqMatrix)

dim(DA5HT_Cq_AllSubjects)
#[1]  48 160
DA5HT_Cq_AllSubjects_QCed<-DA5HT_Cq_AllSubjects[,SubjectInfo_OrderedForDA5HTCqMatrix$LowQualityRNA==F]
dim(DA5HT_Cq_AllSubjects_QCed)
#[1]  48 154

dim(SubjectInfo_OrderedForDA5HTCqMatrix)
#[1] 160  47
SubjectInfo_OrderedForDA5HTCqMatrix_QCed<-SubjectInfo_OrderedForDA5HTCqMatrix[SubjectInfo_OrderedForDA5HTCqMatrix$LowQualityRNA==F,]
dim(SubjectInfo_OrderedForDA5HTCqMatrix_QCed)
#[1] 154  47

#Sanity check:
colnames(DA5HT_Cq_AllSubjects_QCed)
SubjectInfo_OrderedForDA5HTCqMatrix_QCed$ID
#Yes - sample order matches


#What about bad genes?

dim(DA5HT_Cq_AllSubjects_QCed)
#[1]  48 154

DA5HT_NumberNA_PerGene<-apply(DA5HT_Cq_AllSubjects_QCed, 1, function(y) sum(is.na(y)))
hist(DA5HT_NumberNA_PerGene)

#This is useful information - let's output it, along with average Cq:

DA5HT_AverageCq_PerGene<-apply(DA5HT_Cq_AllSubjects_QCed, 1, function(y) mean(y, na.rm=T))
pdf("DA5HT_Scatterplot_AverageCq_vs_NumberNAPerGene.pdf", width=5, heigh=5)
plot(DA5HT_AverageCq_PerGene~DA5HT_NumberNA_PerGene, ylab="Average Cq per probe", xlab= "# of bad measurements per probe")
dev.off()
#The two genes with a large number of NAs also have really high Cq (close to 35)

write.csv(cbind(DA5HT_AverageCq_PerGene, DA5HT_NumberNA_PerGene), "DA5HT_AverageCq_NumberNA_PerGene.csv")


#Which genes are they?
apply(DA5HT_Cq_AllSubjects_QCed, 1, function(y) sum(is.na(y)))

# ADRB1   ADRB2    COMT     DBH     DDC    DRD1    DRD2    DRD3    DRD4    DRD5   GAPDH   HTR1A   HTR1B   HTR1D   HTR1E   HTR1F   HTR2A 
# 0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0 
# HTR2B   HTR2C   HTR3A   HTR3B    HTR4   HTR5A    HTR6    HTR7    MAOA    MAOB SLC18A1 SLC18A2  SLC6A3  SLC6A4      TH   HPRT1    GUSB 
# 0       0       0       0       0       0       0       0       0       0      27       0      20       0       0       0       0 
# B2M    IPO8    TFRC   YWHAZ    PPIA  POLR2A   CASC3 GADD45A    PUM1   PSMC4    EMC7     GPI   RAB7A   REEP5 
# 0       0       0       0       0       0       0       0       0       0       0       0       0       0 

#Yes!  All of the bad measurements come from two genes: SLC18A1 and SLC6A3
#SLC18A1:
27/154
#[1] 0.1753247
#Over 17% of the measurements are bad for SLC18A1
20/154
#[1] 0.1298701
# 13% of the measurements are bad for SLC6A3
#Hmmm... those aren't quite as extreme as for the DA5HT cards. Worth throwing out? Probably. For consistency sake, let's use the same filter as last time:


DA5HT_Cq_AllSubjects_QCed2<-DA5HT_Cq_AllSubjects_QCed[DA5HT_NumberNA_PerGene<20,]
dim(DA5HT_Cq_AllSubjects_QCed2)
#[1]  46 154

#Determining which indices are the housekeeping genes now:
row.names(DA5HT_Cq_AllSubjects_QCed)[c(11,33:48)]
# [1] "GAPDH"   "HPRT1"   "GUSB"    "B2M"     "IPO8"    "TFRC"    "YWHAZ"   "PPIA"    "POLR2A"  "CASC3"   "GADD45A" "PUM1"    "PSMC4"  
# [14] "EMC7"    "GPI"     "RAB7A"   "REEP5" 

row.names(DA5HT_Cq_AllSubjects_QCed2)[c(11,31:46)]
# [1] "GAPDH"   "HPRT1"   "GUSB"    "B2M"     "IPO8"    "TFRC"    "YWHAZ"   "PPIA"    "POLR2A"  "CASC3"   "GADD45A" "PUM1"    "PSMC4"  
# [14] "EMC7"    "GPI"     "RAB7A"   "REEP5" 

#At this point, no samples or genes have NAs.


#######################

#Looking at sample-sample correlations

DA5HT_SubjectSubjectCorMatrix<-cor(DA5HT_Cq_AllSubjects_QCed2, use="pairwise.complete.obs")

pdf("DA5HT_SubjectSubjectCorMatrix.pdf", height=14, width=14)
heatmap(DA5HT_SubjectSubjectCorMatrix)
dev.off()
#Interesting - many of the replicates most closely resemble each other, but certainly not all (maybe 50%?) - more than for the DA5HT card. The correlation blocks are also less obvious than in the DA5HT dataset, although still visible.


plot(DA5HT_Cq_AllSubjects_QCed2[,1]~DA5HT_Cq_AllSubjects_QCed2[,2])
plot(DA5HT_Cq_AllSubjects_QCed2[,4]~DA5HT_Cq_AllSubjects_QCed2[,5])
#Unlike the DA5HT dataset, there aren't any super high leverage points, but to keep the analyses consistent with the other dataset I will try redoing the correlation analysis after z-scoring the data first.

DA5HT_Cq_AllSubjects_Scaled_QCed<-t(scale(t(DA5HT_Cq_AllSubjects_QCed2), center=T, scale=T))
head(DA5HT_Cq_AllSubjects_Scaled_QCed)
#Sanity check:
apply(DA5HT_Cq_AllSubjects_Scaled_QCed, 1, function(y) mean(y, na.rm=T))
apply(DA5HT_Cq_AllSubjects_Scaled_QCed, 1, function(y) sd(y, na.rm=T))
#Looks good.

DA5HT_SubjectSubjectCorMatrix_Scaled_QCed<-cor(DA5HT_Cq_AllSubjects_Scaled_QCed, use="pairwise.complete.obs")
DA5HT_SubjectSubjectCorMatrix_Scaled_QCed

pdf("DA5HT_SubjectSubjectCorMatrix_Scaled_QCed.pdf", height=14, width=14)
heatmap(DA5HT_SubjectSubjectCorMatrix_Scaled_QCed)
dev.off()
#Both the correlation between replicates and the correlation blocks are more obvious now. I would say that at this point the vast majority of the replicates cluster together - this is a very pretty dataset.


#Dealing with the NA values and looking at overall variation in Cq values across samples:

#Just to confirm that there aren't any NA values left:
sum(is.na(DA5HT_Cq_AllSubjects_Scaled_QCed))
#[1] 0

#This code is no longer necessary:
# #Since this is scale data, the mean value is 0 for all genes. Let's replace NAs with that:
# DA5HT_Cq_AllSubjects_Scaled_QCed[is.na(DA5HT_Cq_AllSubjects_Scaled_QCed)]<-0

pdf("DA5HT_Boxplot_CqZscore_perSample_AllGenes.pdf", width=20, height=5)
boxplot(data.frame(DA5HT_Cq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score")
dev.off()
#Huge differences in relative Cq across samples. Seems to partially correlate with card. 
#So we definitely need to do something to normalize for overall signal levels on the cards (which is why PCR normally normalizes by housekeeping genes)
#In a previous analysis, both replicate samples for subjects 43 and 12 had high Cq, but they were later removed during basic QC once we had the full RNA metrics available because they had bad RNA metrics)

#For comparison with later -DeltaCq plot:
pdf("DA5HT_Boxplot_CqZscore_perSample_AllGenes_ylim10.pdf", width=20, height=5)
boxplot(data.frame(DA5HT_Cq_AllSubjects_Scaled_QCed), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for all genes", xlab="Sample ID", ylab="Cq Z-score", ylim=c(-10,10))
dev.off()

#Let's look at the "housekeeping genes" specifically:

row.names(DA5HT_Cq_AllSubjects_Scaled_QCed)[c(11,31:46)]
# [1] "GAPDH"   "HPRT1"   "GUSB"    "B2M"     "IPO8"    "TFRC"    "YWHAZ"   "PPIA"    "POLR2A"  "CASC3"   "GADD45A" "PUM1"    "PSMC4"  
# [14] "EMC7"    "GPI"     "RAB7A"   "REEP5"  

pdf("DA5HT_Boxplot_CqZscore_perSample_HousekeepingGenes.pdf", width=20, height=5)
boxplot(data.frame(DA5HT_Cq_AllSubjects_Scaled_QCed)[c(11,31:46),], cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of Cq Z-scores per sample for housekeeping genes", xlab="Sample ID", ylab="Cq Z-score")
dev.off()
#Yes, the housekeeping genes are tightly variable across samples. Unlike the GabaGlu cards, the housekeeping genes in the DA5HT cards seem to often show variability that is similar in the two replicate samples per subject (i.e., more likely to reflect something biological and/or technical at the level of dissection/extraction)


DA5HT_TrimmedMeanHousekeepingZScore<-apply(DA5HT_Cq_AllSubjects_Scaled_QCed[c(11,31:46),], 2, function(y) mean(y, trim=0.2, na.rm=T))
hist(DA5HT_TrimmedMeanHousekeepingZScore, breaks=20)
#The outliers this time are not as extreme (>2 average z-score for housekeeping genes)
DA5HT_TrimmedMeanHousekeepingZScore[DA5HT_TrimmedMeanHousekeepingZScore>2]
# 39       39       67       67       69       75       69       75 
# 2.370568 2.373525 2.653692 2.521647 2.195683 2.105624 2.161373 2.079592 
#Note - in a previous analysis, subjects 12 and 43 also fell in this group  (and were much more extreme), but were removed during basic QC due to bad quality RNA metrics for this analysis.

DA5HT_MeanHousekeepingZScore<-apply(DA5HT_Cq_AllSubjects_Scaled_QCed[c(11,31:46),], 2, function(y) mean(y, na.rm=T))
hist(DA5HT_MeanHousekeepingZScore, breaks=20)
#The outliers this time are not as extreme (>2 average z-score for housekeeping genes)
DA5HT_MeanHousekeepingZScore[DA5HT_MeanHousekeepingZScore>2]
# 39       39       67       67       69       69 
# 2.220699 2.210241 2.637728 2.551943 2.099785 2.114460 

DA5HT_TrimmedMeanNonHousekeepingZScore<-apply(DA5HT_Cq_AllSubjects_Scaled_QCed[-c(11,31:46),], 2, function(y) mean(y, trim=0.2,  na.rm=T))
DA5HT_MeanNonHousekeepingZScore<-apply(DA5HT_Cq_AllSubjects_Scaled_QCed[-c(11,31:46),], 2, function(y) mean(y, na.rm=T))

plot(DA5HT_TrimmedMeanNonHousekeepingZScore~DA5HT_TrimmedMeanHousekeepingZScore)
cor(DA5HT_TrimmedMeanNonHousekeepingZScore,DA5HT_TrimmedMeanHousekeepingZScore)
#[1] 0.949618
#Nice.

plot(DA5HT_MeanNonHousekeepingZScore~DA5HT_MeanHousekeepingZScore)
cor(DA5HT_MeanNonHousekeepingZScore,DA5HT_MeanHousekeepingZScore)
#[1] 0.936664

#Housekeeping gene and non-housekeeping gene measurements are strongly correlated.


#Trying PCA:

pca_DA5HT<-prcomp(t(DA5HT_Cq_AllSubjects_Scaled_QCed))
tmp<-pca_DA5HT$x[,1:10]
dim(tmp)
#[1] 154  10
rownames(tmp)<-colnames(DA5HT_Cq_AllSubjects_Scaled_QCed)
write.csv(tmp, "PCA_DA5HT.csv")

tmp<-pca_DA5HT$rotation[,1:10]
write.csv(tmp, "pca_DA5HT_Eigenvectors.csv")

png("PCA_ScreePlot_DA5HT.png")
plot(summary(pca_DA5HT)$importance[2,]~(c(1:length(summary(pca_DA5HT)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Proportion of Variance Explained", col=2)
dev.off()
#PC1 is a huge amount of the variance (>50%)

png("PCA_ScreePlot2_DA5HT.png")
plot(summary(pca_DA5HT)$importance[3,]~(c(1:length(summary(pca_DA5HT)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC#", ylab="Cumulative Proportion of Variance Explained", col=2)
dev.off()

png("PC1vsPC2_DA5HT_byDiagnosis.png")
plot(pca_DA5HT$x[,1]~pca_DA5HT$x[,2], xlab="PC2", ylab="PC1", main="PC1 vs. PC2", col=as.factor(SubjectInfo_OrderedForDA5HTCqMatrix_QCed$Diagnosis))
dev.off()
#PC2 seems likely to correlate with diagnosis.
# not driven by a single major outlier.
#Note, in a previous analysis, samples from 43 and 12 were pretty far out on PC1, but after obtaining RNA quality metrics were removed during basic QC.
hist(pca_DA5HT$x[,1], breaks=20)
pca_DA5HT$x[pca_DA5HT$x[,1]<(-10),1]
# 39        39        67        67        69        75        69        75 
# -13.16618 -12.41963 -15.26264 -14.00059 -13.22647 -11.93320 -12.84505 -12.03870
#It's the same samples that had the highest HK Cq

pdf("PC1_DA5HT_byTrimmedMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~DA5HT_TrimmedMeanHousekeepingZScore)
dev.off()
cor(pca_DA5HT$x[,1],DA5HT_TrimmedMeanHousekeepingZScore)
#[1] -0.9872029
#Yep, definitely looks like it is just tracking an overall variability in Cq 
# although note that this dataset contains a large % of HK - 17/46 remaining genes, or ~1/3 HK - so it makes sense that PC1 would closely track HK expression

pdf("PC1_DA5HT_byMeanHousekeepingZScore.pdf", width=5, height=5)
plot(pca_DA5HT$x[,1]~DA5HT_MeanHousekeepingZScore)
dev.off()
cor(pca_DA5HT$x[,1],DA5HT_MeanHousekeepingZScore)
#-0.9836505

#It doesn't even seem to make sense to examine the correlation with other variables, since almost all variation is predicted by HK

#But I probably should examine the relationship between PC2 and PC3 and other variables to determine what should be included in a model examining the influence of diagnosis on HK expression.


#What about some of the other numeric variables, esp. those related to RNA quality?

colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed)

cor(cbind(pca_DA5HT$x[,c(1:4)], DA5HT_TrimmedMeanHousekeepingZScore, DA5HT_MeanHousekeepingZScore, SubjectInfo_OrderedForDA5HTCqMatrix_QCed[,c(7, 9, 11, 16, 24:28, 38:42)]))

Temp<-cor(cbind(pca_DA5HT$x[,c(1:4)], DA5HT_TrimmedMeanHousekeepingZScore, DA5HT_MeanHousekeepingZScore, SubjectInfo_OrderedForDA5HTCqMatrix_QCed[,c(7, 9, 11, 16, 24:28, 38:42)]))

write.csv(Temp, "DA5HT_CorMatrix_PCA_vs_NumericVar_Cq.csv")

#PC1 is HK, RNAConc, RIN, 260/280 (i.e. overall Cq)
#PC2 is Hours final, pH, and RIN, 28s/18s (i.e. degradation)
#PC3/PC4 I'm not sure...


#**************************************************************