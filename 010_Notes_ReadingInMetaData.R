

#Adriana's Frontal Pole PCR Data
## Goal: Figure out why the "housekeeping genes" are showing large amounts of variability
## Megan Hagenauer
## Initial analysis: August 3, 2020
#Re-run 8/25 and 9/1 with data updates
#Re-run 9/21 to fix minor bug (dissection group coded as factor vs. numeric)
#Re-run again 10/13 after receiving additional potential co-variates (RIN, RNAconc)
#Re-run again because ID was being treated as numeric in the DA5HT dataset (but not the GabaGlu dataset???)... but then realized it didn't effect the results.
#Rerun again in early november because new RNA purity/integrity information became available (TTP & TZP) for the samples and I wanted to take a peek at how RNA measures (concentration, RIN, purity) differed in the data before being cleaned up (TTP vs. TZP)
#And re-run one more time (hopefully!) on 11/25 to fix the analysis to exclude subject 43 (bad concentration & purity 260/230) for the entire analysis and not just the final outlier-removed analyses & to fix a bug (GAPDH in DA5HT dataset coded as a target gene), and add annotation.

## Most recent workspace saved as...


####################

#General notes:

##This was re-run 2020 08 25 because I discovered that we had been missing data from card 22.eds (sample IDs 41-44) on the server. Adriana reoutputted it from the proprietary website and uploaded it, so it is now included.
##It was re-run again 2020 09 01 to add a few more early QC steps and at that point I also cleaned up the code a little bit.
## And then Re-run 9/21 one more time to fix minor bug (dissection group coded as factor vs. numeric). At that point I also ran some additional analyses to double-check proper model choice (in particular, to determine whether the addition of Card as a technical co-variate might be overfitting the data) and to double-check that dissection wasn't as useful as a co-variate. The results from this analysis were approved over e-mail and formally written up (FrontalPole_AnalysisMethods_Formal_20200917).
##... but after discussing the write-up, we realized that some minor changes might be useful (e.g., evaluation of Evan's RIN and RNAconc values as potential co-variates), so the analysis was run again 10/13/2020.  In the process, I found what I believe is a decimal place error in Evan's rRNA Ratio values (28s/18s). I fixed this, but then did not re-run all code again because it seemed redundant with RIN.
#Re-run again because ID was being treated as numeric in the DA5HT dataset (but not the GabaGlu dataset???)... but then realized it didn't effect the results.
#Rerun again in early november because new RNA purity/integrity information became available (TTP & TZP) for the samples and I wanted to take a peek at how RNA measures (concentration, RIN, purity) differed in the data before being cleaned up (TTP vs. TZP)
#And re-run one more time (hopefully!) on 11/25 to fix the analysis to exclude subject 43 (bad concentration & purity 260/230) for the entire analysis and not just the final outlier-removed analyses & to fix a bug (GAPDH in DA5HT dataset coded as a target gene)

##################################

#I ran this code after running all analyses to pull up the information for all of the packages:

print(sessionInfo())

# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] multtest_2.32.0     Biobase_2.36.2      BiocGenerics_0.22.0 plyr_1.8.4          lme4_1.1-13        
# [6] Matrix_1.2-10       nlme_3.1-131       
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.4.6       compiler_3.4.1     pillar_1.3.1       nloptr_1.0.4       tools_3.4.1       
# [6] tibble_2.1.1       lattice_0.20-35    mgcv_1.8-17        pkgconfig_2.0.2    rlang_0.3.1       
# [11] rstudioapi_0.6     SparseM_1.77       dplyr_0.8.0.1      MatrixModels_0.4-1 stats4_3.4.1      
# [16] grid_3.4.1         nnet_7.3-12        tidyselect_0.2.5   glue_1.3.1         R6_2.2.2          
# [21] survival_2.41-3    minqa_1.2.4        purrr_0.3.2        car_2.1-5          magrittr_1.5      
# [26] MASS_7.3-47        splines_3.4.1      assertthat_0.2.1   pbkrtest_0.4-7     quantreg_5.33     
# [31] crayon_1.3.4   

#Most recent version:
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.4.6       rstudioapi_0.6     magrittr_1.5       splines_3.4.1      MASS_7.3-47        tidyselect_0.2.5  
# [7] lattice_0.20-35    R6_2.2.2           rlang_0.3.1        minqa_1.2.4        car_2.1-5          dplyr_0.8.0.1     
# [13] tools_3.4.1        nnet_7.3-12        parallel_3.4.1     pbkrtest_0.4-7     grid_3.4.1         mgcv_1.8-17       
# [19] quantreg_5.33      MatrixModels_0.4-1 assertthat_0.2.1   tibble_2.1.1       crayon_1.3.4       purrr_0.3.2       
# [25] nloptr_1.0.4       glue_1.3.1         pillar_1.3.1       compiler_3.4.1     SparseM_1.77       pkgconfig_2.0.2   


#*************************************

#Reading in Files:

#Reading in general experiment info:

setwd("~/Documents/Microarray Gen/FrontalPole/Sample_MetaData")

#Updated 11/25/20 to include GAPDH:
HousekeepingGenes<-read.csv("HouseKeepingGenes_UsedInPCR_FrontalPole.csv", header=T, stringsAsFactors = F)
str(HousekeepingGenes)
# 'data.frame':	25 obs. of  1 variable:
#   $ GeneSymbol: chr  "18S" "ACTB" "B2M" "CASC3" ...

#Earlier versions:
#SubjectInfo<-read.csv("Balancing RNA Extraction04162018_UseThis_Simplified.csv", header=T, stringsAsFactors = F)
#SubjectInfo<-read.csv("Balancing RNA Extraction04162018_UseThis_Simplified_wRINconc.csv", header=T, stringsAsFactors = F)

SubjectInfo<-read.csv("Balancing RNA Extraction04162018_UseThis_Simplified_wAllTTP_TZP.csv", header=T, stringsAsFactors = F)

str(SubjectInfo)
# 'data.frame':	72 obs. of  41 variables:
#   $ Barcode                             : chr  "B007375A" "B010169A" "B000014A" "B010842A" ...
# $ Subject.Number                      : int  5000 5313 2169 4087 5066 2950 4383 4619 4819 3711 ...
# $ Cohort                              : chr  "Cohort 11" "Cohort 13" "Dep Cohort 1" "Dep Cohort 6" ...
# $ Diagnosis                           : chr  "Control" "Schiz" "Control" "BP" ...
# $ Age                                 : int  77 42 18 51 47 45 40 34 79 59 ...
# $ Gender                              : chr  "M" "M" "M" "F" ...
# $ pH                                  : num  6.79 6.77 6.97 6.81 6.6 7.05 6.77 6.7 6.75 6.55 ...
# $ AFS                                 : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Hours.Final                         : num  16.5 23.7 22 16.8 21 20 26.3 17.2 21 27 ...
# $ Slab.Format                         : chr  "SLAB" "SLAB" "SLAB" "SLAB" ...
# $ Slab.Number                         : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Dissecton.Group                     : int  1 1 1 1 1 1 2 2 2 2 ...
# $ ID                                  : int  1 2 3 4 5 6 7 8 9 10 ...
# $ TTP.ID                              : chr  "TTP010001" "TTP010002" "TTP010003" "TTP010004" ...
# $ Block.Weight..g.                    : num  0.72 0.66 0.76 0.56 0.7 0.89 0.85 1.02 0.66 0.87 ...
# $ TZP.ID                              : chr  "TZP010001" "TZP010002" "TZP010003" "TZP010004" ...
# $ TZP_MBNI_Nanodrop_RNAConc..ng.uL.   : num  436 343 394 417 399 ...
# $ TZP_MBNI_Nanodrop_260.280           : num  1.96 2.01 1.99 1.97 2 1.99 2 2 2.02 2 ...
# $ TZP_MBNI_Nanodrop_260.230           : num  1.93 2.18 2.02 2.1 2.19 2.23 2.3 2.3 2.25 2.28 ...
# $ TZP_DSC_Nanodrop_RNAConc..ng.uL.    : num  562 354 451 449 435 ...
# $ TZP_DSC_Nanodrop_260.280            : num  2.14 1.92 2.03 2.03 2.05 2.05 2.05 2.06 2.05 2.06 ...
# $ TZP_DSC_Nanodrop_260.230            : num  1.9 2.06 2.01 2.13 2.22 2.27 2.29 2.28 2.22 2.3 ...
# $ TZP_BioAnalyzer_RIN                 : num  7.7 8.1 7.9 8 8.2 8.1 6.6 7.7 8.6 8 ...
# $ TZP_BioAnalyzer_rRNA.Ratio..28s.18s.: num  2.4 2.3 3.1 2 2.4 2.7 14.4 2.4 2.3 2.3 ...
# $ TZP_Average.RNAConc..ng.uL.         : num  499 349 422 433 417 ...
# $ TZP_Average.260.280                 : num  2.05 1.97 2.01 2 2.03 2.02 2.03 2.03 2.04 2.03 ...
# $ TZP_Average.260.230                 : num  1.92 2.12 2.02 2.12 2.21 2.25 2.3 2.29 2.24 2.29 ...
# $ TTP.ID.1                            : chr  "TTP010001" "TTP010002" "TTP010003" "TTP010004" ...
# $ Initial.Volume..uL.                 : int  720 660 760 560 700 890 850 1020 660 870 ...
# $ Final.Volume..uL.                   : int  620 560 660 460 600 790 750 920 560 770 ...
# $ TTP_MBNI_Nanodrop_RNAConc..ng.uL.   : num  528 398 452 523 463 ...
# $ TTP_MBNI_Nanodrop_260.280           : num  2 1.82 1.77 1.96 1.78 1.8 1.78 1.79 1.8 1.8 ...
# $ TTP_MBNI_Nanodrop_260.230           : num  1.87 2.02 1.99 1.99 1.9 2.04 2.24 1.94 1.69 1.89 ...
# $ TTP_DSC_Nanodrop_RNAConc..ng.uL.    : num  631 429 555 644 551 ...
# $ TTP_DSC_Nanodrop_260.280            : num  2.02 1.92 2.01 2.02 2.03 1.91 1.97 1.98 2.01 1.99 ...
# $ TTP_DSC_Nanodrop_260.230            : num  1.83 2.04 2.05 1.96 1.9 2.09 2.29 1.95 1.62 1.86 ...
# $ TTP_BioAnalyzer_RIN                 : num  7.4 7.8 7.7 7.7 8 7.9 7 7.2 8.2 7.5 ...
# $ TTP_BioAnalyzer_rRNA.Ratio..28s.18s.: num  2.2 2.2 2.8 2.1 2.4 2.5 1.9 2.3 2.3 2.3 ...
# $ TTP_Average.RNAConc..ng.uL.         : num  579 414 503 583 507 ...
# $ TTP_Average.260.280                 : num  2.01 1.87 1.89 1.99 1.91 ...
# $ TTP_Average.260.230                 : num  1.85 2.03 2.02 1.98 1.9 ...

SubjectInfo$Diagnosis<-as.factor(SubjectInfo$Diagnosis)
levels(SubjectInfo$Diagnosis)
SubjectInfo$Diagnosis<-relevel(SubjectInfo$Diagnosis, ref="Control")
levels(SubjectInfo$Diagnosis)
#[1] "Control" "BP"      "Schiz" 

SubjectInfo$Gender<-as.factor(SubjectInfo$Gender)
levels(SubjectInfo$Gender)
SubjectInfo$Gender<-relevel(SubjectInfo$Gender, ref="M")
levels(SubjectInfo$Gender)
#[1] "M" "F"

SubjectInfo$Dissecton.Group<-as.factor(SubjectInfo$Dissecton.Group)
levels(SubjectInfo$Dissecton.Group)
#[1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12"

SubjectInfo$ID<-as.character(SubjectInfo$ID)


#################################


#Dealing With Missing MetaData:


#Note: One of the subjects is missing pH. I'm going to use average pH for right now, but that is a problem that should really be remedied if we want to control for pH as a co-variate.
#I double checked with the 2019 version of the Pritzker database that I have, and it is missing there too. I'll need to double check on the server, and then probably contact UCIrvine.
#For now, I went back to the original .csv document and replaced with NA (it was written as 0!), and then in the code below replaced that NA with mean pH.

SubjectInfo[is.na(SubjectInfo$pH),]
# Barcode Subject.Number   Cohort Diagnosis Age Gender pH AFS Hours.Final Slab.Format Slab.Number Dissecton.Group ID TTP.ID
# 35 B010783A           4463 Cohort 7     Schiz  31      M NA   0        23.5        SLAB           1               6 75       
# Block.Weight..g. TZP.ID MBNI_Nanodrop_RNAConc..ng.uL. MBNI_Nanodrop_260.280 MBNI_Nanodrop_260.230 DSC_Nanodrop_RNAConc..ng.uL.
# 35               NA                                   NA                    NA                    NA                           NA
# DSC_Nanodrop_260.280 DSC_Nanodrop_260.230 TZP_BioAnalyzer_RIN TZP_BioAnalyzer_rRNA.Ratio..28s.18s.TZP_Average.RNAConc..ng.uL.
# 35                   NA                   NA              NA                               NA                      NA
#TZP_Average.260.280TZP_Average.260.230
# 35              NA              NA

SubjectInfo$MissingPH<-is.na(SubjectInfo$pH)

SubjectInfo$pH[SubjectInfo$MissingPH]<-mean(SubjectInfo$pH, na.rm=T)

#Interestingly, this same subject is one of the two subjects missing RNA quality (RIN/RNAconc) information (i.e., one of the two subjects where the RNAextraction was re-done):

SubjectInfo[is.na(SubjectInfo$Block.Weight..g.),]
# Barcode Subject.Number    Cohort Diagnosis Age Gender       pH AFS Hours.Final Slab.Format Slab.Number Dissecton.Group ID TTP.ID
# 34 B008974A           5162 Cohort 12        BP  44      M 6.780000   0        34.3        SLAB           1               6 74       
# 35 B010783A           4463  Cohort 7     Schiz  31      M 6.809296   0        23.5        SLAB           1               6 75       
# Block.Weight..g. TZP.ID MBNI_Nanodrop_RNAConc..ng.uL. MBNI_Nanodrop_260.280 MBNI_Nanodrop_260.230 DSC_Nanodrop_RNAConc..ng.uL.
# 34               NA                                   NA                    NA                    NA                           NA
# 35               NA                                   NA                    NA                    NA                           NA
# DSC_Nanodrop_260.280 DSC_Nanodrop_260.230 TZP_BioAnalyzer_RIN TZP_BioAnalyzer_rRNA.Ratio..28s.18s.TZP_Average.RNAConc..ng.uL.
# 34                   NA                   NA              NA                               NA                      NA
# 35                   NA                   NA              NA                               NA                      NA
#TZP_Average.260.280TZP_Average.260.230
# 34              NA              NA
# 35              NA              NA

SubjectInfo$MissingExtractionData<-is.na(SubjectInfo$Block.Weight..g.)


#Note - it occurs to me that these subjects should probably also have their dissection group changed, because they needed re-extraction
levels(SubjectInfo$Dissecton.Group) <- c(levels(SubjectInfo$Dissecton.Group), "6redone")
SubjectInfo$Dissecton.Group[SubjectInfo$MissingExtractionData]<-"6redone"

table(SubjectInfo$Dissecton.Group)
# 1       2       3       4       5       6       7       8       9      10      11      12 6redone 
# 6       5       6       6       6       4       6       6       6       6       6       7       2 

#I'm going to replace those missing measurements with mean measurements also:

SubjectInfo$Block.Weight..g.[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$Block.Weight..g., na.rm=T)

SubjectInfo$TTP_BioAnalyzer_RIN[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TTP_BioAnalyzer_RIN, na.rm=T)

SubjectInfo$TZP_BioAnalyzer_RIN[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TZP_BioAnalyzer_RIN, na.rm=T)

SubjectInfo$TTP_BioAnalyzer_rRNA.Ratio..28s.18s.[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TTP_BioAnalyzer_rRNA.Ratio..28s.18s., na.rm=T)  

SubjectInfo$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TZP_BioAnalyzer_rRNA.Ratio..28s.18s., na.rm=T)  

SubjectInfo$TTP_Average.RNAConc..ng.uL.[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TTP_Average.RNAConc..ng.uL., na.rm=T)

SubjectInfo$TZP_Average.RNAConc..ng.uL.[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TZP_Average.RNAConc..ng.uL., na.rm=T)

SubjectInfo$TTP_Average.260.280[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TTP_Average.260.280, na.rm=T)

SubjectInfo$TZP_Average.260.280[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TZP_Average.260.280, na.rm=T)

SubjectInfo$TTP_Average.260.230[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TTP_Average.260.230, na.rm=T)

SubjectInfo$TZP_Average.260.230[SubjectInfo$MissingExtractionData]<-mean(SubjectInfo$TZP_Average.260.230, na.rm=T)

#I didn't bother replacing the missing data in the individual MBNI and DSC measurements


SubjectInfo[is.na(SubjectInfo$Hours.Final),]
#<0 rows> (or 0-length row.names)

SubjectInfo[is.na(SubjectInfo$Age),]
#<0 rows> (or 0-length row.names)

#No other commonly missing co-variates


#I found what I believe is a decimal place error in Evan's rRNA Ratio values (28s/18s) (14.4 vs. 1.44). I replaced it with the corrected value. 
SubjectInfo$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.[SubjectInfo$TZP_BioAnalyzer_rRNA.Ratio..28s.18s.>14]<-1.44

#########################

#########################

#QC based on basic sample quality metrics:


#Outputting histograms for the original full dataset:

setwd("~/Documents/Microarray Gen/FrontalPole/Output/Histograms/Histograms_FullSample")

colnames(SubjectInfo)

NumericColumnsInSubjectInfo<-c(5,7,9,15,17:27,29:41)

Temp<-SubjectInfo[,NumericColumnsInSubjectInfo]

for(i in c(1:ncol(Temp))){
  pdf(paste("Histogram_", colnames(Temp)[i], "_FullSample.pdf", sep=""), width=4, height=4)
  hist(Temp[,i], xlab=colnames(Temp)[i], main="", col="grey", breaks=20)
  dev.off()
}


#There were three subjects that I ended up tossing out because their RNA quality metrics were strikingly bad in the final (purified) TZP sample: really low RIN in two subjects (<3) and a really low 260/230 purity (<1) as measured by two different nanodrops (MBNI & DSC) accompanied by unusually low concentration.
#Notably, in earlier versions of our analyses, these subjects were also thrown out for having outlier gene expression measurements, typically in both replicates (e.g., very extreme PCA scores)
#Also note that we did not have detailed RNA extraction data from two subjects - we will have to depend on basic QC to determine if something is wrong for them.


setwd("~/Documents/Microarray Gen/FrontalPole/Output/Histograms/Histograms_AfterBasicQC")

Temp<-SubjectInfo[SubjectInfo$LowQualityRNA==F,NumericColumnsInSubjectInfo]

for(i in c(1:ncol(Temp))){
  pdf(paste("Histogram_", colnames(Temp)[i], "_AfterBasicQC.pdf", sep=""), width=4, height=4)
  hist(Temp[,i], xlab=colnames(Temp)[i], main="", col="grey", breaks=20)
  dev.off()
}

setwd("~/Documents/Microarray Gen/FrontalPole/Sample_MetaData")

SubjectInfo$LowQualityRNA<-(SubjectInfo$TZP_BioAnalyzer_RIN<3 | SubjectInfo$TZP_Average.260.230<1) 

write.csv(SubjectInfo, "SubjectInfo_FollowingBasicQC.csv")

#For looking at demographics in more detail:

SubjectInfo_SurvivedBasicQC<-SubjectInfo[SubjectInfo$LowQualityRNA==FALSE,]
SubjectInfo_SurvivedBasicQC_NoMissingExtractionData<-SubjectInfo[(SubjectInfo$LowQualityRNA|SubjectInfo$MissingExtractionData)==FALSE,]

############################

