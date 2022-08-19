#170_BothDatasets_GeneVsDrugPlots
#Megan Hagenauer, Evelyn Richardson
#2021-06-22
#updated 2021-10-07 to fix the ambiguous records about suicide vs. overdose, and to add in ordinal ratings for presence/absence of drugs based on toxicology reports (+1 = present, 0= unknown, -1=absent)
#At that point I also changed the gene list to remove MAOB
#And removed a few variables that had too small of a sample size to be meaningful (<3):  easily_distracted, abdominal discomfort, mood stabilizers
#And added interactions with the legal system (because of HTR2B's story in Finland)


#I did this code first:

#Set working directory to folder with relevant data files:
#E.g., on Megan's computer:
setwd("~/Documents/Microarray Gen/FrontalPole/Output/Intermediate_forDrugPlots")
list.files()
# [1] "All info subjects in cards_Drugs_MaybeIsYes_V2wOverdose.csv"                
# [2] "All info subjects in cards_Drugs_MaybeIsYes.csv"                            
# [3] "DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2.csv"                                
# [4] "GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2.csv"                              
# [5] "Original"                                                                   
# [6] "SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes.csv"  
# [7] "SubjectInfo_OrderedForDA5HTCqMatrix_QCed3.csv"                              
# [8] "SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes.csv"
# [9] "SubjectInfo_OrderedForGabaGluCqMatrix_QCed3.csv"                            
# [10] "Workspace_ForFrontalPolePCRDrugPlots_2021.RData" 

#Read in the relevant objects:

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3<-read.csv("SubjectInfo_OrderedForDA5HTCqMatrix_QCed3.csv", header=TRUE, stringsAsFactors = FALSE)

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3<-read.csv("SubjectInfo_OrderedForGabaGluCqMatrix_QCed3.csv", header=TRUE, stringsAsFactors = FALSE)

#Original Analysis
#DrugsAndSymptoms_MaybeIsYes<-read.csv("All info subjects in cards_Drugs_MaybeIsYes.csv", header=TRUE, stringsAsFactors = FALSE)


DrugsAndSymptoms_MaybeIsYes<-read.csv("All info subjects in cards_Drugs_MaybeIsYes_V2wOverdose.csv", header=TRUE, stringsAsFactors = FALSE)

#Double checking something weird:
table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$suicidal_ideations)
#         no yes
# BP       2  19
# Control 26   1
# Schiz   15   9

table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$manner_of_death)
#           Accident Sudden medical condition Suicide Undetermined Undetermined/Accidental overdose
# BP             4                        4      11            1                                1
# Control        2                       25       0            0                                0
# Schiz          3                       11       6            0                                4


library(plyr)

str(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
str(DrugsAndSymptoms_MaybeIsYes)

colnames(DrugsAndSymptoms_MaybeIsYes)[1]<-"Subject.Number"
  
SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes<-join(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, DrugsAndSymptoms_MaybeIsYes, by="Subject.Number", type="left") 

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes<-join(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, DrugsAndSymptoms_MaybeIsYes, by="Subject.Number", type="left") 

#After thinking about the toxicology-based variables, I realized that it was a little silly to place any CTRL subjects in the unknown (0) category for antipsychotics and antidepressants (opioids are still a reasonable possibility, given their clinical usage)  

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Tox_Antipsychotics[SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Diagnosis=="Control"]<-(-1)

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Tox_Antidepressants[SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Diagnosis=="Control"]<-(-1)

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Tox_Antipsychotics[SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Diagnosis=="Control"]<-(-1)

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Tox_Antidepressants[SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes$Diagnosis=="Control"]<-(-1)


write.csv(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes, "SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes.csv")

write.csv(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes, "SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes.csv")

 

####################


#Code for actually making some plots and running exploratory analyses:

#Set working directory to folder with relevant data files:
#E.g., on Megan's computer:
setwd("~/Documents/Microarray Gen/FrontalPole/Output/Intermediate_forDrugPlots")
list.files()

#Read in the relevant qPCR data objects:

DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2<-read.csv("DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2.csv", header=TRUE, stringsAsFactors = FALSE)

GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2<-read.csv("GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2.csv", header=TRUE, stringsAsFactors = FALSE)

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug<-read.csv("SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes.csv")

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug<-read.csv("SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_DrugsAndSymptoms_MaybeIsYes.csv")


#Only some of those variables need re-coding for running an exploratory analysis:

levels(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis)
#[1] "BP"      "Control" "Schiz"
SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis<-relevel(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, ref="Control")
levels(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis)
#[1] "Control" "BP"      "Schiz" 

levels(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Diagnosis)
#[1] "BP"      "Control" "Schiz"  
SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Diagnosis<-relevel(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Diagnosis, ref="Control")
levels(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Diagnosis)
#[1] "Control" "BP"      "Schiz" 

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$provisional_axis_i_substance_related) 
#this probably would need to be consolidated in some way - several categories are too small to be useful. Probably better to use Evelyn's analysis.

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$psychosis_freq_or_dur) 
#               Associated with Substance Use Chronic Cyclical Not applicable Single episode
# Control                             0       0        0             60              0
# BP                                  2       0       18             20              6
# Schiz                               0      44        4              0              0

#This once could be interesting - the Associated with Substance Use category is too small though.

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Psychosis<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$psychosis_freq_or_dur!="Not applicable"

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Psychosis<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$psychosis_freq_or_dur!="Not applicable"


table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Psychosis)
#         FALSE TRUE
# Control    60    0
# BP         20   26
# Schiz       0   48

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$negative_sxs)
#         no yes
# Control 60   0
# BP      46   0
# Schiz   16  32

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$affective_freq_or_dur)
#           Chronic Cyclical Not applicable Single episode
# Control  0       0        0             60              0
# BP       0      26       18              2              0
# Schiz    2      28        6              8              4

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Affective<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$affective_freq_or_dur!="Not applicable"& SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$affective_freq_or_dur!=""

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Affective<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$affective_freq_or_dur!="Not applicable"& SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$affective_freq_or_dur!=""


table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Affective)
#         FALSE TRUE
# Control    60    0
# BP          2   44
# Schiz      10   38

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$anxiety_freq_or_dur)

#           Chronic Cyclical Not applicable Single episode
# Control       0        0             60              0
# BP            8       16             22              0
# Schiz        24        4             18              2

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Anxiety<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$anxiety_freq_or_dur!="Not applicable"

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Anxiety<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$anxiety_freq_or_dur!="Not applicable"

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Anxiety)

#           FALSE TRUE
# Control    60    0
# BP         22   24
# Schiz      18   30


table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$syndrome_at_time_of_death)  
#               Depressed Intoxicated None indicated Psychotic
# Control  0         0           0             60         0
# BP       6        28           6              6         0
# Schiz   14         6           0             26         2

#hmm... there are a lot of unknowns there, but depressed at time of death might be worth noting... although that may just overlap with suicide.

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$manner_of_death)  
# Accident Sudden medical condition Suicide Undetermined Undetermined/Accidental overdose
# Control        4                       56       0            0                                0
# BP             8                        8      26            2                                2
# Schiz          4                       26      10            0                                8

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$syndrome_at_time_of_death, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$manner_of_death)  
#               Accident Sudden medical condition Suicide Undetermined Undetermined/Accidental overdose
#                       2                        8       4            2                                4
# Depressed             0                        0      32            0                                2
# Intoxicated           4                        2       0            0                                0
# None indicated       10                       78       0            0                                4
# Psychotic             0                        2       0            0                                0
#Yes, looks like basically all depressed at TOD are suicide or overdose

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$syndrome_at_time_of_death, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Suicide) 

#                 0  1
#                 16  4
# Depressed       2 32
# Intoxicated     6  0
# None indicated 92  0
# Psychotic       2  0

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$syndrome_at_time_of_death, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Overdose) 
#                 0  1
#                 12  8
# Depressed      24 10
# Intoxicated     2  4
# None indicated 82 10
# Psychotic       2  0

# Old analysis:
# SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Suicide<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$manner_of_death=="Suicide"
# 
# SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$Suicide<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug$manner_of_death=="Suicide"


#Other true/false variables that might be interesting but are already coded properly:
# $ delusions
# $ hallucinations
# [69] "fatigue"
# [71] "trouble_concentrating" 
# [79] "reckless" 
# "irritable" 
# "syndrome_at_time_of_death"
# [109] "SubstancesOrNot"                                           "Antipsychotics"                                           
# [111] "Antidepressants..Serotonergic.Norandrenergic.."            "Stimulants..Dopaminergic.."                               
# [113] "Mood.Stabilizers..Lithium."                                "Opioids"                                                  
# [115] "Alcohol"                                                   "Depressants..Benzodiazepines."                            
# [117] "Cannabinoid"                                               "Hallucinogens"                                            
# [119] "Tobacco"    


table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$suicidal_ideations)
#         no yes
# Control 56   4
# BP       4  42
# Schiz   34  14
#Something isn't right. There is only 1 control subject with suicidal ideation.
#Ah - it must be one of the subjects represented by 4 samples.

table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$Diagnosis, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug$diagnosis)
#         BP Control Schiz
# Control  0      60     0
# BP      46       0     0
# Schiz    0       0    48



#Which genes to emphasize:

GenesWFDR10<-c("HTR2B", "DRD4", "SST", "ABAT", "GPHN", "MAPK1")

colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug)

#removed due to small sample size:
#easily_distracted, abdominal discomfort, mood stabilizers

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug[,c(7, 50:54, 61:65, 68:80, 82:83, 86:93, 95:96, 104, 113:116, 118:126)]

colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)

SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory<-SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug[,c(7, 53:57, 64:68, 71:83, 85:86, 89:96, 98:99, 107, 116:119, 121:129)]

colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)

str(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)



#*****************************

#I would like to quickly examine the multicollinearity between these variables - it would be easier if they were coded as 0/1

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric<-SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory

# SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric[SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric=="yes"]<-"1"
# SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric[SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric=="no"]<-"0"

str(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric)

TempBinaryVar<-matrix(as.character(as.matrix(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,c(7:47)])), 154, 41)

TempBinaryVar[TempBinaryVar=="yes"]<-"1"
TempBinaryVar[TempBinaryVar=="no"]<-"0"
str(TempBinaryVar)

TempBinaryVar2<-matrix(as.numeric(TempBinaryVar), 154,41)

colnames(TempBinaryVar2)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,c(7:47)])

str(TempBinaryVar2)

TempLogicalVar<-matrix(as.numeric(as.matrix(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,c(48:50)])), 154, 3)

str(TempLogicalVar)
colnames(TempLogicalVar)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,c(48:50)])

TempIntegerVar<-matrix(as.numeric(as.matrix(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,c(2:6)])), 154, 5)

colnames(TempIntegerVar)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,c(2:6)])

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric<-cbind(TempIntegerVar,TempBinaryVar2,TempLogicalVar)

#Well that was a stupid amount of effort. Now I have a matrix of numeric binary variables, but there are still multiple entries for each subject. I probably should have just worked with the original subject database. Meh. 

#Ah well, good enough for getting a preliminary taste of how these variables travel together.

write.csv(cor(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric), "CorrelationMatrix_SubjectInfoForDA5HT.csv")
#That's actually pretty interesting. I guess I should put in the extra legwork to reduce the redundancies by subject.
#It might be worthwhile to throw in some of the other post-mortem subject variables as well - might be interesting. 
#And I should probably make a column with binary diagnosis.

SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric

str(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)
#Numeric Variables worth adding:
#Make numeric version: 6, 8

SubjectInfo_For_Drug_ForExploratory_NumericMore<-cbind(as.matrix(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3[,c(7, 9,11, 16, 24:28)]), SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory_Numeric, as.numeric(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis=="BP"), as.numeric(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Diagnosis=="Schiz"), as.numeric(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Gender=="F"))

colnames(SubjectInfo_For_Drug_ForExploratory_NumericMore)[59]<-"Bipolar"  

colnames(SubjectInfo_For_Drug_ForExploratory_NumericMore)[60]<-"Schiz"

colnames(SubjectInfo_For_Drug_ForExploratory_NumericMore)[61]<-"Female"                  
                                                                                                                   #Sorry the naming in this code is ridiculously stupid because I ended up changing what I was doing midflow
   #Probably not any dumber than the order that I am doing things though LOL

SubjectIDs<-names(table(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Subject.Number))
SubjectInfo_For_Drug_ForExploratory_Numeric<-matrix(0, length(SubjectIDs), 61)
row.names(SubjectInfo_For_Drug_ForExploratory_Numeric)<-SubjectIDs
colnames(SubjectInfo_For_Drug_ForExploratory_Numeric)<-colnames(SubjectInfo_For_Drug_ForExploratory_NumericMore)

for(i in c(1:ncol(SubjectInfo_For_Drug_ForExploratory_NumericMore))){
  SubjectInfo_For_Drug_ForExploratory_Numeric[,i]<-tapply(SubjectInfo_For_Drug_ForExploratory_NumericMore[,i], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3$Subject.Number, function(y) mean(y, na.rm=TRUE))
}

str(SubjectInfo_For_Drug_ForExploratory_Numeric)
# num [1:69, 1:61] 18 23 58 26 45 53 36 42 56 36 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:69] "2169" "2311" "2316" "2466" ...
# ..$ : chr [1:61] "Age" "pH" "Hours.Final" "Block.Weight..g." ...

write.csv(SubjectInfo_For_Drug_ForExploratory_Numeric, "SubjectInfo_For_Drug_ForExploratory_NumericBinary.csv")

write.csv(cor(SubjectInfo_For_Drug_ForExploratory_Numeric), "CorrelationMatrix_SubjectInfo_For_Drug_ForExploratory_NumericBinary.csv")

pdf("Heatmap_CorrelationMatrix_SubjectInfo_For_Drug_ForExploratory_NumericBinary.pdf", height=10, width=10)
heatmap(cor(SubjectInfo_For_Drug_ForExploratory_Numeric))
dev.off()

write.csv(cor(SubjectInfo_For_Drug_ForExploratory_Numeric[SubjectInfo_For_Drug_ForExploratory_Numeric[,59]==1,]), "CorrelationMatrix_SubjectInfo_For_Drug_ForExploratory_NumericBinary_Bipolar.csv")


write.csv(cor(SubjectInfo_For_Drug_ForExploratory_Numeric[SubjectInfo_For_Drug_ForExploratory_Numeric[,60]==1,]), "CorrelationMatrix_SubjectInfo_For_Drug_ForExploratory_NumericBinary_Schiz.csv")

pdf("Heatmap_CorrelationMatrix_SubjectInfo_For_Drug_ForExploratory_NumericBinary_BPSchiz.pdf", height=10, width=10)
heatmap(cor(SubjectInfo_For_Drug_ForExploratory_Numeric[SubjectInfo_For_Drug_ForExploratory_Numeric[,60]==1|SubjectInfo_For_Drug_ForExploratory_Numeric[,59]==1,]))
dev.off()

#I should loop sample sizes for all of these variables too:

colnames(SubjectInfo_For_Drug_ForExploratory_Numeric)

Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric<-rep("CTRL", length=nrow(SubjectInfo_For_Drug_ForExploratory_Numeric))
Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric[SubjectInfo_For_Drug_ForExploratory_Numeric[,59]==1]<-"Bipolar"
Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric[SubjectInfo_For_Drug_ForExploratory_Numeric[,60]==1]<-"Schiz"

Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric
col 10-58 skip 12-14

Diagnosis_SampleSizes_ForExploratoryVariables<-matrix(0,58,3)
colnames(Diagnosis_SampleSizes_ForExploratoryVariables)<-colnames(table(SubjectInfo_For_Drug_ForExploratory_Numeric[,10], Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric))
row.names(Diagnosis_SampleSizes_ForExploratoryVariables)<-colnames(SubjectInfo_For_Drug_ForExploratory_Numeric)[1:58]

#Just grabbing the counts for the "yes"
for(i in c(10:58)){
  Diagnosis_SampleSizes_ForExploratoryVariables[i,]<-table(SubjectInfo_For_Drug_ForExploratory_Numeric[,i], Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric)[2,]
}

#And then cutting out the empty rows while outputting:
write.csv(Diagnosis_SampleSizes_ForExploratoryVariables[c(10:11, 15:58),], "Diagnosis_SampleSizes_ForExploratoryVariables.csv")

#Getting the full sample sizes for effective comparison:
table(Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric)

# Diagnosis_For_SubjectInfo_For_Drug_ForExploratory_Numeric
# Bipolar    CTRL   Schiz 
# 21      26      22 
  
#************************
  
  #Differential expression results - graphs first

for(i in c(2:ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory))){

pdf(paste("Boxplot_DA5HT_HTR2B_ByDiagnosis_", colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4), ylim=c(min(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]), max(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"])), main="HTR2B", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("HTR2B Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf(paste("Boxplot_DA5HT_DRD4_ByDiagnosis_", colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4), ylim=c(min(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"]), max(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"])),  main="DRD4", outline=FALSE)
stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("DRD4 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

#In original:
# pdf(paste("Boxplot_DA5HT_MAOB_ByDiagnosis_", colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
# boxplot(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOB"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4), main="MAOB", outline=FALSE)
# stripchart(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOB"]~SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
#            method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
# mtext(expression(paste("MAOB Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
# dev.off()

}

#I need to come back to this - my boxplots are getting different yaxes because outline=FALSE

for(i in c(2:ncol(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory))){
  
pdf(paste("Boxplot_GabaGlu_ABAT_ByDiagnosis_", colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4),  ylim=c(min(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"]), max(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"])), main="ABAT", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("ABAT Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()


pdf(paste("Boxplot_GabaGlu_SST_ByDiagnosis_", colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4),  ylim=c(min(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"]), max(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"])), main="SST", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("SST Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf(paste("Boxplot_GabaGlu_GPHN_ByDiagnosis_", colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4),  ylim=c(min(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"]), max(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"])), main="GPHN", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("GPHN Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

pdf(paste("Boxplot_GabaGlu_MAPK1_ByDiagnosis_", colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)[i], ".pdf", sep=""), width=10, height=6)
boxplot(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, xlab="Diagnosis", cex.axis=1.3, cex.lab=1.3, pch=20, cex=1.7, col=c(2,2,3,3,4,4),  ylim=c(min(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"]), max(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"])), main="MAPK1", outline=FALSE)
stripchart(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"]~SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory[,i]+SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory$Diagnosis, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 21, cex=1, col = 'black')
mtext(expression(paste("MAPK1 Log(2) Fold Change (-", Delta, Delta, "Cq)", sep="")), cex=1.3, side=2, line=2.5)
dev.off()

}


#I should loop some statistical output:

#GabaGlu dataset:

ncol(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
#[1] 50

#Test:
library(lme4)
i<-2
Temp<-data.frame(y=GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
Model<-lmer(y~Diagnosis+Temp[,51+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)

ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
colnames(ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
colnames(ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
row.names(ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)

for(i in c(2:ncol(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="ABAT"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
  Model<-lmer(y~Diagnosis+Temp[,51+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}
write.csv(ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
colnames(SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
colnames(SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
row.names(SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)

for(i in c(2:ncol(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="SST"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
  Model<-lmer(y~Diagnosis+Temp[,51+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}
write.csv(SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
colnames(MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
colnames(MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
row.names(MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)

for(i in c(2:ncol(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAPK1"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
  Model<-lmer(y~Diagnosis+Temp[,51+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}
write.csv(MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
colnames(GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
colnames(GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
row.names(GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)

for(i in c(2:ncol(SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)=="GPHN"], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
  Model<-lmer(y~Diagnosis+Temp[,51+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}
write.csv(GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


#DA5HT candidate genes:

ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3)

HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
colnames(HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
colnames(HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
row.names(HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)

for(i in c(2:ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="HTR2B"], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
  Model<-lmer(y~Diagnosis+Temp[,48+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}
write.csv(HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
colnames(DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
colnames(DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
row.names(DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)

for(i in c(2:ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="DRD4"], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
  Model<-lmer(y~Diagnosis+Temp[,48+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
  rm(Temp)
  rm(Model)
}
write.csv(DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


#From original analysis:

# MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 50, 10)
# colnames(MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-row.names(car::Anova(Model, type="III"))
# colnames(MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"ExploratoryVariable"
# row.names(MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
# 
# for(i in c(2:ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory))){
#   Temp<-data.frame(y=DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)=="MAOB"], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
#   Model<-lmer(y~Diagnosis+Temp[,48+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
#   MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i,]<-car::Anova(Model, type="III")[[3]]
#   rm(Temp)
#   rm(Model)
# }
# write.csv(MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")

#For interpreting p-values:
#We are only able to observe effects that are the same in BP/Schiz (and/or massive enough to still average out across both diagnosis groups as an effect)
#Some of these variables do not have an adequate sample size - look at the plots
#For considering multiple comparisons correction:
#49 exploratory variables
#6 genes
6*49
#[1] 294
0.05/294
#[1] 0.000170068
#If we went straight Bonferroni-style correction, we would only consider results with <0.000170068
#If we're willing to tolerate more false-detection, we could double this: <0.00034
#I am inclined to pay more attention to what looks suggestive (e.g., p<0.01), see whether it seems to "explain" or "reverse" the diagnosis effect, and then determine whether there is additional evidence to support the relationship (either documented relationship or evidence from Iwamoto et al.), and then still consider things exploratory.

#I may just run a formal FDR correction, even though it really is *still* very exploratory - these variables are, by definition, collinear because they are all related to diagnosis.

TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-rbind.data.frame(ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1, c(2:3)], SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)])


#Original code:
# TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-rbind.data.frame(ABAT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1, c(2:3)], SST_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], GPHN_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], MAPK1_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], DRD4_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], HTR2B_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)], MAOB_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-1,c(2:3)])

row.names(TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)
#R added numbers after each row.name when it repeats. 

str(TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)
# 'data.frame':	294 obs. of  2 variables:
#   $ Diagnosis          : num  0.02218 0.00137 0.00116 0.00113 0.00173 ...
# $ ExploratoryVariable: num  0.741 0.4585 0.0874 0.3091 0.3738 ...

TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval$Gene<-c(rep("ABAT", 49), rep("SST", 49), rep("GPHN", 49), rep("MAPK1", 49), rep("DRD4", 49), rep("HTR2B", 49))

#Calculating FDR:
library(multtest)

colnames(TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)

TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval$ExploratoryVar_FDR<-TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,2]

  tempPvalAdj<-mt.rawp2adjp(TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,2], proc="BH")
  TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval$ExploratoryVar_FDR<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]


write.csv(TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_PvalwFDR.csv")

#Here are the relationships that survive (somewhat sketchy) false detection correction:

TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[TopGenes_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval$ExploratoryVar_FDR<0.10,]

#                               Diagnosis ExploratoryVariable  Gene ExploratoryVar_FDR
# Overdose1                   9.380086e-02        8.785656e-04   SST       3.689975e-02
# Opioids1                    1.105059e-02        1.456939e-04   SST       1.372682e-02
# aggitated3                  2.548195e-01        1.145413e-04 MAPK1       1.372682e-02
# Stimulants..Dopaminergic..3 2.012785e-01        2.430972e-03 MAPK1       8.933822e-02
# Hallucinogens3              5.102575e-03        2.334494e-04 MAPK1       1.372682e-02
# Opioids4                    1.899793e-02        2.081675e-04  DRD4       1.372682e-02
# aggitated5                  2.087423e-02        8.660248e-04 HTR2B       3.689975e-02
# Opioids5                    9.936015e-05        9.419662e-09 HTR2B       2.769381e-06

#Opioids pops up a lot - what is the sample size for that?

table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$Opioids)
#         no yes
# BP      15   6
# Control 27   0
# Schiz   21   3
#Small, but not negligable. Interesting. (note that this doesn't take into account samples removed for low quality from the dataset - I should probably double-check that)

#Aggitated pops up a lot too:
table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$aggitated)
#         no yes
# BP      12   9
# Control 27   0
# Schiz   16   8
#Pretty decent.

#Easily distracted is worth ignoring:
# table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$easily_distracted)
#         no yes
# BP      20   1
# Control 27   0
# Schiz   24   0

table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$Stimulants..Dopaminergic..)
#         no yes
# BP      11  10
# Control 27   0
# Schiz   19   5
#Solid sample size too.

table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$Hallucinogens)
#         no yes
# BP      19   2
# Control 27   0
# Schiz   23   1
#Probably too small of a sample size to be meaningful, but since there is also a relationship with stimulants, it may just indicate a general relationship with drugs of abuse. Did all of the subjects that used hallucinogens also use stimulants?

table(DrugsAndSymptoms_MaybeIsYes$Stimulants..Dopaminergic.., DrugsAndSymptoms_MaybeIsYes$Hallucinogens)
#     no yes
# no  56   1
# yes 13   2
#Not quite.

table(DrugsAndSymptoms_MaybeIsYes$diagnosis, DrugsAndSymptoms_MaybeIsYes$Antidepressants..Serotonergic.Norandrenergic..)
#         no yes
# BP      17   4
# Control 27   0
# Schiz   20   4
#small, but doable.


#Another question: 
#Is it possible that some of these exploratory variables explain diagnosis effects for some of the weaker (nominal) relationships that have replicated across datasets?
#Used genes with nominal relationships in the frontal pole in either the qPCR or microarray meta-analysis + at least one other nominal relationship (frontal pole or DLPFC)

#Note - I came back to this and updated it to include the top genes as well, since they would need to be included in the FDR correction. 

#I also realized this output would be waaaaay more useful if I had betas/SE/tstats to go with it:


WeakerResults_Genes<-c("BDNF", "DRD2", "DRD3", "HOMER1", "SLC6A12", "COMT", "SLC6A11", "GABRB1", "GRM5", "CALB1", "NSF", "SNCA", "MAOB", "GFAP", "HTR2B", "DRD4", "MAPK1", "ABAT", "GPHN", "SST")

length(WeakerResults_Genes)
#[1] 20


WeakerResults_Genes_DA5HT<-colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)[colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)%in%WeakerResults_Genes]
#[1] "COMT"  "DRD2"  "DRD3"  "DRD4"  "HTR2B" "MAOB" 

WeakerResults_Genes_GABAGLU<-colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)[colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)%in%WeakerResults_Genes]
# [1] "ABAT"    "SST"     "BDNF"    "CALB1"   "GABRB1"  "GPHN"    "GRM5"    "HOMER1"  "GFAP"    "MAPK1"   "NSF"    
# [12] "SLC6A11" "SLC6A12" "SNCA" 

20*50
#[1] 1000
6*50
#[1] 300
14*50
#[1] 700


WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 300, 10)
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-paste("AnovaPval_", row.names(car::Anova(Model, type="III")), sep="")
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"AnovaPval_ExploratoryVariable"
row.names(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 6)

WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta<-matrix(0, 300, 10)
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta)<-paste("LMBeta_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta)[4]<-"LMBeta_ExploratoryVariable"
row.names(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 6)

WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE<-matrix(0, 300, 10)
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE)<-paste("LMSE_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE)[4]<-"LMSE_ExploratoryVariable"
row.names(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 6)

WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat<-matrix(0, 300, 10)
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat)<-paste("LMTstat_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat)[4]<-"LMTstat_ExploratoryVariable"
row.names(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 6)

WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval<-matrix(0, 300, 10)
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)<-paste("LMpval_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)[4]<-"LMpval_ExploratoryVariable"
row.names(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 6)


for(j in c(1:length(WeakerResults_Genes_DA5HT))){
for(i in c(2:ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory))){
  Temp<-data.frame(y=DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)==WeakerResults_Genes_DA5HT[j]], SubjectInfo_OrderedForDA5HTCqMatrix_QCed3, SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory)
  Temp$Diagnosis<-relevel(as.factor(Temp$Diagnosis), ref="Control")
  Model<-lmer(y~Diagnosis+Temp[,48+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
  WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i+((j-1)*50),]<-car::Anova(Model, type="III")[[3]]
  WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),1]
  WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),2]
  WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),3]
  WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),4]
  rm(Temp)
  rm(Model)
}
}

temp<-vector(mode="character", length=0)

for(i in c(1:length(WeakerResults_Genes_DA5HT))){
temp<-c(temp,rep(WeakerResults_Genes_DA5HT[i], 50))
}

str(temp)


WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-data.frame(Gene=temp, ExploratoryVariable=row.names(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval), WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta, WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE, WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat, WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)

write.csv(WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")

#I need to go in and remove the empty rows for the diagnosis model before running FDR corrections


WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-matrix(0, 700, 10)
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-paste("AnovaPval_", row.names(car::Anova(Model, type="III")), sep="")
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)[3]<-"AnovaPval_ExploratoryVariable"
row.names(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 14)

WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta<-matrix(0, 700, 10)
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta)<-paste("LMBeta_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta)[4]<-"LMBeta_ExploratoryVariable"
row.names(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 14)

WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE<-matrix(0, 700, 10)
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE)<-paste("LMSE_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE)[4]<-"LMSE_ExploratoryVariable"
row.names(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 14)

WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat<-matrix(0, 700, 10)
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat)<-paste("LMTstat_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat)[4]<-"LMTstat_ExploratoryVariable"
row.names(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 14)

WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval<-matrix(0, 700, 10)
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)<-paste("LMpval_", row.names(summary.lm(lm(Model))$coefficients)[c(1:10)], sep="")
colnames(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)[4]<-"LMpval_ExploratoryVariable"
row.names(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)<-rep(colnames(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory), 14)


for(j in c(1:length(WeakerResults_Genes_GABAGLU))){
  for(i in c(2:ncol(SubjectInfo_OrderedForDA5HTCqMatrix_QCed3_Drug_ForExploratory))){
    Temp<-data.frame(y=GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)==WeakerResults_Genes_GABAGLU[j]], SubjectInfo_OrderedForGabaGluCqMatrix_QCed3, SubjectInfo_OrderedForGabaGluCqMatrix_QCed3_Drug_ForExploratory)
    Temp$Diagnosis<-relevel(as.factor(Temp$Diagnosis), ref="Control")
    Model<-lmer(y~Diagnosis+Temp[,51+i]+Age+Gender+pH+Hours.Final+TZP_BioAnalyzer_RIN+TZP_Average.RNAConc..ng.uL.+Card + (1 | ID), data = Temp, REML=F)
    WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[i+((j-1)*50),]<-car::Anova(Model, type="III")[[3]]
    WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),1]
    WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),2]
    WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),3]
    WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval[i+((j-1)*50),]<-summary.lm(lm(Model))$coefficients[c(1:10),4]
    rm(Temp)
    rm(Model)
  }
}

temp<-vector(mode="character", length=0)

for(i in c(1:length(WeakerResults_Genes_GABAGLU))){
  temp<-c(temp,rep(WeakerResults_Genes_GABAGLU[i], 50))
}

str(temp)

WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-data.frame(Gene=temp, ExploratoryVariable=row.names(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval), WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Beta, WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_SE, WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Tstat, WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_LMpval)

write.csv(WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


WeakerResults_GABAGLU_MLM_ForBind<-WeakerResults_GABAGLU_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(1,51,101,151,201,251,301,351,401,451, 501, 551, 601, 651),]   

WeakerResults_DA5HT_MLM_ForBind<-WeakerResults_DA5HT_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[-c(1,51,101,151, 201, 251),]

colnames(WeakerResults_GABAGLU_MLM_ForBind)

TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval<-rbind.data.frame(WeakerResults_GABAGLU_MLM_ForBind,WeakerResults_DA5HT_MLM_ForBind)

#Calculating FDR:
library(multtest)

colnames(TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval)

TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval$ExploratoryVar_FDR<-TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,5]

tempPvalAdj<-mt.rawp2adjp(TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval[,5], proc="BH")
TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval$ExploratoryVar_FDR<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]


write.csv(TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval, "TopGenes_AndWeakerResults_MLM_Exploratory_ByTheUsualSuspectsAndCard_NegDeltaCq_Pval.csv")


#************

#Making scatterplots for the DrugBudger results:

setwd("~/Documents/Microarray Gen/FrontalPole/Output/Intermediate_forDrugPlots")

DrugBudger_Antipsychotics<-read.csv("DrugBudger_Antipsychotics_forR.csv", header=TRUE, stringsAsFactors = FALSE)
DrugBudger_Corticosteroids<-read.csv("DrugBudger_Corticosteroid_forR.csv", header=TRUE, stringsAsFactors = FALSE)

pdf("DrugBudger_Antipsychotics_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(DrugBudger_Antipsychotics$DrugBudger_Log2FC~DrugBudger_Antipsychotics$Pritzker_qPCR_Bipolar, col=3, pch=16, ylab="DrugBudger: Antipsychotic Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(DrugBudger_Antipsychotics$DrugBudger_Log2FC~DrugBudger_Antipsychotics$Pritzker_qPCR_Bipolar)
abline(BP_Line, col=3, lwd=3)

points(DrugBudger_Antipsychotics$DrugBudger_Log2FC~DrugBudger_Antipsychotics$Pritzker_qPCR_Schiz, col=4, pch=16)
Schiz_Line<-lm(DrugBudger_Antipsychotics$DrugBudger_Log2FC~DrugBudger_Antipsychotics$Pritzker_qPCR_Schiz)
abline(Schiz_Line, col=4, lwd=3)
dev.off()


summary.lm(BP_Line)
# Call:
#   lm(formula = DrugBudger_Antipsychotics$DrugBudger_Log2FC ~ DrugBudger_Antipsychotics$Pritzker_qPCR_Bipolar)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.45019 -0.45801  0.05954  0.50077  2.49676 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                      0.01952    0.05524   0.353  0.72422   
# DrugBudger_Antipsychotics$Pritzker_qPCR_Bipolar -0.68955    0.26409  -2.611  0.00985 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7178 on 167 degrees of freedom
# Multiple R-squared:  0.03922,	Adjusted R-squared:  0.03347 
# F-statistic: 6.818 on 1 and 167 DF,  p-value: 0.009847

summary.lm(Schiz_Line)
# Call:
#   lm(formula = DrugBudger_Antipsychotics$DrugBudger_Log2FC ~ DrugBudger_Antipsychotics$Pritzker_qPCR_Schiz)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.46075 -0.43147  0.06832  0.49515  2.50252 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                    -0.0168     0.0574  -0.293   0.7701  
# DrugBudger_Antipsychotics$Pritzker_qPCR_Schiz  -0.5198     0.1998  -2.602   0.0101 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7179 on 167 degrees of freedom
# Multiple R-squared:  0.03896,	Adjusted R-squared:  0.0332 
# F-statistic:  6.77 on 1 and 167 DF,  p-value: 0.0101

#hmmm.... that stats should probably be a multilevel model actually (since the relationships are nested by gene):

library(nlme)

BP_Model<-lme(DrugBudger_Log2FC~Pritzker_qPCR_Bipolar, random=~1|Gene,data=DrugBudger_Antipsychotics, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(BP_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: DrugBudger_Antipsychotics 
# AIC     BIC    logLik
# 363.0594 375.579 -177.5297
# 
# Random effects:
#   Formula: ~1 | Gene
# (Intercept)  Residual
# StdDev:   0.2831751 0.6563007
# 
# Fixed effects: DrugBudger_Log2FC ~ Pritzker_qPCR_Bipolar 
#                           Value Std.Error  DF    t-value p-value
# (Intercept)            0.0443440 0.0847087 149  0.5234887  0.6014
# Pritzker_qPCR_Bipolar -0.8189027 0.3747548  18 -2.1851691  0.0423
# Correlation: 
#   (Intr)
# Pritzker_qPCR_Bipolar 0.044 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.22692410 -0.51594251  0.07777773  0.62697398  2.89466938 
# 
# Number of Observations: 169
# Number of Groups: 20 


Schiz_Model<-lme(DrugBudger_Log2FC~Pritzker_qPCR_Schiz, random=~1|Gene,data=DrugBudger_Antipsychotics, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Schiz_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: DrugBudger_Antipsychotics 
# AIC      BIC    logLik
# 362.3151 374.8347 -177.1576
# 
# Random effects:
#   Formula: ~1 | Gene
# (Intercept)  Residual
# StdDev:   0.2839143 0.6546303
# 
# Fixed effects: DrugBudger_Log2FC ~ Pritzker_qPCR_Schiz 
# Value  Std.Error  DF    t-value p-value
# (Intercept)          0.0001706 0.08756261 149  0.0019486  0.9984
# Pritzker_qPCR_Schiz -0.6742297 0.28710518  18 -2.3483717  0.0305
# Correlation: 
#   (Intr)
# Pritzker_qPCR_Schiz 0.255 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.24002370 -0.52692257  0.05793599  0.63348003  2.88081263 
# 
# Number of Observations: 169
# Number of Groups: 20 

library(scales)

pdf("DrugBudger_Antipsychotics_Vs_Pritzker_Antipsychotics.pdf", width=5, height=5.5)
plot(DrugBudger_Antipsychotics$DrugBudger_Log2FC~DrugBudger_Antipsychotics$Pritzker_Antipsychotics, col=alpha(1, 0.5), pch=16, ylab="DrugBudger: Antipsychotic Log2FC", xlab="Pritzker qPCR: Antipsychotic Log2FC")
Temp_Line<-lm(DrugBudger_Antipsychotics$DrugBudger_Log2FC~DrugBudger_Antipsychotics$Pritzker_Antipsychotics)
abline(Temp_Line, col=2, lwd=3)
abline(a=0, b=0)
dev.off()

summary.lm(Temp_Line)
# Call:
#   lm(formula = DrugBudger_Antipsychotics$DrugBudger_Log2FC ~ DrugBudger_Antipsychotics$Pritzker_Antipsychotics)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.3641 -0.4760  0.1169  0.4989  2.5352 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                       -0.00311    0.05850  -0.053    0.958
# DrugBudger_Antipsychotics$Pritzker_Antipsychotics  0.53493    0.34037   1.572    0.118
# 
# Residual standard error: 0.727 on 167 degrees of freedom
# Multiple R-squared:  0.01457,	Adjusted R-squared:  0.008673 
# F-statistic:  2.47 on 1 and 167 DF,  p-value: 0.1179

#Proper stats:
Drug_Model<-lme(DrugBudger_Log2FC~Pritzker_Antipsychotics, random=~1|Gene, data=DrugBudger_Antipsychotics, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Drug_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: DrugBudger_Antipsychotics 
# AIC      BIC    logLik
# 365.9815 378.5011 -178.9908
# 
# Random effects:
#   Formula: ~1 | Gene
# (Intercept)  Residual
# StdDev:   0.3172633 0.6567847
# 
# Fixed effects: DrugBudger_Log2FC ~ Pritzker_Antipsychotics 
#                             Value Std.Error  DF   t-value p-value
# (Intercept)             0.0106103 0.0975166 149 0.1088047  0.9135
# Pritzker_Antipsychotics 0.6745108 0.5395012  18 1.2502489  0.2272
# Correlation: 
#   (Intr)
# Pritzker_Antipsychotics -0.363
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.16980156 -0.54639091  0.07579826  0.57167525  2.81630455 
# 
# Number of Observations: 169
# Number of Groups: 20 

#Positively correlated but not significant.
#Hmm... I'm still not sure how much I trust these results.

pdf("DrugBudger_Corticosteroids_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(DrugBudger_Corticosteroids$DrugBudger_Log2FC~DrugBudger_Corticosteroids$Pritzker_qPCR_Bipolar, col=3, pch=16, ylab="DrugBudger: Corticosteroid Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(DrugBudger_Corticosteroids$DrugBudger_Log2FC~DrugBudger_Corticosteroids$Pritzker_qPCR_Bipolar)
abline(BP_Line, col=3, lwd=3)

points(DrugBudger_Corticosteroids$DrugBudger_Log2FC~DrugBudger_Corticosteroids$Pritzker_qPCR_Schiz, col=4, pch=16)
Schiz_Line<-lm(DrugBudger_Corticosteroids$DrugBudger_Log2FC~DrugBudger_Corticosteroids$Pritzker_qPCR_Schiz)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

BP_Model<-lme(DrugBudger_Log2FC~Pritzker_qPCR_Bipolar, random=~1|Gene,data=DrugBudger_Corticosteroids, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(BP_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: DrugBudger_Corticosteroids 
# AIC      BIC    logLik
# 184.0003 193.2166 -88.00015
# 
# Random effects:
#   Formula: ~1 | Gene
# (Intercept) Residual
# StdDev:   0.1107385 0.787322
# 
# Fixed effects: DrugBudger_Log2FC ~ Pritzker_qPCR_Bipolar 
# Value Std.Error DF    t-value p-value
# (Intercept)           -0.0138369 0.0978868 57 -0.1413561  0.8881
# Pritzker_qPCR_Bipolar  1.2443415 0.4646108 15  2.6782448  0.0172
# Correlation: 
#   (Intr)
# Pritzker_qPCR_Bipolar 0.018 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.10711338 -0.72225640  0.09129824  0.74566194  2.51486849 
# 
# Number of Observations: 74
# Number of Groups: 17 

Schiz_Model<-lme(DrugBudger_Log2FC~Pritzker_qPCR_Schiz, random=~1|Gene,data=DrugBudger_Corticosteroids, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Schiz_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: DrugBudger_Corticosteroids 
# AIC      BIC    logLik
# 183.4893 192.7056 -87.74467
# 
# Random effects:
#   Formula: ~1 | Gene
# (Intercept)  Residual
# StdDev:   0.1028116 0.7855511
# 
# Fixed effects: DrugBudger_Log2FC ~ Pritzker_qPCR_Schiz 
# Value Std.Error DF   t-value p-value
# (Intercept)         0.0666029 0.1017992 57 0.6542574  0.5156
# Pritzker_qPCR_Schiz 1.0064662 0.3601976 15 2.7942054  0.0136
# Correlation: 
#   (Intr)
# Pritzker_qPCR_Schiz 0.304 
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.0662591 -0.6917951  0.2400279  0.6445472  2.5298479 
# 
# Number of Observations: 74
# Number of Groups: 17 


#Making other parallel plots for the Pritzker Exploratory Variables:

ExploratoryDrugResults_TopGenes_ForPlotting<-read.csv("Pritzker_ExploratoryDrugResults_TopGenes_ForPlotting_R.csv", header=TRUE, stringsAsFactors = FALSE)
colnames(ExploratoryDrugResults_TopGenes_ForPlotting)

# [1] "Gene.Symbol"                       "Pritzker_qPCR_Bipolar"             "Pritzker_qPCR_Schiz"              
# [4] "Microarray_MetaAnalysis_Schiz"     "Pritzker_qPCR_Alcohol"             "Gandal_Microarray_Meta_AAD"       
# [7] "Pritzker_qPCR_Cannabinoid"         "Pritzker_qPCR_Tobacco"             "Pritzker_qPCR_Stimulants"         
# [10] "Pritzker_qPCR_Opioids"             "Pritzker_qPCR_Overdose"            "Pritzker_Antipsychotics"          
# [13] "Iwamoto_Microarray_GreyMatter_MDD" "Gandal_Microarray_Meta_MDD"        "Pritzker.qPCR_Antidepressants"    
# [16] "Pritzker_qPCR_Suicide"  


pdf("Gandal_Alcohol_Vs_Pritzker_Alcohol.pdf", width=5, height=5.5)
plot(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Alcohol, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=alpha(1, 0.5), pch=16, ylab="Gandal et al.: AAD Log2FC", xlab="Pritzker qPCR: Alcohol Log2FC")
Temp_Line<-lm(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Alcohol, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Temp_Line, col=1, lwd=3)
dev.off()

summary.lm(Temp_Line)

# Call:
#   lm(formula = Gandal_Microarray_Meta_AAD ~ Pritzker_qPCR_Alcohol, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32333 -0.14713 -0.01891  0.06039  0.65675 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.03581    0.05314   0.674    0.509
# Pritzker_qPCR_Alcohol  0.06198    0.34492   0.180    0.860
# 
# Residual standard error: 0.2283 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.001896,	Adjusted R-squared:  -0.05682 
# F-statistic: 0.03229 on 1 and 17 DF,  p-value: 0.8595

sqrt(0.001896)
#0.04354308


pdf("Gandal_Alcohol_Vs_Pritzker_Cannabinoid.pdf", width=5, height=5.5)
plot(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Cannabinoid, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=alpha(1, 0.5), pch=16, ylab="Gandal et al.: AAD Log2FC", xlab="Pritzker qPCR: Cannabinoid Log2FC")
Temp_Line<-lm(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Cannabinoid, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Temp_Line, col=1, lwd=3)
dev.off()

summary.lm(Temp_Line)
# Call:
#   lm(formula = Gandal_Microarray_Meta_AAD ~ Pritzker_qPCR_Cannabinoid, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.25497 -0.11997 -0.01425  0.07075  0.47525 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                0.04522    0.04474   1.011   0.3263  
# Pritzker_qPCR_Cannabinoid  0.63407    0.24747   2.562   0.0202 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1941 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.2786,	Adjusted R-squared:  0.2362 
# F-statistic: 6.565 on 1 and 17 DF,  p-value: 0.02019

sqrt(0.2786)
#[1] 0.5278257

summary.lm(lm(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Opioids, data=ExploratoryDrugResults_TopGenes_ForPlotting))
# Call:
#   lm(formula = Gandal_Microarray_Meta_AAD ~ Pritzker_qPCR_Opioids, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.33856 -0.06203 -0.02782  0.03580  0.52977 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.05780    0.05100   1.133    0.273
# Pritzker_qPCR_Opioids  0.15402    0.09574   1.609    0.126
# 
# Residual standard error: 0.2129 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.1321,	Adjusted R-squared:  0.08108 
# F-statistic: 2.588 on 1 and 17 DF,  p-value: 0.1261

sqrt(0.1321)
#0.3634556

summary.lm(lm(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Tobacco, data=ExploratoryDrugResults_TopGenes_ForPlotting))

# Call:
#   lm(formula = Gandal_Microarray_Meta_AAD ~ Pritzker_qPCR_Tobacco, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.29559 -0.07437 -0.01934  0.06629  0.47567 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            0.03684    0.03812   0.967  0.34733   
# Pritzker_qPCR_Tobacco  0.62509    0.16045   3.896  0.00116 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1661 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.4717,	Adjusted R-squared:  0.4406 
# F-statistic: 15.18 on 1 and 17 DF,  p-value: 0.001162

sqrt(0.4717)
#[1] 0.6868042

summary.lm(lm(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Stimulants, data=ExploratoryDrugResults_TopGenes_ForPlotting))
# Call:
#   lm(formula = Gandal_Microarray_Meta_AAD ~ Pritzker_qPCR_Stimulants, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.30502 -0.09693 -0.02258  0.05006  0.56464 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)              0.009736   0.044426   0.219   0.8291  
# Pritzker_qPCR_Stimulants 0.528419   0.190937   2.768   0.0132 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1898 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3106,	Adjusted R-squared:   0.27 
# F-statistic: 7.659 on 1 and 17 DF,  p-value: 0.01318

sqrt(0.3106)
#[1] 0.557315

summary.lm(lm(Gandal_Microarray_Meta_AAD~Pritzker_qPCR_Overdose, data=ExploratoryDrugResults_TopGenes_ForPlotting))

# Call:
#   lm(formula = Gandal_Microarray_Meta_AAD ~ Pritzker_qPCR_Overdose, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.33191 -0.06892 -0.02677  0.06900  0.46701 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)             0.06369    0.04964   1.283   0.2167  
# Pritzker_qPCR_Overdose  0.31120    0.15819   1.967   0.0657 .
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2063 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.1854,	Adjusted R-squared:  0.1375 
# F-statistic:  3.87 on 1 and 17 DF,  p-value: 0.06569

sqrt(0.1854)
#[1] 0.430581

pdf("PritzkerAntipsychotics_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Pritzker_Antipsychotics~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=3, pch=16, ylab="Pritzker qPCR: Antipsychotic Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Pritzker_Antipsychotics~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(BP_Line, col=3, lwd=3)

points(Pritzker_Antipsychotics~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=4, pch=16)
Schiz_Line<-lm(Pritzker_Antipsychotics~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.23225 -0.12361 -0.02563  0.06570  0.30226 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.06892    0.03516   1.960   0.0657 .
# Pritzker_qPCR_Bipolar -0.32843    0.14756  -2.226   0.0390 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1568 on 18 degrees of freedom
# Multiple R-squared:  0.2158,	Adjusted R-squared:  0.1723 
# F-statistic: 4.954 on 1 and 18 DF,  p-value: 0.03905

summary.lm(Schiz_Line)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.20584 -0.13068 -0.03384  0.08385  0.32407 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)          0.05963    0.03855   1.547    0.139
# Pritzker_qPCR_Schiz -0.18944    0.12464  -1.520    0.146
# 
# Residual standard error: 0.1667 on 18 degrees of freedom
# Multiple R-squared:  0.1137,	Adjusted R-squared:  0.0645 
# F-statistic:  2.31 on 1 and 18 DF,  p-value: 0.1459

colnames(ExploratoryDrugResults_TopGenes_ForPlotting)

pdf("PritzkerOverdose_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Pritzker_qPCR_Overdose~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=3, pch=16, ylab="Pritzker qPCR: Overdose Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Pritzker_qPCR_Overdose~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(BP_Line, col=3, lwd=3)

points(Pritzker_qPCR_Overdose~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=4, pch=16)
Schiz_Line<-lm(Pritzker_qPCR_Overdose~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Overdose ~ Pritzker_qPCR_Bipolar, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.23611 -0.11486 -0.01392  0.07612  0.41918 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.08793    0.03873  -2.270   0.0357 *  
#   Pritzker_qPCR_Bipolar  1.03364    0.16255   6.359 5.45e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1728 on 18 degrees of freedom
# Multiple R-squared:  0.692,	Adjusted R-squared:  0.6748 
# F-statistic: 40.43 on 1 and 18 DF,  p-value: 5.451e-06

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Overdose ~ Pritzker_qPCR_Schiz, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.21162 -0.10524 -0.01067  0.05814  0.39971 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -0.03971    0.03804  -1.044     0.31    
# Pritzker_qPCR_Schiz  0.83807    0.12302   6.813 2.23e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1646 on 18 degrees of freedom
# Multiple R-squared:  0.7206,	Adjusted R-squared:  0.705 
#F-statistic: 46.41 on 1 and 18 DF,  p-value: 2.229e-06

pdf("PritzkerOpioids_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Pritzker_qPCR_Opioids~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=3, pch=16, ylab="Pritzker qPCR: Opioids Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Pritzker_qPCR_Opioids~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(BP_Line, col=3, lwd=3)

points(Pritzker_qPCR_Opioids~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=4, pch=16)
Schiz_Line<-lm(Pritzker_qPCR_Opioids~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Opioids ~ Pritzker_qPCR_Bipolar, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.31416 -0.16072 -0.03656  0.08528  0.45339 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.12831    0.05059  -2.536   0.0207 *  
#   Pritzker_qPCR_Bipolar  1.89343    0.21230   8.919 5.04e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2257 on 18 degrees of freedom
# Multiple R-squared:  0.8155,	Adjusted R-squared:  0.8052 
#F-statistic: 79.54 on 1 and 18 DF,  p-value: 5.042e-08

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Opioids ~ Pritzker_qPCR_Schiz, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.22296 -0.14038 -0.06905  0.06714  0.54125 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -0.04125    0.04985  -0.827    0.419    
# Pritzker_qPCR_Schiz  1.51913    0.16119   9.424 2.21e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2156 on 18 degrees of freedom
# Multiple R-squared:  0.8315,	Adjusted R-squared:  0.8221 
# F-statistic: 88.82 on 1 and 18 DF,  p-value: 2.207e-08


colnames(ExploratoryDrugResults_TopGenes_ForPlotting)

BP_Line<-lm(Pritzker_qPCR_Alcohol~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
Schiz_Line<-lm(Pritzker_qPCR_Alcohol~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)

summary.lm(BP_Line)
# Call:
#   lm(formula = Pritzker_qPCR_Alcohol ~ Pritzker_qPCR_Bipolar, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.27712 -0.08529 -0.00393  0.07095  0.44492 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)           -0.02603    0.03469  -0.750    0.463
# Pritzker_qPCR_Bipolar  0.08661    0.14556   0.595    0.559
# 
# Residual standard error: 0.1547 on 18 degrees of freedom
# Multiple R-squared:  0.01929,	Adjusted R-squared:  -0.0352 
# F-statistic: 0.354 on 1 and 18 DF,  p-value: 0.5593

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Alcohol ~ Pritzker_qPCR_Schiz, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.28672 -0.07754  0.00091  0.07099  0.42383 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)         -0.02434    0.03600  -0.676    0.508
# Pritzker_qPCR_Schiz  0.04030    0.11640   0.346    0.733
# 
# Residual standard error: 0.1557 on 18 degrees of freedom
# Multiple R-squared:  0.006616,	Adjusted R-squared:  -0.04857 
# F-statistic: 0.1199 on 1 and 18 DF,  p-value: 0.7332

BP_Line<-lm(Pritzker_qPCR_Cannabinoid~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
Schiz_Line<-lm(Pritzker_qPCR_Cannabinoid~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)

summary.lm(BP_Line)
# Call:
#   lm(formula = Pritzker_qPCR_Cannabinoid ~ Pritzker_qPCR_Bipolar, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.22473 -0.06678  0.02567  0.05269  0.36154 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.01995    0.03111  -0.641 0.529401    
# Pritzker_qPCR_Bipolar  0.53233    0.13055   4.078 0.000707 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1388 on 18 degrees of freedom
# Multiple R-squared:  0.4802,	Adjusted R-squared:  0.4513 
# F-statistic: 16.63 on 1 and 18 DF,  p-value: 0.0007066

summary.lm(Schiz_Line)
# Call:
#   lm(formula = Pritzker_qPCR_Cannabinoid ~ Pritzker_qPCR_Schiz, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.224873 -0.057074 -0.008483  0.087783  0.272770 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.006765   0.029610   0.228 0.821863    
# Pritzker_qPCR_Schiz 0.455602   0.095743   4.759 0.000157 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1281 on 18 degrees of freedom
# Multiple R-squared:  0.5571,	Adjusted R-squared:  0.5325 
# F-statistic: 22.64 on 1 and 18 DF,  p-value: 0.000157

BP_Line<-lm(Pritzker_qPCR_Tobacco~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
Schiz_Line<-lm(Pritzker_qPCR_Tobacco~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)

summary.lm(BP_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Tobacco ~ Pritzker_qPCR_Bipolar, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.25629 -0.11966 -0.06113  0.08577  0.65739 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           -0.008288   0.047479  -0.175   0.8634  
# Pritzker_qPCR_Bipolar  0.541900   0.199252   2.720   0.0141 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2118 on 18 degrees of freedom
# Multiple R-squared:  0.2912,	Adjusted R-squared:  0.2519 
# F-statistic: 7.397 on 1 and 18 DF,  p-value: 0.01405

sqrt(0.2912)
#0.5396295

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Tobacco ~ Pritzker_qPCR_Schiz, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.25699 -0.10150 -0.04837  0.03760  0.53385 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          0.02393    0.04361   0.549  0.58998   
# Pritzker_qPCR_Schiz  0.52777    0.14103   3.742  0.00149 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1887 on 18 degrees of freedom
# Multiple R-squared:  0.4376,	Adjusted R-squared:  0.4063 
# F-statistic: 14.01 on 1 and 18 DF,  p-value: 0.001491


BP_Line<-lm(Pritzker_qPCR_Stimulants~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
Schiz_Line<-lm(Pritzker_qPCR_Stimulants~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)

summary.lm(BP_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Stimulants ~ Pritzker_qPCR_Bipolar, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.26704 -0.11768 -0.01593  0.07443  0.54765 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            0.04354    0.04038   1.078  0.29522   
# Pritzker_qPCR_Bipolar  0.64926    0.16948   3.831  0.00122 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1801 on 18 degrees of freedom
# Multiple R-squared:  0.4491,	Adjusted R-squared:  0.4185 
# F-statistic: 14.68 on 1 and 18 DF,  p-value: 0.001224

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Stimulants ~ Pritzker_qPCR_Schiz, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.26734 -0.10864  0.00942  0.07369  0.43198 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.07724    0.03771   2.048 0.055415 .  
# Pritzker_qPCR_Schiz  0.56996    0.12194   4.674 0.000189 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1631 on 18 degrees of freedom
# Multiple R-squared:  0.5483,	Adjusted R-squared:  0.5232 
# F-statistic: 21.85 on 1 and 18 DF,  p-value: 0.0001889

##########


pdf("IwamotoMDD_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Iwamoto_Microarray_GreyMatter_MDD~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=3, pch=16, ylab="Iwamoto Microarray: MDD Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Iwamoto_Microarray_GreyMatter_MDD~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(BP_Line, col=3, lwd=3)

points(Iwamoto_Microarray_GreyMatter_MDD~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=4, pch=16)
Schiz_Line<-lm(Iwamoto_Microarray_GreyMatter_MDD~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Iwamoto_Microarray_GreyMatter_MDD ~ Pritzker_qPCR_Bipolar, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.33884 -0.04558  0.01236  0.06910  0.14001 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           0.008234   0.025256   0.326  0.74838   
# Pritzker_qPCR_Bipolar 0.352898   0.104489   3.377  0.00358 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1094 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.4015,	Adjusted R-squared:  0.3663 
# F-statistic: 11.41 on 1 and 17 DF,  p-value: 0.003579

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Iwamoto_Microarray_GreyMatter_MDD ~ Pritzker_qPCR_Schiz, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32484 -0.03499  0.01946  0.05252  0.11747 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          0.02520    0.02482   1.015   0.3242   
# Pritzker_qPCR_Schiz  0.29864    0.07841   3.809   0.0014 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1039 on 17 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.4604,	Adjusted R-squared:  0.4287 
# F-statistic:  14.5 on 1 and 17 DF,  p-value: 0.001405


#This may be partially redundant - if I remember correctly, Gandal included Iwamoto's data in their microarray meta-analysis
#Nope - double-checked, all independent datasets, although from a variety of frontal cortical regions (BA9, BA46, BA25 (subgenual))

pdf("GandalMDD_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Gandal_Microarray_Meta_MDD~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=3, pch=16, ylab="Gandal Microarray Meta-Analysis: MDD Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Gandal_Microarray_Meta_MDD~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(BP_Line, col=3, lwd=3)

points(Gandal_Microarray_Meta_MDD~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=4, pch=16)
Schiz_Line<-lm(Gandal_Microarray_Meta_MDD~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Gandal_Microarray_Meta_MDD ~ Pritzker_qPCR_Bipolar, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.19623 -0.01783  0.01017  0.03165  0.11922 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)           -0.02386    0.01555  -1.535    0.142
# Pritzker_qPCR_Bipolar  0.09650    0.06525   1.479    0.156
# 
# Residual standard error: 0.06935 on 18 degrees of freedom
# Multiple R-squared:  0.1084,	Adjusted R-squared:  0.05883 
# F-statistic: 2.188 on 1 and 18 DF,  p-value: 0.1564

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Gandal_Microarray_Meta_MDD ~ Pritzker_qPCR_Schiz, 
#      data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.189966 -0.018999  0.004755  0.038724  0.099449 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         -0.01858    0.01571  -1.182   0.2526  
# Pritzker_qPCR_Schiz  0.08821    0.05081   1.736   0.0996 .
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06797 on 18 degrees of freedom
# Multiple R-squared:  0.1434,	Adjusted R-squared:  0.09583 
# F-statistic: 3.014 on 1 and 18 DF,  p-value: 0.09965


pdf("PritzkerSuicide_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Pritzker_qPCR_Suicide~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=3, pch=16, ylab="Pritzker qPCR: Suicide Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Pritzker_qPCR_Suicide~Pritzker_qPCR_Bipolar, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(BP_Line, col=3, lwd=3)

points(Pritzker_qPCR_Suicide~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting, col=4, pch=16)
Schiz_Line<-lm(Pritzker_qPCR_Suicide~Pritzker_qPCR_Schiz, data=ExploratoryDrugResults_TopGenes_ForPlotting)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Suicide ~ Pritzker_qPCR_Bipolar, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.40087 -0.08278  0.01547  0.09001  0.33297 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.04049    0.03971   1.020   0.3214  
# Pritzker_qPCR_Bipolar -0.47115    0.16666  -2.827   0.0112 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1771 on 18 degrees of freedom
# Multiple R-squared:  0.3075,	Adjusted R-squared:  0.269 
# F-statistic: 7.992 on 1 and 18 DF,  p-value: 0.01117

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Pritzker_qPCR_Suicide ~ Pritzker_qPCR_Schiz, data = ExploratoryDrugResults_TopGenes_ForPlotting)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.44743 -0.07201 -0.02201  0.10215  0.28933 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          0.01725    0.03974   0.434   0.6695   
# Pritzker_qPCR_Schiz -0.39815    0.12850  -3.098   0.0062 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1719 on 18 degrees of freedom
# Multiple R-squared:  0.3478,	Adjusted R-squared:  0.3116 
# F-statistic:   9.6 on 1 and 18 DF,  p-value: 0.006201


#Making a large correlation matrix for the exploratory analyses involving the top 20 genes:

ExploratoryVariablesLog2FC<-read.csv("ExploratoryVariableLog2FC_forR.csv", header=TRUE, stringsAsFactors = FALSE)
str(ExploratoryVariablesLog2FC)

write.csv(cor(as.matrix(ExploratoryVariablesLog2FC[,c(2:58)]), use="pairwise.complete.obs"), "CorrelationMatrix_Log2FC_AllExploratoryVar_Diagnosis.csv")

pdf("Heatmap_CorrMatrix_Log2FC_AllExploratoryVar_Diagnosis.pdf", width=10, height=10)
heatmap(cor(as.matrix(ExploratoryVariablesLog2FC[,c(2:58)]), use="pairwise.complete.obs"))
dev.off()



#More plots from the exploratory analysis:

colnames(ExploratoryVariablesLog2FC)

pdf("PritzkerDisorganized_Speech_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Disorganized_Speech~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=3, pch=16, ylab="Pritzker qPCR: Disorganized Speech Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Disorganized_Speech~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(BP_Line, col=3, lwd=3)

points(Disorganized_Speech~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=4, pch=16)
Schiz_Line<-lm(Disorganized_Speech~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Disorganized_Speech ~ BP_Pritzker_qPCR_GreyMatter, 
#      data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38931 -0.07317 -0.02352  0.07995  0.38811 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 -0.01609    0.03768  -0.427    0.675    
# BP_Pritzker_qPCR_GreyMatter  0.95476    0.15813   6.038 1.04e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1681 on 18 degrees of freedom
# Multiple R-squared:  0.6695,	Adjusted R-squared:  0.6511 
# F-statistic: 36.46 on 1 and 18 DF,  p-value: 1.043e-05

summary.lm(Schiz_Line)
# Call:
#   lm(formula = Disorganized_Speech ~ SCHIZ_Pritzker_qPCR_GreyMatter, 
#      data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.253721 -0.092876 -0.004094  0.112088  0.242764 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     0.02973    0.03533   0.841    0.411    
# SCHIZ_Pritzker_qPCR_GreyMatter  0.79039    0.11423   6.919 1.81e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1528 on 18 degrees of freedom
# Multiple R-squared:  0.7267,	Adjusted R-squared:  0.7116 
# F-statistic: 47.87 on 1 and 18 DF,  p-value: 1.815e-06

pdf("Fatigue_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Fatigue~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=3, pch=16, ylab="Pritzker qPCR: Fatigue Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Fatigue~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(BP_Line, col=3, lwd=3)

points(Fatigue~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=4, pch=16)
Schiz_Line<-lm(Fatigue~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Fatigue ~ BP_Pritzker_qPCR_GreyMatter, data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32009 -0.05525 -0.02913  0.05580  0.30856 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                 -0.03332    0.02983  -1.117  0.27858   
# BP_Pritzker_qPCR_GreyMatter  0.47127    0.12517   3.765  0.00142 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.133 on 18 degrees of freedom
# Multiple R-squared:  0.4406,	Adjusted R-squared:  0.4095 
# F-statistic: 14.17 on 1 and 18 DF,  p-value: 0.001418

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Fatigue ~ SCHIZ_Pritzker_qPCR_GreyMatter, data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32124 -0.04325 -0.00711  0.03699  0.37606 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                    -0.01128    0.03020  -0.374    0.713   
# SCHIZ_Pritzker_qPCR_GreyMatter  0.38286    0.09766   3.920    0.001 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1306 on 18 degrees of freedom
# Multiple R-squared:  0.4606,	Adjusted R-squared:  0.4306 
# F-statistic: 15.37 on 1 and 18 DF,  p-value: 0.001003

pdf("WorthOrHelpless_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Worth_or_Helpless~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=3, pch=16, ylab="Pritzker qPCR: Decreased Worth or Helpless Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Worth_or_Helpless~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(BP_Line, col=3, lwd=3)

points(Worth_or_Helpless~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=4, pch=16)
Schiz_Line<-lm(Worth_or_Helpless~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Worth_or_Helpless ~ BP_Pritzker_qPCR_GreyMatter, 
#      data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.311112 -0.095395  0.002545  0.086386  0.250261 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 -0.06247    0.03142  -1.988 0.062234 .  
# BP_Pritzker_qPCR_GreyMatter  0.57739    0.13187   4.379 0.000362 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1402 on 18 degrees of freedom
# Multiple R-squared:  0.5158,	Adjusted R-squared:  0.4889 
# F-statistic: 19.17 on 1 and 18 DF,  p-value: 0.0003622

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Worth_or_Helpless ~ SCHIZ_Pritzker_qPCR_GreyMatter, 
#      data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32064 -0.09479  0.01881  0.08308  0.24319 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                    -0.03978    0.03546  -1.122    0.277   
# SCHIZ_Pritzker_qPCR_GreyMatter  0.41407    0.11465   3.612    0.002 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1534 on 18 degrees of freedom
# Multiple R-squared:  0.4202,	Adjusted R-squared:  0.3879 
# F-statistic: 13.04 on 1 and 18 DF,  p-value: 0.001995

pdf("Weight_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Weight~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=3, pch=16, ylab="Pritzker qPCR: Weight Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Weight~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(BP_Line, col=3, lwd=3)

points(Weight~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=4, pch=16)
Schiz_Line<-lm(Weight~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Weight ~ BP_Pritzker_qPCR_GreyMatter, data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.238004 -0.031254 -0.007428  0.052634  0.219145 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                 -0.01045    0.02201  -0.475  0.64077   
# BP_Pritzker_qPCR_GreyMatter  0.33208    0.09236   3.595  0.00207 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09818 on 18 degrees of freedom
# Multiple R-squared:  0.418,	Adjusted R-squared:  0.3856 
# F-statistic: 12.93 on 1 and 18 DF,  p-value: 0.002069

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Weight ~ SCHIZ_Pritzker_qPCR_GreyMatter, data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.21149 -0.03772  0.01142  0.05452  0.20245 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                    0.004497   0.022793   0.197  0.84579   
# SCHIZ_Pritzker_qPCR_GreyMatter 0.262275   0.073701   3.559  0.00224 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.0986 on 18 degrees of freedom
# Multiple R-squared:  0.413,	Adjusted R-squared:  0.3804 
# F-statistic: 12.66 on 1 and 18 DF,  p-value: 0.002244


pdf("InteractionsWLegalSystem_Vs_PritzkerDiagnosis.pdf", width=5, height=5.5)
plot(Interactions_w_LegalSystem~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=3, pch=16, ylab="Pritzker qPCR: Interactions w/ Legal System Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Interactions_w_LegalSystem~BP_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(BP_Line, col=3, lwd=3)

points(Interactions_w_LegalSystem~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC, col=4, pch=16)
Schiz_Line<-lm(Interactions_w_LegalSystem~SCHIZ_Pritzker_qPCR_GreyMatter, data=ExploratoryVariablesLog2FC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

summary.lm(BP_Line)

# Call:
#   lm(formula = Interactions_w_LegalSystem ~ BP_Pritzker_qPCR_GreyMatter, 
#      data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.17841 -0.06971  0.00521  0.08531  0.17862 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                  0.06292    0.02426   2.594   0.0183 *
#   BP_Pritzker_qPCR_GreyMatter  0.26726    0.10180   2.625   0.0172 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1082 on 18 degrees of freedom
# Multiple R-squared:  0.2769,	Adjusted R-squared:  0.2367 
# F-statistic: 6.892 on 1 and 18 DF,  p-value: 0.01716

summary.lm(Schiz_Line)

# Call:
#   lm(formula = Interactions_w_LegalSystem ~ SCHIZ_Pritzker_qPCR_GreyMatter, 
#      data = ExploratoryVariablesLog2FC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.169748 -0.073269 -0.007074  0.095987  0.165816 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                     0.07811    0.02301   3.394  0.00323 **
#   SCHIZ_Pritzker_qPCR_GreyMatter  0.25142    0.07441   3.379  0.00334 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09954 on 18 degrees of freedom
# Multiple R-squared:  0.3881,	Adjusted R-squared:  0.3541 
# F-statistic: 11.42 on 1 and 18 DF,  p-value: 0.003343


#I just realized that I have the full Martin et al. antipsychotic results that I could plot too. I averaged the results from the two outputs provided from the paper. I have not averaged them by gene yet - there are multiple probe entries for each gene:
#I have not filtered them yet to focus on our top (most reliable) diagnosis effects similar to the other antipsychotic graphs

setwd("~/Documents/Microarray Gen/FrontalPole/Drugs")

PritzkerqPCRvsMartin<-read.csv("FrontalPole_qPCRGenes_vsMartin_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
colnames(PritzkerqPCRvsMartin)
str(PritzkerqPCRvsMartin)


pdf("MartinHaloperidol_Vs_PritzkerDiagnosis_allTargets.pdf", width=5, height=5.5)
plot(Ratio.HAL.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin, col=3, pch=16, ylab="Martin et al. Primates: Haloperidol vs. Vehicle Ratio", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Ratio.HAL.vs..VEH_AVE~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin)
abline(BP_Line, col=3, lwd=3)

points(Ratio.HAL.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin, col=4, pch=16)
Schiz_Line<-lm(Ratio.HAL.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin)
abline(Schiz_Line, col=4, lwd=3)
dev.off()

#using mlm to evaluate results so that the (supposed) covariance between results targeting the same gene is properly modeled

library(nlme)

BP_Model<-lme(Ratio.HAL.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta, random=~1|Gene.Symbol,data=PritzkerqPCRvsMartin, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(BP_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: PritzkerqPCRvsMartin 
# AIC       BIC   logLik
# -620.5382 -608.1878 314.2691
# 
# Random effects:
#   Formula: ~1 | Gene.Symbol
# (Intercept)   Residual
# StdDev:  0.03841325 0.02282189
# 
# Fixed effects: Ratio.HAL.vs..VEH_AVE ~ Diagnosis_BP_PostHocSummary_MLM_Beta 
# Value  Std.Error DF   t-value p-value
# (Intercept)                           1.0349301 0.00489279 84 211.52146  0.0000
# Diagnosis_BP_PostHocSummary_MLM_Beta -0.0460861 0.03395454 76  -1.35729  0.1787
# Correlation: 
#   (Intr)
# Diagnosis_BP_PostHocSummary_MLM_Beta -0.078
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.49815799 -0.42683976 -0.07530992  0.51580898  2.57662135 
# 
# Number of Observations: 162
# Number of Groups: 78

Schiz_Model<-lme(Ratio.HAL.vs..VEH_AVE ~Diagnosis_Schiz_PostHocSummary_MLM_Beta, random=~1|Gene.Symbol,data=PritzkerqPCRvsMartin, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Schiz_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: PritzkerqPCRvsMartin 
# AIC       BIC   logLik
# -622.2842 -609.9338 315.1421
# 
# Random effects:
#   Formula: ~1 | Gene.Symbol
# (Intercept)   Residual
# StdDev:   0.0378123 0.02284805
# 
# Fixed effects: Ratio.HAL.vs..VEH_AVE ~ Diagnosis_Schiz_PostHocSummary_MLM_Beta 
# Value   Std.Error DF   t-value p-value
# (Intercept)                              1.0335670 0.004837369 84 213.66307  0.0000
# Diagnosis_Schiz_PostHocSummary_MLM_Beta -0.0535837 0.028088879 76  -1.90765  0.0602
# Correlation: 
#   (Intr)
# Diagnosis_Schiz_PostHocSummary_MLM_Beta 0.091 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.50729916 -0.41478535 -0.05654293  0.50579613  2.57126509 
# 
# Number of Observations: 162
# Number of Groups: 78 


pdf("MartinOlanzapine_Vs_PritzkerDiagnosis_allTargets.pdf", width=5, height=5.5)
plot(Ratio.OLANZ.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin, col=3, pch=16, ylab="Martin et al. Primates: Olanzapine vs. Vehicle Ratio", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Ratio.OLANZ.vs..VEH_AVE~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin)
abline(BP_Line, col=3, lwd=3)

points(Ratio.OLANZ.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin, col=4, pch=16)
Schiz_Line<-lm(Ratio.OLANZ.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin)
abline(Schiz_Line, col=4, lwd=3)
dev.off()


BP_Model<-lme(Ratio.OLANZ.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta, random=~1|Gene.Symbol,data=PritzkerqPCRvsMartin, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(BP_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: PritzkerqPCRvsMartin 
# AIC       BIC   logLik
# -659.1128 -646.7624 333.5564
# 
# Random effects:
#   Formula: ~1 | Gene.Symbol
# (Intercept)   Residual
# StdDev:  0.02488765 0.02394299
# 
# Fixed effects: Ratio.OLANZ.vs..VEH_AVE ~ Diagnosis_BP_PostHocSummary_MLM_Beta 
# Value   Std.Error DF   t-value p-value
# (Intercept)                           1.0110260 0.003602469 84 280.64808  0.0000
# Diagnosis_BP_PostHocSummary_MLM_Beta -0.0275484 0.025452043 76  -1.08236  0.2825
# Correlation: 
#   (Intr)
# Diagnosis_BP_PostHocSummary_MLM_Beta -0.077
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.21073774 -0.54145016 -0.02618467  0.49952339  3.28905664 
# 
# Number of Observations: 162
# Number of Groups: 78 

Schiz_Model<-lme(Ratio.OLANZ.vs..VEH_AVE ~Diagnosis_Schiz_PostHocSummary_MLM_Beta, random=~1|Gene.Symbol,data=PritzkerqPCRvsMartin, na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Schiz_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: PritzkerqPCRvsMartin 
# AIC       BIC   logLik
# -658.3684 -646.0181 333.1842
# 
# Random effects:
#   Formula: ~1 | Gene.Symbol
# (Intercept)   Residual
# StdDev:  0.02510277 0.02393294
# 
# Fixed effects: Ratio.OLANZ.vs..VEH_AVE ~ Diagnosis_Schiz_PostHocSummary_MLM_Beta 
# Value  Std.Error DF   t-value p-value
# (Intercept)                              1.0104799 0.00363085 84 278.30395  0.0000
# Diagnosis_Schiz_PostHocSummary_MLM_Beta -0.0140418 0.02145228 76  -0.65456  0.5147
# Correlation: 
#   (Intr)
# Diagnosis_Schiz_PostHocSummary_MLM_Beta 0.103 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.21286720 -0.54991862 -0.01642241  0.51078669  3.28965557 
# 
# Number of Observations: 162
# Number of Groups: 78 


#pulling out results for just our top 20 (most reliable) diagnosis-related genes:

Top20Genes<-c("HTR2B", "SST", "DRD4", "BDNF", "DRD2", "DRD3", "HOMER1", "SLC6A12", "COMT", "SLC6A11", "GABRB1", "GRM5", "CALB1", "ABAT", "MAPK1", "NSF", "SNCA", "GPHN", "MAOB", "GFAP")
length(Top20Genes)
#[1] 20


pdf("MartinHaloperidol_Vs_PritzkerDiagnosis_Top20Genes.pdf", width=5, height=5.5)
plot(Ratio.HAL.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,], col=3, pch=16, ylab="Martin et al. Primates: Haloperidol vs. Vehicle Ratio", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Ratio.HAL.vs..VEH_AVE~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,])
abline(BP_Line, col=3, lwd=3)

points(Ratio.HAL.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,], col=4, pch=16)
Schiz_Line<-lm(Ratio.HAL.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,])
abline(Schiz_Line, col=4, lwd=3)
dev.off()

#using mlm to evaluate results so that the (supposed) covariance between results targeting the same gene is properly modeled

library(nlme)

BP_Model<-lme(Ratio.HAL.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta, random=~1|Gene.Symbol,data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,], na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(BP_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol %in% Top20Genes,      ] 
# AIC       BIC   logLik
# -148.9069 -142.1513 78.45343
# 
# Random effects:
#   Formula: ~1 | Gene.Symbol
# (Intercept)   Residual
# StdDev:  0.03361962 0.02476711
# 
# Fixed effects: Ratio.HAL.vs..VEH_AVE ~ Diagnosis_BP_PostHocSummary_MLM_Beta 
# Value  Std.Error DF   t-value p-value
# (Intercept)                           1.0545008 0.00974637 23 108.19424  0.0000
# Diagnosis_BP_PostHocSummary_MLM_Beta -0.0154309 0.04018042 15  -0.38404  0.7063
# Correlation: 
#   (Intr)
# Diagnosis_BP_PostHocSummary_MLM_Beta 0.119 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.45970321 -0.61479989  0.06121382  0.56342090  2.31760484 
# 
# Number of Observations: 40
# Number of Groups: 17 

Schiz_Model<-lme(Ratio.HAL.vs..VEH_AVE ~Diagnosis_Schiz_PostHocSummary_MLM_Beta, random=~1|Gene.Symbol,data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,], na.action = na.omit, control = lmeControl(opt = 'optim'), method="ML")
summary(Schiz_Model)

# Linear mixed-effects model fit by maximum likelihood
# Data: PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol %in% Top20Genes,      ] 
# AIC       BIC   logLik
# -148.9681 -142.2126 78.48405
# 
# Random effects:
#   Formula: ~1 | Gene.Symbol
# (Intercept)   Residual
# StdDev:  0.03350634 0.02477819
# 
# Fixed effects: Ratio.HAL.vs..VEH_AVE ~ Diagnosis_Schiz_PostHocSummary_MLM_Beta 
# Value  Std.Error DF   t-value p-value
# (Intercept)                              1.053559 0.01012287 23 104.07700  0.0000
# Diagnosis_Schiz_PostHocSummary_MLM_Beta -0.014413 0.03167062 15  -0.45509  0.6556
# Correlation: 
#   (Intr)
# Diagnosis_Schiz_PostHocSummary_MLM_Beta 0.301 
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.46009909 -0.60284128  0.05692978  0.59066004  2.31551915 
# 
# Number of Observations: 40
# Number of Groups: 17 

#I might have just killed our statistical power. :( sigh.


pdf("MartinOlanzapine_Vs_PritzkerDiagnosis_Top20Genes.pdf", width=5, height=5.5)
plot(Ratio.OLANZ.vs..VEH_AVE ~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,], col=3, pch=16, ylab="Martin et al. Primates: Olanzapine vs. Vehicle Ratio", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Ratio.OLANZ.vs..VEH_AVE~Diagnosis_BP_PostHocSummary_MLM_Beta , data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,])
abline(BP_Line, col=3, lwd=3)

points(Ratio.OLANZ.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,], col=4, pch=16)
Schiz_Line<-lm(Ratio.OLANZ.vs..VEH_AVE~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=PritzkerqPCRvsMartin[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes,])
abline(Schiz_Line, col=4, lwd=3)
dev.off()

#I'm going to average by gene and output these for our top gene table:

temp<-cbind(tapply(X=PritzkerqPCRvsMartin$Ratio.HAL.vs..VEH_AVE[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes], INDEX = PritzkerqPCRvsMartin$Gene.Symbol[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes], function(y) mean(y, na.rm=TRUE)),
tapply(X=PritzkerqPCRvsMartin$Ratio.OLANZ.vs..VEH_AVE[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes], INDEX = PritzkerqPCRvsMartin$Gene.Symbol[PritzkerqPCRvsMartin$Gene.Symbol%in%Top20Genes], function(y) mean(y, na.rm=TRUE)))
write.csv(temp, "Martin_AverageRatio_forTop20Genes.csv")




