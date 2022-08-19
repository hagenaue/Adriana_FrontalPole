#Adding in functional MetaData & comparisons w/ previous DLPFC results:
#Re-done on July 27, 2021 to include the large Gandal Meta-analysis (instead of some of the smaller previous studies)

####################

library(plyr)

###################

setwd("~/Documents/Microarray Gen/FrontalPole/Sample_MetaData")

#Reading in basic annotation for the qPCR genes:

HousekeepingGenes<-read.csv("HouseKeepingGenes_UsedInPCR_FrontalPole.csv", header=TRUE, stringsAsFactors = FALSE)
str(HousekeepingGenes)
# 'data.frame':	25 obs. of  1 variable:
#   $ GeneSymbol: chr  "18S" "ACTB" "B2M" "CASC3" ...


#Reading in Adriana's custom gene annotation (i.e., the associated functional categories that made these targets interesting to her)

#Note: Annotation for ITPR1 was missing - I added it in (part of GRM1 complex)
GeneCategories<-read.csv("GeneCategories_BYADRIANA.csv", header=T, stringsAsFactors = FALSE)
str(GeneCategories)
# 'data.frame':	115 obs. of  3 variables:
#   $ GeneSymbol                : chr  "GFAP" "AQP4" "GJA1" "S100B" ...
# $ FunctionalCategory_Main   : chr  "ASTROCYTE" "ASTROCYTE" "ASTROCYTE" "ASTROCYTE" ...
# $ FunctionalCategory_Adriana: chr  "ASTROCYTE MARKER" "ASTROCYTE MARKER" "ASTROCYTE MARKER" "ASTROCYTE MARKER" ...

table(GeneCategories$FunctionalCategory_Main)
# ASTROCYTE    DOPAMINE        GABA   GLUTAMATE   MONOAMINE       OTHER PRESYNAPTIC   SEROTONIN 
# 4           5          30          44           9           6           2          15 

table(GeneCategories$FunctionalCategory_Adriana)

# ADENOSINE RECEPTOR                                                    ASTROCYTE MARKER 
# 2                                                                   4 
# BETA-1 ADRENERGIC RECEPTOR                                                    CALCIUM CHANNELS 
# 2                                                                   2 
# CORTICAL INTERNEURON MARKER                                                  DOPAMINE RECEPTORS 
# 3                                                                   5 
# DOPAMINE TRANSPORTER                          ELEMENTS ON THE GABA POSTSYNAPTIC TERMINAL 
# 1                                                                   1 
# ELEMENTS ON THE GLUTAMATERGIC POSTSYNAPTIC DENSITY (PSD)                  ELEMENTS ON THE GLUTAMATERGIC PRESYNAPTIC TERMINAL 
# 5                                                                   4 
# ELEMENTS ON THE PRESYNAPTIC GABAERGIC TERMINAL AND GABA TRASPORTERS                                  ENZYMES RELATED TO GABA METABOLISM 
# 5                                                                   3 
# ENZYMES RELATED TO GLUTAMATE METABOLISM         ENZYMES RELATED TO METABOTROPIC GLUTAMATE RECEPTOR ACTIVITY 
# 3                                                                   5 
# ENZYMES RELATED TO MONOAMINE METABOLISM                                           GABA IONOTROPIC RECEPTORS 
# 6                                                                  15 
# GABA METABOTROPIC RECEPORS          GLUTAMATE TRANSPORTERS, GLUTAMINE AND CYSTEIN TRANSPORTERS 
# 2                                                                   6 
# IONOROPIC GLUTAMATE RECEPTORS NMDA ACTIVATOR                         IONOROPIC GLUTAMATE RECEPTORS NMDA SUBUNITS 
# 1                                                                   4 
# IONOTROPIC GLUTAMATE RECEPORS KAINATE SUBUNITS                        IONOTROPIC GLUTAMATE RECEPTORS AMPA SUBUNITS 
# 4                                                                   4 
# METABOTROPIC GLUAMATE RECEPTORS                                               MONOAMINE TRANSPORTER 
# 8                                                                   1 
# OTHER                                 PRESYNAPTIC MONOAMINERGIC ELEMENTS  
# 1                                                                   2 
# PURINOCEPTOR FOR ATP                                                SEROTONIN RECEPTORS  
# 1                                                                  14 
# SEROTONIN TRANSPORTER 
# 1 


##############################

#Comparing with my other databases of previous knowledge from the DLPFC:

setwd("~/Documents/Microarray Gen/FrontalPole/CompareWDLPFC")

list.files(".")
# [1] "Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_MicroarrayMetaAnalysis.csv"
# [2] "Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_RNASeqMetaAnalysis.csv"    
# [3] "Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_RNASeqMetaAnalysis.xlsx"   
# [4] "Lanz_GeneByCellTypeSubjVar2DF.csv"                                                
# [5] "LanzDLPFC_LM4Basic_Psych.csv"                                                     
# [6] "Maycox_41380_2009_BFmp200918_MOESM434_ESM.xls"                                    
# [7] "Maycox_SchizMicroarray.pdf"                                                       
# [8] "S01_Table_CellTypeSpecificGenes.csv"                                              
# [9] "S07_CellTypeSpecificEffects_ChoiMetaAnalysis_BP.csv"                              
# [10] "S07_CellTypeSpecificEffects_MistryMetaAnalysis_Schiz.csv"                         
# [11] "S07_CellTypeSpecificEffects_ParticularCellTypes.csv"                              
# [12] "S08_Table_Pritzker_PsychResultsDiffModels_M2UsualSuspects.csv"                    
# [13] "S08_Table_Pritzker_PsychResultsDiffModels_M4UsualSuspectsAndCellType.csv"         
# [14] "S09_Table_CMC_PsychEffectsinDifferentModels_BP_M2UsualSuspects.csv"               
# [15] "S09_Table_CMC_PsychEffectsinDifferentModels_Schiz_M2UsualSuspects.csv"            
# [16] "S10_Table_Barnes_PsychEffectsinDifferentModels_Schiz_M2UsualSuspects.csv"         
# [17] "S11_Table_Lanz_PsychEffectsinDifferentModels_BP_M2_UsualSuspects.csv"             
# [18] "S11_Table_Lanz_PsychEffectsinDifferentModels_Schiz_M2_UsualSuspects.csv"          
# [19] "S12_Table_Narayan_PsychEffectsinDifferentModels_Schiz_M2UsualSuspects.csv"     

#Note: 

#Gandal et al. 2018 is a really excellent, highly-powered set of meta-analyses of either cortical microarray data or DLPFC RNA-Seq data - top priority for comparison.
   # A note of caution: The Gandal RNA-Seq meta-analysis is independent of our datasets, but the cortical microarray meta-analysis actually incorporates the Iwamoto et al. 2004 and Maycox et al. 2009 datasets as part of it. A comparison may still be interesting, but we shouldn't emphasize it in our paper.

#Although the Gandal results should be better than the intermediate CMC RNA-Seq findings, the Gandal findings do not include an estimate of the intercept and/or average expression levels, so we will still need to use the CMC RNA-Seq results too if we want to get an estimate of sensitivity for highly vs. lowly expressed genes.

#We should include Lanz and Narayan datasets for a grey-matter only comparison.
    # A note of caution: Lanz was derived from PITT brain bank, and therefore may have subjects/tissue samples that overlap with the Gandal RNA-Seq meta-analysis.
    # A note of caution: Lanz and Narayan were also included in the Gandal microarray meta-analysis.

#We will also need the Pritzker results to have estimates of the effects of other variables (pH, PMI, age, gender)
    # Due to the batch confound, the Pritzker SCHIZ/BP results should not be used as comparison.
    # The Pritzker MDD results are reasonably well-powered and may also be an interesting point of comparison.

#I also have the BrainInABlender database for information about cell type specificity

#...and estimates of cell-type specificity from Pritzker DLPFC that may or may not be useful.

#####################

#Reading in Gandal et al. 2018:

GandalCortex_MicroarrayMeta<-read.csv("Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_MicroarrayMetaAnalysis.csv", header=T, stringsAsFactors = FALSE)
GandalFrontalCortex_RNASeqMeta<-read.csv("Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_RNASeqMetaAnalysis.csv", header=T, stringsAsFactors = FALSE)


str(GandalCortex_MicroarrayMeta)

GandalCortex_MicroarrayMeta_SchizEffects<-GandalCortex_MicroarrayMeta[,c(1:7,11:13)]
GandalCortex_MicroarrayMeta_BPEffects<-GandalCortex_MicroarrayMeta[,c(1:7,14:16)]
GandalCortex_MicroarrayMeta_MDDEffects<-GandalCortex_MicroarrayMeta[,c(1:7,17:19)]
str(GandalCortex_MicroarrayMeta_SchizEffects)
str(GandalCortex_MicroarrayMeta_BPEffects)
str(GandalCortex_MicroarrayMeta_MDDEffects)

str(GandalFrontalCortex_RNASeqMeta)
#'data.frame':	15823 obs. of  9 variables:
#That's annoying - they don't have anything except ensembl annotation. 
#I wonder how many of the ensembl genes from the RNASeqMeta were in the microarray meta (which does have other annotation):
sum(GandalFrontalCortex_RNASeqMeta$X%in%GandalCortex_MicroarrayMeta$ensembl_gene_id)
#[1] 14785
#Most of them, but not all. 

length(GandalFrontalCortex_RNASeqMeta$X)
#[1] 15823

library(org.Hs.eg.db)

x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

xx[1]
# $`1`
# [1] "A1BG"
names(xx[1])
#[1] "1"
names(xx[800])
#[1] "100128946"
#Hmm...
#Apparently this package only maps to Entrez, not Ensembl. :(

columns(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT" 

#Ah - I think I need to map both Symbol and Ensembl back to Entrez to align them.
#Or maybe not?

keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"     

uniKeys <- head(keys(org.Hs.eg.db, keytype="ENSEMBL"))
cols <- c("SYMBOL", "ENTREZID")
select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")
# 'select()' returned 1:1 mapping between keys and columns
# ENSEMBL   SYMBOL ENTREZID
# 1 ENSG00000121410     A1BG        1
# 2 ENSG00000175899      A2M        2
# 3 ENSG00000256069    A2MP1        3
# 4 ENSG00000171428     NAT1        9
# 5 ENSG00000156006     NAT2       10
# 6 ENSG00000196136 SERPINA3       12

#Lovely.

colnames(GandalFrontalCortex_RNASeqMeta)[1]<-"ENSEMBL"
str(GandalFrontalCortex_RNASeqMeta)
#'data.frame':	15823 obs. of  9 variables:

sum(is.na(GandalFrontalCortex_RNASeqMeta$ENSEMBL))
#[1] 0
sum(duplicated(GandalFrontalCortex_RNASeqMeta$ENSEMBL))
#[1] 0

EnsemblVsGeneSymbol_Human<-select(org.Hs.eg.db, keys=GandalFrontalCortex_RNASeqMeta$ENSEMBL, columns=cols, keytype="ENSEMBL")
#'select()' returned 1:many mapping between keys and columns
#'#Note: Ensembl gene annotation is mapping to either multiple Entrez genes, multiple gene symbols, or both.

str(EnsemblVsGeneSymbol_Human)
# 'data.frame':	16188 obs. of  3 variables:
#  $ ENSEMBL : chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
# $ SYMBOL  : chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...
# $ ENTREZID: chr  "7105" "8813" "57147" "55732" ...

sum(is.na(EnsemblVsGeneSymbol_Human$ENSEMBL))
#[1] 0
sum(is.na(EnsemblVsGeneSymbol_Human$SYMBOL))
#[1] 1614
sum(is.na(EnsemblVsGeneSymbol_Human$ENTREZID))
#[1] 1614

sum(duplicated(EnsemblVsGeneSymbol_Human$ENSEMBL))
#[1] 365
sum(duplicated(EnsemblVsGeneSymbol_Human$SYMBOL))
#[1] 1641
sum(duplicated(EnsemblVsGeneSymbol_Human$ENTREZID))
#[1] 1640

#Important for later joins:
#Many of the ENSEMBL genes don't have associated symbols or ENTREZIDs (1614)
#Some of the ENSEMBL genes map to more than one symbol or ENTREZID (365)
#A small handful of ENSEMBL genes map to the same symbol or ENTREZID (36)

GandalFrontalCortex_RNASeqMeta_wSymbols<-join(GandalFrontalCortex_RNASeqMeta, EnsemblVsGeneSymbol_Human, by="ENSEMBL", type="left", match="all")
str(GandalFrontalCortex_RNASeqMeta_wSymbols)
# 'data.frame':	16188 obs. of  11 variables:
#   $ ENSEMBL      : chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
# $ SCZ.logFC    : num  0.04095 -0.03274 -0.00353 -0.01482 0.07601 ...
# $ SCZ.z        : num  1.235 -1.42 -0.185 -0.537 1.228 ...
# $ SCZ.P.Value  : num  0.217 0.156 0.853 0.591 0.22 ...
# $ SCZ.adj.P.Val: num  0.483 0.4 0.938 0.795 0.486 ...
# $ BD.logFC     : num  0.05429 -0.05279 -0.0432 -0.07828 0.00349 ...
# $ BD.z         : num  1.5091 -1.6643 -1.6483 -1.9198 0.0387 ...
# $ BD.P.Value   : num  0.1313 0.0961 0.0993 0.0549 0.9691 ...
# $ BD.adj.P.Val : num  0.567 0.509 0.514 0.4 0.996 ...
# $ SYMBOL       : chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...
# $ ENTREZID     : chr  "7105" "8813" "57147" "55732" ...

#So that I have this for later:
write.csv(GandalFrontalCortex_RNASeqMeta_wSymbols, "GandalFrontalCortex_RNASeqMeta_wSymbols.csv")

#For the purposes of joining with other databases, I think I should pull out the results from the Ensembl genes that lack other forms of annotation or that map to more than one gene symbol:

GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA<-GandalFrontalCortex_RNASeqMeta_wSymbols[is.na(GandalFrontalCortex_RNASeqMeta_wSymbols$SYMBOL)==FALSE,]
GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped<-GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA[duplicated(GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA$ENSEMBL)==FALSE,]

str(GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped)
#' #'data.frame':	14209 obs. of  11 variables:
#' $ ENSEMBL      : chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
#' $ SCZ.logFC    : num  0.04095 -0.03274 -0.00353 -0.01482 0.07601 ...
#' $ SCZ.z        : num  1.235 -1.42 -0.185 -0.537 1.228 ...
#' $ SCZ.P.Value  : num  0.217 0.156 0.853 0.591 0.22 ...
#' $ SCZ.adj.P.Val: num  0.483 0.4 0.938 0.795 0.486 ...
#' $ BD.logFC     : num  0.05429 -0.05279 -0.0432 -0.07828 0.00349 ...
#' $ BD.z         : num  1.5091 -1.6643 -1.6483 -1.9198 0.0387 ...
#' $ BD.P.Value   : num  0.1313 0.0961 0.0993 0.0549 0.9691 ...
#' $ BD.adj.P.Val : num  0.567 0.509 0.514 0.4 0.996 ...
#' $ SYMBOL       : chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...
#' $ ENTREZID     : chr  "7105" "8813" "57147" "55732" ...

#That data.frame is quite a bit smaller.

write.csv(GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped, "GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped.csv")

GandalFrontalCortex_RNASeqMeta_SchizEffects<-GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped[,c(1:5, 10:11)]
GandalFrontalCortex_RNASeqMeta_BPEffects<-GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped[,c(1,6:9, 10:11)]

str(GandalFrontalCortex_RNASeqMeta_SchizEffects)
# 'data.frame':	14209 obs. of  7 variables:
# $ ENSEMBL      : chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
# $ SCZ.logFC    : num  0.04095 -0.03274 -0.00353 -0.01482 0.07601 ...
# $ SCZ.z        : num  1.235 -1.42 -0.185 -0.537 1.228 ...
# $ SCZ.P.Value  : num  0.217 0.156 0.853 0.591 0.22 ...
# $ SCZ.adj.P.Val: num  0.483 0.4 0.938 0.795 0.486 ...
# $ SYMBOL       : chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...
# $ ENTREZID     : chr  "7105" "8813" "57147" "55732" ...

str(GandalFrontalCortex_RNASeqMeta_BPEffects)
# $ ENSEMBL     : chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
# $ BD.logFC    : num  0.05429 -0.05279 -0.0432 -0.07828 0.00349 ...
# $ BD.z        : num  1.5091 -1.6643 -1.6483 -1.9198 0.0387 ...
# $ BD.P.Value  : num  0.1313 0.0961 0.0993 0.0549 0.9691 ...
# $ BD.adj.P.Val: num  0.567 0.509 0.514 0.4 0.996 ...
# $ SYMBOL      : chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...
# $ ENTREZID    : chr  "7105" "8813" "57147" "55732" ...


##############

#Reading in other comparison datasets:

#Information about Cell Type Specificity from the BrainInABlender database (Hagenauer et al. 2018, suppl table 1)
CellTypeSpecificity<-read.csv("S01_Table_CellTypeSpecificGenes.csv", header=T, stringsAsFactors = FALSE)
PritzkerDLPFC_CorWCellTypeIndex<-read.csv("S08_Table_Pritzker_PsychResultsDiffModels_M4UsualSuspectsAndCellType.csv", header=T, stringsAsFactors = FALSE)


#In my original (much earlier) analysis, I also compared with these two meta-analyses, but they are more weakly powered and less anatomically defined than Gandal et al. 2018:
#Information about Schiz and Bipolar effects in the cortex:
# ChoiMetaAnalysis_BPEffects<-read.csv("S07_CellTypeSpecificEffects_ChoiMetaAnalysis_BP.csv", header=T, stringsAsFactors = FALSE)
# MistryMetaAnalysis_SchizEffects<-read.csv("S07_CellTypeSpecificEffects_MistryMetaAnalysis_Schiz.csv", header=T, stringsAsFactors = FALSE)


#My database of cell type specific diagnosis effects used for validation in Hagenauer et al. 2018 (Suppl table 7):
CellTypeSpecific_DiagnosisEffects<-read.csv("S07_CellTypeSpecificEffects_ParticularCellTypes.csv", header=T, stringsAsFactors = FALSE)


#These results were all from DLPFC analyses run as part of our PLOS cell type paper (Hagenauer et al. 2018, similar pipeline for each dataset, already published):

# Effects for diagnosis and all co-variates from the Pritzker DLPFC microarray dataset (Suppl Table 8):
PritzkerDLPFC_BPSchizEffects<-read.csv("S08_Table_Pritzker_PsychResultsDiffModels_M2UsualSuspects.csv", header=T, stringsAsFactors = FALSE)

# Effects for diagnosis (and baseline expression/intercept) from the CommonMind Consortium RNA-Seq dataset (Suppl Table 9):
CMCDLPFC_BPEffects<-read.csv("S09_Table_CMC_PsychEffectsinDifferentModels_BP_M2UsualSuspects.csv", header=T, stringsAsFactors = FALSE)
CMCDLPFC_SchizEffects<-read.csv("S09_Table_CMC_PsychEffectsinDifferentModels_Schiz_M2UsualSuspects.csv", header=T, stringsAsFactors = FALSE)

#Effects for diagnosis from the Lanz et al. dataset (Suppl Table 11):
LanzDLPFC_BPEffects<-read.csv("S11_Table_Lanz_PsychEffectsinDifferentModels_BP_M2_UsualSuspects.csv", header=T, stringsAsFactors = FALSE)
LanzDLPFC_SchizEffects<-read.csv("S11_Table_Lanz_PsychEffectsinDifferentModels_Schiz_M2_UsualSuspects.csv", header=T, stringsAsFactors = FALSE)

#Effects for diagnosis from the Narayann et al. dataset (Suppl Table 12):
NarayanDLPFC_SchizEffects<-read.csv("S12_Table_Narayan_PsychEffectsinDifferentModels_Schiz_M2UsualSuspects.csv", header=T, stringsAsFactors = FALSE)

#################

#I would like to grab annotation for all genes originally included in the frontal pole qPCR analysis (I'm curious about the ones that were thrown out during QC too)

setwd("~/Documents/Microarray Gen/FrontalPole/Output/GabaGlu")

GabaGlu_Cq_AllSubjects<-read.csv("GabaGlu_AverageCq_NumberNA_PerGene.csv", header=TRUE, stringsAsFactors = FALSE)

setwd("~/Documents/Microarray Gen/FrontalPole/Output/DA5HT")

DA5HT_Cq_AllSubjects<-read.csv("DA5HT_AverageCq_NumberNA_PerGene.csv", header=TRUE, stringsAsFactors = FALSE)

setwd("~/Documents/Microarray Gen/FrontalPole/CompareWDLPFC")

length(row.names(DA5HT_Cq_AllSubjects))
#[1] 48

head(DA5HT_Cq_AllSubjects)
# X DA5HT_AverageCq_PerGene DA5HT_NumberNA_PerGene
# 1 ADRB1                24.09584                      0
# 2 ADRB2                24.94944                      0
# 3  COMT                24.54420                      0
# 4   DBH                28.31647                      0
# 5   DDC                32.98314                      0
# 6  DRD1                25.02894                      0

length(row.names(GabaGlu_Cq_AllSubjects))
#[1] 96

head(GabaGlu_Cq_AllSubjects)
# X GabaGlu_AverageCq_PerGene GabaGlu_NumberNA_PerGene
# 1    ABAT                  21.05391                        0
# 2   ADCY7                  25.50221                        0
# 3  ADORA1                  24.08014                        0
# 4 ADORA2A                  27.24140                        0
# 5 ALDH5A1                  23.20638                        0
# 6     SST                  22.70235                        0

FrontalPole_qPCRGenes<-c(GabaGlu_Cq_AllSubjects$X, DA5HT_Cq_AllSubjects$X)
FrontalPole_qPCRGenes
#[1] "ABAT"    "ADCY7"   "ADORA1"  "ADORA2A" "ALDH5A1" "SST"    ...

str(FrontalPole_qPCRGenes)

FrontalPole_qPCRGenes_MasterAnnotation<-data.frame(GeneSymbol=FrontalPole_qPCRGenes, HousekeepingGenes=FrontalPole_qPCRGenes%in%HousekeepingGenes$GeneSymbol, stringsAsFactors = FALSE)

setwd("~/Documents/Microarray Gen/FrontalPole/Sample_MetaData")


FrontalPole_qPCRGenes_MasterAnnotation<-join(FrontalPole_qPCRGenes_MasterAnnotation, GeneCategories, by="GeneSymbol", type="left")

str(FrontalPole_qPCRGenes_MasterAnnotation)
# 'data.frame':	144 obs. of  4 variables:
#   $ GeneSymbol                : chr  "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# $ HousekeepingGenes         : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ FunctionalCategory_Main   : chr  "GABA" "GLUTAMATE" "OTHER" "OTHER" ...
# $ FunctionalCategory_Adriana: chr  "ENZYMES RELATED TO GABA METABOLISM" "ENZYMES RELATED TO METABOTROPIC GLUTAMATE RECEPTOR ACTIVITY" "ADENOSINE RECEPTOR" "ADENOSINE RECEPTOR" ...

write.csv(FrontalPole_qPCRGenes_MasterAnnotation, "FrontalPole_qPCRGenes_MasterAnnotation.csv")

############

#Adding information about predicted cell type specificity of expression (correlation with cell type index from the Pritzker microarray dataset in Hagenauer et al. 2018):

colnames(PritzkerDLPFC_CorWCellTypeIndex)[2]<-"GeneSymbol"

FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType<-join(FrontalPole_qPCRGenes_MasterAnnotation, PritzkerDLPFC_CorWCellTypeIndex, by="GeneSymbol", type="left")

write.csv(FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType, "FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType.csv")

# Information about Cell Type Specificity from the BrainInABlender database (Hagenauer et al. 2018, suppl table 1)

colnames(CellTypeSpecificity)[4]<-"GeneSymbol"
FrontalPole_qPCRGenes_Annotation_CellType<-join(CellTypeSpecificity, FrontalPole_qPCRGenes_MasterAnnotation, by="GeneSymbol", type="inner")
dim(FrontalPole_qPCRGenes_Annotation_CellType)
#[1] 72 17
write.csv(FrontalPole_qPCRGenes_Annotation_CellType, "FrontalPole_qPCRGenes_Annotation_CellType.csv")

#some of this human annotation seems somewhat funky (the homologies between Human and Mouse in particular). I'm going to try joining directly with mouse using the same name to see if anything more pops up.

colnames(CellTypeSpecificity)[4]<-"GeneSymbol_Human"

CellTypeSpecificity<-data.frame(CellTypeSpecificity, GeneSymbol=toupper(CellTypeSpecificity[,5]))

FrontalPole_qPCRGenes_Annotation_CellTypeMouse<-join(CellTypeSpecificity, FrontalPole_qPCRGenes_MasterAnnotation, by="GeneSymbol", type="inner")
dim(FrontalPole_qPCRGenes_Annotation_CellTypeMouse)
#[1] 67 18
write.csv(FrontalPole_qPCRGenes_Annotation_CellTypeMouse, "FrontalPole_qPCRGenes_Annotation_CellTypeMouse.csv")

#I looked over this output by hand and used it to compile a more abbreviated spreadsheet:

FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType_2<-read.csv("FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType_2.csv", header=T, stringsAsFactors = FALSE)
str(FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType_2)
#'data.frame':	144 obs. of  52 variables:

table(FrontalPole_qPCRGenes_Annotation_PritzkerCovariatesCellType_2$Consensus_CellType)
#                 Astrocyte        Endothelial          Microglia         Neuron_All Neuron_Interneuron  Neuron_Projection 
# 110                  7                  2                  1                  6                 11                  7 

#34 genes are unambiguously related to a particular cell type.


##################################################

#How about relationships with diagnosis and co-variates?

FrontalPole_qPCRGenes_MasterAnnotation

colnames(PritzkerDLPFC_BPSchizEffects)[2]<-"GeneSymbol"
colnames(CMCDLPFC_BPEffects)[2]<-"GeneSymbol"
colnames(CMCDLPFC_SchizEffects)[2]<-"GeneSymbol"
CMCDLPFC_SchizEffects<-CMCDLPFC_SchizEffects[,-c(10:28)]
str(NarayanDLPFC_SchizEffects)

#If I'm joining these all together, I should probably annotate the colnames better:
colnames(CellTypeSpecific_DiagnosisEffects)<-paste("CellTypeSpecific", colnames(CellTypeSpecific_DiagnosisEffects), sep="_")
colnames(CellTypeSpecific_DiagnosisEffects)[2]<-"GeneSymbol"

str(PritzkerDLPFC_BPSchizEffects)
# 'data.frame':	11979 obs. of  29 variables:
#   $ Pritzker_Probe                                              : chr  "10_at" "100_at" "1000_at" "10000_at" ...
# $ GeneSymbol                                                  : chr  "NAT2" "ADA" "CDH2" "AKT3" ...
colnames(PritzkerDLPFC_BPSchizEffects)<-paste("Pritzker", colnames(PritzkerDLPFC_BPSchizEffects), sep="_")
colnames(PritzkerDLPFC_BPSchizEffects)[2]<-"Pritzker_GeneSymbol"
temp<-data.frame(do.call('rbind', strsplit(as.character(PritzkerDLPFC_BPSchizEffects$Pritzker_Probe),'_',fixed=TRUE)), stringsAsFactors = FALSE)
# X1 X2    X3
# 1    10 at    10
# 2   100 at   100
# 3  1000 at  1000
# 4 10000 at 10000
# 5 10001 at 10001
# 6 10002 at 10002

str(temp)
# 'data.frame':	11979 obs. of  3 variables:
#   $ X1: chr  "10" "100" "1000" "10000" ...
# $ X2: chr  "at" "at" "at" "at" ...
# $ X3: chr  "10" "100" "1000" "10000" ...

PritzkerDLPFC_BPSchizEffects$ENTREZID<-temp$X1
head(PritzkerDLPFC_BPSchizEffects)

# colnames(ChoiMetaAnalysis_BPEffects)<-paste("Choi_BP", colnames(ChoiMetaAnalysis_BPEffects), sep="_")
# colnames(ChoiMetaAnalysis_BPEffects)[2]<-"GeneSymbol"
# 
# colnames(MistryMetaAnalysis_SchizEffects)<-paste("Mistry_Schiz", colnames(MistryMetaAnalysis_SchizEffects), sep="_")
# colnames(MistryMetaAnalysis_SchizEffects)[2]<-"GeneSymbol"


#These individual datasets were included in the Gandal Meta-analysis but have other useful information associated with them:

#Includes estimates of average expression levels:

colnames(CMCDLPFC_BPEffects)<-paste("CMC_BP", colnames(CMCDLPFC_BPEffects), sep="_")
colnames(CMCDLPFC_BPEffects)[2]<-"CMC_BP_GeneSymbol"
colnames(CMCDLPFC_BPEffects)[1]<-"ENSEMBL"

colnames(CMCDLPFC_SchizEffects)<-paste("CMC_Schiz", colnames(CMCDLPFC_SchizEffects), sep="_")
colnames(CMCDLPFC_SchizEffects)[2]<-"CMC_BP_GeneSymbol"
colnames(CMCDLPFC_SchizEffects)[1]<-"ENSEMBL"

sum(is.na(CMCDLPFC_BPEffects$ENSEMBL))
#[1] 0

#From Grey Matter only dissections:

colnames(LanzDLPFC_BPEffects)<-paste("Lanz_BP", colnames(LanzDLPFC_BPEffects), sep="_")
colnames(LanzDLPFC_BPEffects)[3]<-"Lanz_BP_GeneSymbol"
colnames(LanzDLPFC_BPEffects)[2]<-"ENTREZID"

colnames(LanzDLPFC_SchizEffects)<-paste("Lanz_Schiz", colnames(LanzDLPFC_SchizEffects), sep="_")
colnames(LanzDLPFC_SchizEffects)[3]<-"Lanz_Schiz_GeneSymbol"
colnames(LanzDLPFC_SchizEffects)[2]<-"ENTREZID"

colnames(NarayanDLPFC_SchizEffects)<-paste("Narayan_Schiz", colnames(NarayanDLPFC_SchizEffects), sep="_")
colnames(NarayanDLPFC_SchizEffects)[3]<-"Narayan_Schiz_GeneSymbol"
colnames(NarayanDLPFC_SchizEffects)[2]<-"ENTREZID"


colnames(GandalCortex_MicroarrayMeta_BPEffects)<-paste("Gandal_MicroarrayMeta_BP", colnames(GandalCortex_MicroarrayMeta_BPEffects), sep="_")
colnames(GandalCortex_MicroarrayMeta_BPEffects)[3]<-"GeneSymbol"
colnames(GandalCortex_MicroarrayMeta_BPEffects)[1]<-"ENSEMBL"
colnames(GandalCortex_MicroarrayMeta_BPEffects)[4]<-"ENTREZID"

colnames(GandalCortex_MicroarrayMeta_SchizEffects)<-paste("Gandal_MicroarrayMeta_Schiz", colnames(GandalCortex_MicroarrayMeta_SchizEffects), sep="_")
colnames(GandalCortex_MicroarrayMeta_SchizEffects)[3]<-"GeneSymbol"
colnames(GandalCortex_MicroarrayMeta_SchizEffects)[1]<-"ENSEMBL"
colnames(GandalCortex_MicroarrayMeta_SchizEffects)[4]<-"ENTREZID"

colnames(GandalCortex_MicroarrayMeta_MDDEffects)<-paste("Gandal_MicroarrayMeta_MDD", colnames(GandalCortex_MicroarrayMeta_MDDEffects), sep="_")
colnames(GandalCortex_MicroarrayMeta_MDDEffects)[3]<-"GeneSymbol"
colnames(GandalCortex_MicroarrayMeta_MDDEffects)[1]<-"ENSEMBL"
colnames(GandalCortex_MicroarrayMeta_MDDEffects)[4]<-"ENTREZID"

#On second thought, maybe I should just do this:

colnames(GandalCortex_MicroarrayMeta)<-paste("Gandal_MicroarrayMeta", colnames(GandalCortex_MicroarrayMeta), sep="_")
colnames(GandalCortex_MicroarrayMeta)[1]<-"ENSEMBL"
colnames(GandalCortex_MicroarrayMeta)[3]<-"Gandal_MicroarrayMeta_GeneSymbol"
colnames(GandalCortex_MicroarrayMeta)[4]<-"Gandal_MicroarrayMeta_ENTREZID"

sum(is.na(GandalCortex_MicroarrayMeta$ENSEMBL))
#[1] 0

colnames(GandalFrontalCortex_RNASeqMeta_BPEffects)<-paste("Gandal_RNASeqMeta_BP", colnames(GandalFrontalCortex_RNASeqMeta_BPEffects), sep="_")
colnames(GandalFrontalCortex_RNASeqMeta_BPEffects)[6]<-"GeneSymbol"
colnames(GandalFrontalCortex_RNASeqMeta_BPEffects)[1]<-"ENSEMBL"
colnames(GandalFrontalCortex_RNASeqMeta_BPEffects)[7]<-"ENTREZID"

colnames(GandalFrontalCortex_RNASeqMeta_SchizEffects)<-paste("Gandal_RNASeqMeta_Schiz", colnames(GandalFrontalCortex_RNASeqMeta_SchizEffects), sep="_")
colnames(GandalFrontalCortex_RNASeqMeta_SchizEffects)[6]<-"GeneSymbol"
colnames(GandalFrontalCortex_RNASeqMeta_SchizEffects)[1]<-"ENSEMBL"
colnames(GandalFrontalCortex_RNASeqMeta_SchizEffects)[7]<-"ENTREZID"

#On second thought - let's use the original data.frame (before removing multi-mapped and NA)
colnames(GandalFrontalCortex_RNASeqMeta_wSymbols)<-paste("Gandal_RNASeqMeta", colnames(GandalFrontalCortex_RNASeqMeta_wSymbols), sep="_")
colnames(GandalFrontalCortex_RNASeqMeta_wSymbols)[1]<-"ENSEMBL"
colnames(GandalFrontalCortex_RNASeqMeta_wSymbols)[10]<-"Gandal_RNASeqMeta_GeneSymbol"
colnames(GandalFrontalCortex_RNASeqMeta_wSymbols)[11]<-"Gandal_RNASeqMeta_ENTREZID"

#Oh wait, for joining purposes, I still need to remove the ENSEMBL IDs that are duplicated. It should be o.k. for now, because this is the only ENSEMBL dataset that has those duplications.

#Test case
# Temp<-join_all(dfs=list(FrontalPole_qPCRGenes_MasterAnnotation, CellTypeSpecific_DiagnosisEffects), by="GeneSymbol", type="left", match="first")
# str(Temp)
# rm(Temp)

#Note - this code definitely isn't as efficient as it could have been (since I broke things down into diagnosis-specific data.frames!)

length(CMCDLPFC_BPEffects$ENSEMBL)
#[1] 22053

length(GandalFrontalCortex_RNASeqMeta_wSymbols$ENSEMBL)
#[1] 16188

sum(GandalFrontalCortex_RNASeqMeta_wSymbols$ENSEMBL%in%CMCDLPFC_BPEffects$ENSEMBL)
#[1] 15738

#hmmm... the CMC is older annotation, so I don't want to throw things out that are missing in CMC
#But Gandal has a smaller number of ENSEMBL annotated genes, and I suspect that is because they used a more conservative filter. Since we are targeting low-level expressed transcripts, I don't want to use that as a filter either.
#But match "all" is going to produce problems because there are duplicated (multimapped) ENSEMBL genes in Gandal.
#Maybe we should throw those out after all. Let's see which genes they are:

GandalFrontalCortex_RNASeqMeta_wSymbols$GeneSymbol[duplicated(GandalFrontalCortex_RNASeqMeta_wSymbols$ENSEMBL)]
#Looks like a whole bunch of MIR and LOC genes. For the sake of a clean join, I'm going to dump them:

GandalFrontalCortex_RNASeqMeta_wSymbols_NoMultiMapped<-GandalFrontalCortex_RNASeqMeta_wSymbols[duplicated(GandalFrontalCortex_RNASeqMeta_wSymbols$ENSEMBL)==FALSE,]

#I'm going to join everything that has ENSEMBL as its primary annotation first:
Gandal_Microarray_RNASeq_W_CMCExpression<-join_all(dfs=list(GandalFrontalCortex_RNASeqMeta_wSymbols_NoMultiMapped, GandalCortex_MicroarrayMeta, CMCDLPFC_BPEffects, CMCDLPFC_SchizEffects), by="ENSEMBL", type="full", match="all")

str(Gandal_Microarray_RNASeq_W_CMCExpression)
#'data.frame':	31500 obs. of  49 variables:
#hmmm... that seems concerning. 10,000 genes were added. Maybe the annotation for one of the datasets isn't aligning with the others due ot a formatting error?
sum(duplicated(Gandal_Microarray_RNASeq_W_CMCExpression))
#[1] 0

#Nope - outputted and double-checked and there doesn't appear to be a joining problem, there are just different sets of Ensembl genes included in the microarray dataset vs. the two RNA-Seq datasets.  This makes some sense - microarray would not be limited to genes that are brain-expressed.


#hmmm... I wonder though how it is dealing with the fact that the ENTREZID and GeneSymbol columns have the same names - which version is it choosing when it combines them?
#Would it be worthwhile to not use the original ENTREZID and GeneSymbol (which may not match between datasets), and just re-annotate?
#Yeah, let's do that:

#Adding more up-to-date annotation (again)

uniKeys <- head(keys(org.Hs.eg.db, keytype="ENSEMBL"))
cols <- c("SYMBOL", "ENTREZID")
select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")
# 'select()' returned 1:1 mapping between keys and columns
# ENSEMBL   SYMBOL ENTREZID
# 1 ENSG00000121410     A1BG        1
# 2 ENSG00000175899      A2M        2
# 3 ENSG00000256069    A2MP1        3
# 4 ENSG00000171428     NAT1        9
# 5 ENSG00000156006     NAT2       10
# 6 ENSG00000196136 SERPINA3       12

#Lovely.

colnames(Gandal_Microarray_RNASeq_W_CMCExpression)[1]
#[1] "ENSEMBL"

EnsemblVsGeneSymbol_Human<-select(org.Hs.eg.db, keys=Gandal_Microarray_RNASeq_W_CMCExpression$ENSEMBL, columns=cols, keytype="ENSEMBL")
#'select()' returned 1:many mapping between keys and columns


Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols<-join(Gandal_Microarray_RNASeq_W_CMCExpression, EnsemblVsGeneSymbol_Human, by="ENSEMBL", type="left", match="all")
str(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols)

# 'data.frame':	32050 obs. of  51 variables:
#   $ ENSEMBL                               : chr  "ENSG00000000003" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
# $ Gandal_RNASeqMeta_SCZ.logFC           : num  0.04095 -0.03274 -0.00353 -0.01482 0.07601 ...
# $ Gandal_RNASeqMeta_SCZ.z               : num  1.235 -1.42 -0.185 -0.537 1.228 ...
# $ Gandal_RNASeqMeta_SCZ.P.Value         : num  0.217 0.156 0.853 0.591 0.22 ...

colnames(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols)

sum(duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL))
#[1] 550
#A few of the ENSEMBL genes are mapped to more than one gene symbol.

sum(duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL))
#[1] 9724

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL))
#[1] 9541

#Very few ENSEMBL genes mapping to the same gene symbol, but many are mapping to no symbol.

sum(duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID))
#[1] 9723

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID))
#[1] 9541

#How many of these are driven by the CMC data having older annotation?

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$CMC_Schiz_logFC_M2)==FALSE)
#[1] 22508
#That's bigger than the original CMC data frame by 
22508-22053 
#[1] 455
#So there were 455 genes in CMC where the Ensembl annotation mapped to multiple other annotations (Entrez, GeneSymbol)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL) & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$CMC_Schiz_logFC_M2)==FALSE)
#[1] 5056
#And 5056 of the CMC genes don't have up-to date symbol annotation

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$CMC_Schiz_logFC_M2)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 16997
#number of Ensembl annotated genes in CMC dataset that map to unique gene symbols (not multimapped)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$CMC_Schiz_logFC_M2)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 16997
#number of Ensembl annotated genes in CMC dataset that map to unique Entrez IDs (not multimapped)


sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_MDD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 22612
#Number of original Ensembl-annotated results for the MDD meta-analysis

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_MDD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 18981
#number of Ensembl annotated genes in MDD microarray meta-analysis dataset that map to unique Gene symbol (not multimapped)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_MDD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 18981
#number of Ensembl annotated genes in MDD microarray meta-analysis dataset that map to unique Entrez IDs (not multimapped)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_AAD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 22221
#Number of original Ensembl-annotated results for the AAD meta-analysis

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_AAD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 18844
#number of Ensembl annotated genes in AAD microarray meta-analysis dataset that map to unique Gene symbol (not multimapped)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_AAD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 18844
#number of Ensembl annotated genes in AAD microarray meta-analysis dataset that map to unique EntrezID (not multimapped)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_ASD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 15640
#number of Ensembl annotated genes in ASD microarray meta-analysis dataset that map to unique Gene symbol (not multimapped)
sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID)==FALSE & is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$Gandal_MicroarrayMeta_ASD.beta_log2FC)==FALSE & duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL)==FALSE)
#[1] 15640
#number of Ensembl annotated genes in ASD microarray meta-analysis dataset that map to unique EntrezID (not multimapped)


setwd("~/Documents/Microarray Gen/FrontalPole/CompareWDLPFC")

write.csv(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols, "Gandal_Microarray_RNASeq_W_CMCExpression.csv")



#Next, I'm going to join all of the datasets that use ENTREZID as their primary annotation:

sum(is.na(PritzkerDLPFC_BPSchizEffects$ENTREZID))
#[1] 0

sum(duplicated(PritzkerDLPFC_BPSchizEffects$ENTREZID))
#[1] 11
#That's weird.
PritzkerDLPFC_BPSchizEffects$ENTREZID[duplicated(PritzkerDLPFC_BPSchizEffects$ENTREZID)]
# [1] "AFFX-HSAC07/X00351"    "AFFX-HSAC07/X00351"    "AFFX-HUMGAPDH/M33197"  "AFFX-HUMGAPDH/M33197"  "AFFX-HUMISGF3A/M97935"
# [6] "AFFX-HUMISGF3A/M97935" "AFFX-HUMISGF3A/M97935" "AFFX-HUMRGE/M10098"    "AFFX-HUMRGE/M10098"    "AFFX-M27830"          
# [11] "AFFX-M27830"
#Ah - they are affymetrix control probes - not useful anyway.

sum(is.na(LanzDLPFC_BPEffects$ENTREZID))
#[1] 0
LanzDLPFC_SchizEffects

sum(is.na(NarayanDLPFC_SchizEffects$ENTREZID))
#[1] 0

DLPFC_LanzNarayanPrizker<-join_all(dfs=list(LanzDLPFC_BPEffects, LanzDLPFC_SchizEffects, NarayanDLPFC_SchizEffects, PritzkerDLPFC_BPSchizEffects), by="ENTREZID", type="full", match="all")

str(DLPFC_LanzNarayanPrizker)
#'data.frame':	20005 obs. of  48 variables:

sum(duplicated(DLPFC_LanzNarayanPrizker$ENTREZID))
#[1] 11


#Adding more up-to-date gene symbol annotation:

DLPFC_LanzNarayanPrizker_wSymbols<-join(DLPFC_LanzNarayanPrizker, EnsemblVsGeneSymbol_Human, by="ENTREZID", type="left", match="all")
str(DLPFC_LanzNarayanPrizker_wSymbols)
#'data.frame':	20108 obs. of  50 variables:

sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$SYMBOL))
#[1] 1419
sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$ENSEMBL))
#[1] 1419

sum(duplicated(DLPFC_LanzNarayanPrizker_wSymbols$SYMBOL))
#[1] 1521
sum(duplicated(DLPFC_LanzNarayanPrizker_wSymbols$ENSEMBL))
#[1] 1456
#Looks like most of the ENTREZ-ENSEMBL mappings are one-to-one...except when they are missing.

#for documentation: The final size of the annotated datasets for each study:

sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$SYMBOL)==FALSE & duplicated(DLPFC_LanzNarayanPrizker_wSymbols$ENTREZID)==FALSE & is.na(DLPFC_LanzNarayanPrizker_wSymbols$Lanz_BP_logFC_M2)==FALSE)
#[1] 18480
#Number of Entrez IDs mapping to unique symbols in the Lanz dataset

sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$ENSEMBL)==FALSE & duplicated(DLPFC_LanzNarayanPrizker_wSymbols$ENTREZID)==FALSE & is.na(DLPFC_LanzNarayanPrizker_wSymbols$Lanz_BP_logFC_M2)==FALSE)
#[1] 18480
#After doing this multiple times, I'm pretty sure the answer for both symbol and other annotation is the same.

sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$SYMBOL)==FALSE & duplicated(DLPFC_LanzNarayanPrizker_wSymbols$ENTREZID)==FALSE & is.na(DLPFC_LanzNarayanPrizker_wSymbols$Narayan_Schiz_adj.P.Val_M2)==FALSE)
#[1] 18480
#Number of Entrez IDs mapping to unique symbols in the Narayan dataset (same as Lanz - same chip, same previous analysis)

sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$SYMBOL)==FALSE & duplicated(DLPFC_LanzNarayanPrizker_wSymbols$ENTREZID)==FALSE & is.na(DLPFC_LanzNarayanPrizker_wSymbols$Pritzker_Beta_Diagnosis..MDD.vs..Control.baseline.)==FALSE)
#[1] 11668
#Number of Entrez IDs mapping to unique symbols in the Pritzker dataset (older chip)

write.csv(DLPFC_LanzNarayanPrizker_wSymbols, "DLPFC_LanzNarayanPrizker.csv")


#For documention purposes - the Lanz and Narayan datasets are tightly correlated despite having small sample sizes and being from separate brain banks (i.e., very high quality)
cor(DLPFC_LanzNarayanPrizker_wSymbols$Lanz_Schiz_logFC_M2, DLPFC_LanzNarayanPrizker_wSymbols$Narayan_Schiz_logFC_M2, use="pairwise.complete.obs")
#[1] 0.4576384




#It might be better (as far as avoiding duplicates doing repetitive joining), to try to join these two data frames by ENSEMBL before doing a join with the qPCR using GeneSymbol.  Let's snoop:

nrow(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols)
#[1] 32050

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENSEMBL))
#[1] 0

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID))
#[1] 9541

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL))
#[1] 9541

sum(duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL))
#[1] 9724


Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols<-join(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols, DLPFC_LanzNarayanPrizker_wSymbols, by="ENSEMBL", type="full", match="all")

str(Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols)
#'data.frame':	33548 obs. of  98 variables:

sum(duplicated(Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols$ENSEMBL))
#[1] 2047

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols$ENSEMBL))
#[1] 1419

#Not bad. The join went smoothly - since there weren't any NAs in the ENSEMBL genes in the CMC/Gandal data.frame, only the NAs from the other data.frame remain in the joined data.frame (no exponential increase)

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols$SYMBOL))
#[1] 10960

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols$ENTREZID))
#[1] 9541

#Around a third of these entries are missing other forms of annotation. 

write.csv(Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols, "Gandal_Microarray_RNASeq_W_CMCExpression_W_LanzNarayanPrizker_wSymbols.csv")


#Pulling out the results for the Frontal Pole qPCR genes:

colnames(FrontalPole_qPCRGenes_MasterAnnotation)[1]<-"SYMBOL"
colnames(CellTypeSpecific_DiagnosisEffects)[2]<-"SYMBOL"

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects<-join_all(dfs=list(FrontalPole_qPCRGenes_MasterAnnotation, CellTypeSpecific_DiagnosisEffects, Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols, DLPFC_LanzNarayanPrizker_wSymbols), by="SYMBOL", type="left", match="all")


# FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects<-join_all(dfs=list(FrontalPole_qPCRGenes_MasterAnnotation, CellTypeSpecific_DiagnosisEffects, PritzkerDLPFC_BPSchizEffects, ChoiMetaAnalysis_BPEffects, CMCDLPFC_BPEffects, LanzDLPFC_BPEffects, MistryMetaAnalysis_SchizEffects, CMCDLPFC_SchizEffects, LanzDLPFC_SchizEffects, NarayanDLPFC_SchizEffects), by="GeneSymbol", type="left", match="first")
# str(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects)

write.csv(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects, "FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects.csv")


# PreviousData_BipolarPval<-cbind(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Pritzker_Pval_Nominal_Diagnosis_BP, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Choi_BP_NominalP, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$CMC_BP_P.Value_M2,FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Lanz_BP_P.Value_M2)

PreviousData_BipolarPval<-cbind(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_RNASeqMeta_BD.P.Value, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_MicroarrayMeta_BD.P.value, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Lanz_BP_P.Value_M2)
#I should probably change this - it turns out that the Gandal Microarray meta-analysis overlaps with our frontal pole microarray samples and with Lanz
#I should also add our frontal pole microarray results into this.

PreviousData_BipolarPval_Minimum<-apply(PreviousData_BipolarPval, 1, function(y) min(y, na.rm=TRUE))
# Warning message:
#   In min(y, na.rm = TRUE) : no non-missing arguments to min; returning Inf
PreviousData_BipolarPval_Minimum[PreviousData_BipolarPval_Minimum>1]<-NA

pdf("Histogram_Bipolar_PreviousMinimumPval.pdf", width=4, height=4)
hist(PreviousData_BipolarPval_Minimum, breaks=20, main="", xlim=c(0,1))
dev.off()
#Strongly left-skewed - could be an artifact of looking for the min (i.e., it would be left skewed even if the pvals were random)
#I wonder if HK show the same skew?
pdf("Histogram_Bipolar_PreviousMinimumPval_HK.pdf", width=4, height=4)
hist(PreviousData_BipolarPval_Minimum[FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$HousekeepingGenes==TRUE], main="",xlim=c(0,1))
dev.off()
#Yep, looks pretty similar.

sum(PreviousData_BipolarPval_Minimum<0.05, na.rm=TRUE)
#[1] 68
length(PreviousData_BipolarPval_Minimum)
#[1] 144

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$PreviousData_BipolarPval_Minimum<-PreviousData_BipolarPval_Minimum

#Maybe instead I should count the number of p<0.05?  (btw, this is a really inadequate form of meta-anaysis - two datasets could have different direction of effect! - but a good start)

PreviousData_BipolarPval_TotalSigPval<-apply(PreviousData_BipolarPval, 1, function(y) sum(y<0.05, na.rm=TRUE))
length(PreviousData_BipolarPval_TotalSigPval)
#[1] 144
table(PreviousData_BipolarPval_TotalSigPval)
# 0  1  2  3 
# 76 48 17  3
#Note that this is partially circular, since Lanz was included in the microarray meta-analysis.

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$SYMBOL[which(PreviousData_BipolarPval_TotalSigPval>1)]
# [1] "SST"     "PVALB"   "BDNF"    "GABRA1"  "GABRG1"  "GABRG3"  "GAD1"    "GLS"     "GRIA4"   "GRIK1"   "GRIN2C"  "GRM2"    "HOMER1" 
#[14] "NSF"     "AQP4"    "SLC7A11" "PGK1"    "HTR2A"   "MAOA"    "YWHAZ" 

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$PreviousData_BipolarPval_TotalSigPval<-PreviousData_BipolarPval_TotalSigPval


PreviousData_BipolarFDR<-cbind(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_MicroarrayMeta_BD.FDR, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_RNASeqMeta_BD.adj.P.Val,FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Lanz_BP_adj.P.Val_M2)

PreviousData_BipolarFDR_Minimum<-apply(PreviousData_BipolarFDR, 1, function(y) min(y, na.rm=TRUE))
# Warning message:
#   In min(y, na.rm = TRUE) : no non-missing arguments to min; returning Inf
PreviousData_BipolarFDR_Minimum[PreviousData_BipolarFDR_Minimum>1]<-NA
length(PreviousData_BipolarFDR_Minimum)
#[1] 144

sum(PreviousData_BipolarFDR_Minimum<0.05, na.rm=TRUE)
#[1] 13

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$SYMBOL[which(PreviousData_BipolarFDR_Minimum<0.05)]
#[1] "SST"     "PVALB"   "BDNF"    "GABRA4"  "GABRG3"  "GRIK1"   "HOMER1"  "GJA1"    "SLC1A3"  "HTR2A"   "MAOA"    "SLC18A2" "PPIA"  

sum(PreviousData_BipolarFDR_Minimum<0.10, na.rm=TRUE)
#[1] 24

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$SYMBOL[which(PreviousData_BipolarFDR_Minimum<0.10)]
# [1] "SST"     "PVALB"   "BDNF"    "GABBR2"  "GABRA4"  "GABRD"   "GABRG3"  "GABRR1"  "GAD1"    "GLUL"    "GRIA2"   "GRIK1"   "HOMER1" 
#[14] "GJA1"    "SLC1A2"  "SLC1A3"  "SLC6A1"  "SLC7A11" "HMBS"    "GAPDH"   "HTR2A"   "MAOA"    "SLC18A2" "PPIA"   

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$PreviousData_BipolarFDR_Minimum<-PreviousData_BipolarFDR_Minimum


#Same analysis for Schiz:

PreviousData_SchizPval<-cbind(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_RNASeqMeta_SCZ.P.Value, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_MicroarrayMeta_SCZ.P.value, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Lanz_Schiz_P.Value_M2, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Narayan_Schiz_P.Value_M2)

PreviousData_SchizPval_Minimum<-apply(PreviousData_SchizPval, 1, function(y) min(y, na.rm=TRUE))
# Warning message:
#   In min(y, na.rm = TRUE) : no non-missing arguments to min; returning Inf
PreviousData_SchizPval_Minimum[PreviousData_SchizPval_Minimum>1]<-NA

pdf("Histogram_Schiz_PreviousMinimumPval.pdf", width=4, height=4)
hist(PreviousData_SchizPval_Minimum, breaks=20, main="", xlim=c(0,1))
dev.off()
#Strongly left-skewed - could be an artifact of looking for the min (i.e., it would be left skewed even if the pvals were random)
#I wonder if HK show the same skew?
pdf("Histogram_Schiz_PreviousMinimumPval_HK.pdf", width=4, height=4)
hist(PreviousData_SchizPval_Minimum[FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$HousekeepingGenes==TRUE], main="",xlim=c(0,1))
dev.off()
#Yep, looks pretty similar.

sum(PreviousData_SchizPval_Minimum<0.05, na.rm=TRUE)
#[1] 80
length(PreviousData_SchizPval_Minimum)
#[1] 144

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$PreviousData_SchizPval_Minimum<-PreviousData_SchizPval_Minimum

#Maybe instead I should count the number of p<0.05?  (btw, this is a really inadequate form of meta-anaysis - two datasets could have different direction of effect! - but a good start)

PreviousData_SchizPval_TotalSigPval<-apply(PreviousData_SchizPval, 1, function(y) sum(y<0.05, na.rm=TRUE))
length(PreviousData_SchizPval_TotalSigPval)
#[1] 144
table(PreviousData_SchizPval_TotalSigPval)
# 0  1  2  3  4 
# 64 41 25 11  3
#Note that this is partially circular, since Lanz and Narayan were included in the microarray meta-analysis.


FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$SYMBOL[which(PreviousData_SchizPval_TotalSigPval>1)]
# [1] "ADORA1"  "ALDH5A1" "SST"     "PVALB"   "CALB1"   "GABBR1"  "GABBR2"  "GABRA2"  "GABRB1"  "GABRD"   "GABRE"   "GABRG1"  "GABRG3" 
# [14] "GLS"     "GRIA4"   "GRIK1"   "GRIK2"   "GRIN2C"  "NSF"     "PHGDH"   "AQP4"    "GJA1"    "S100B"   "SLC17A6" "SLC17A7" "SLC1A1"
# [27] "SLC1A3"  "SLC6A1"  "SLC6A11" "SLC6A12" "SLC7A11" "HMBS"    "TBP"     "DRD5"    "GAPDH"   "HTR1F"   "MAOA"    "MAOB"    "PPIA" 

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$PreviousData_SchizPval_TotalSigPval<-PreviousData_SchizPval_TotalSigPval
  

PreviousData_SchizFDR<-cbind(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_RNASeqMeta_SCZ.adj.P.Val, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Gandal_MicroarrayMeta_SCZ.FDR,FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Lanz_Schiz_adj.P.Val_M2, FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$Narayan_Schiz_adj.P.Val_M2)

PreviousData_SchizFDR_Minimum<-apply(PreviousData_SchizFDR, 1, function(y) min(y, na.rm=TRUE))
# Warning message:
#   In min(y, na.rm = TRUE) : no non-missing arguments to min; returning Inf
PreviousData_SchizFDR_Minimum[PreviousData_SchizFDR_Minimum>1]<-NA
length(PreviousData_SchizFDR_Minimum)
#[1] 144

sum(PreviousData_SchizFDR_Minimum<0.05, na.rm=TRUE)
#[1] 38

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$SYMBOL[which(PreviousData_SchizFDR_Minimum<0.05)]
# [1] "SST"     "PVALB"   "BDNF"    "CALB1"   "GABBR2"  "GABRA4"  "GABRB1"  "GABRD"   "GABRE"   "GABRG3"  "GLS"     "GLUL"    "GRIA2"  
# [14] "GRIA4"   "GRIK1"   "GRIN2C"  "GRM3"    "HOMER1"  "HOMER2"  "PHGDH"   "AQP4"    "GJA1"    "S100B"   "SLC17A6" "SLC1A2"  "SLC1A3"
#[27] "SLC38A1" "SLC6A1"  "SLC6A11" "SLC6A13" "SLC7A11" "HMBS"    "DRD1"    "HTR1F"   "HTR2A"   "MAOA"    "MAOB"    "PPIA"  

sum(PreviousData_SchizFDR_Minimum<0.10, na.rm=TRUE)
#[1] 52

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$SYMBOL[which(PreviousData_SchizFDR_Minimum<0.10)]
# # [1] "ABAT"    "ADORA1"  "ALDH5A1" "SST"     "PVALB"   "BDNF"    "CALB1"   "GABBR1"  "GABBR2"  "GABRA2"  "GABRA4"  "GABRB1"  "GABRD"  
# [14] "GABRE"   "GABRG3"  "GAD1"    "GLS"     "GLUL"    "GRIA2"   "GRIA4"   "GRIK1"   "GRIN1"   "GRIN2C"  "GRM2"    "GRM3"    "GRM4"   
# [27] "GRM5"    "HOMER1"  "HOMER2"  "NSF"     "PHGDH"   "AQP4"    "GJA1"    "S100B"   "SLC17A6" "SLC17A7" "SLC1A2"  "SLC1A3"  "SLC38A1"
# [40] "SLC6A1"  "SLC6A11" "SLC6A12" "SLC6A13" "SLC7A11" "HMBS"    "DRD1"    "DRD5"    "HTR1F"   "HTR2A"   "MAOA"    "MAOB"    "PPIA"  

FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects$PreviousData_SchizFDR_Minimum<-PreviousData_SchizFDR_Minimum

#Since the direction of effect is recorded using a variety of different mechanisms, it would probably be easiest for me just to look it up by hand.

# FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects<-data.frame(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects, PreviousData_BipolarPval_Minimum, PreviousData_BipolarPval_TotalSigPval, PreviousData_BipolarFDR_Minimum, PreviousData_SchizPval_Minimum, PreviousData_SchizPval_TotalSigPval, PreviousData_SchizFDR_Minimum)

write.csv(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects, "FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects_2.csv")

dim(FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects)

#I annotated this to indicate direction of effect and which datasets were showing it.
#Then went back and outputted graphs for the genes that had strong evidence for diagnosis effects based on previous datasets.

##############################

#I trimmed down the columns to the ones that were most relevant, added back in the qPCR results, and then made a copy of the spreadsheet in a format I could easily read back in to R:

setwd("~/Documents/Microarray Gen/FrontalPole/CompareWDLPFC")
FrontalPolePCRResults_vsDLPFCResults_forGraphing<-read.csv("FrontalPolePCRResults_vsDLPFCResults_forR.csv", header=TRUE, stringsAsFactors = FALSE)

dim(FrontalPolePCRResults_vsDLPFCResults_forGraphing)
#[1] 111 162


#Schiz in Frontal Pole qPCR vs. DLPFC RNA-Seq meta-analysis:

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_RNASeqMeta.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_SCZ.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32401 -0.03596  0.01001  0.04409  0.30428 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.032804   0.009602   3.416 0.000943 ***
#   Gandal_RNASeqMeta_SCZ.logFC 0.769006   0.104229   7.378  6.6e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08787 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.3692,	Adjusted R-squared:  0.3624 
# F-statistic: 54.44 on 1 and 93 DF,  p-value: 6.604e-11


#Filtered by p-value: - this code is wrong and needs to be re-done (uses microarray meta-analysis p-value as the filter!)

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_RNASeqMeta_P05.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_SCZ.P.Value<0.05,], main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_SCZ.P.Value<0.05,])
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)
# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_SCZ.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_SCZ.P.Value < 
#                                                                0.05, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18137 -0.03371 -0.00972  0.01520  0.32824 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                  0.04345    0.01848   2.351  0.02629 * 
#   Gandal_RNASeqMeta_SCZ.logFC  0.51016    0.14446   3.531  0.00151 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09421 on 27 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.316,	Adjusted R-squared:  0.2906 
# F-statistic: 12.47 on 1 and 27 DF,  p-value: 0.001506

#It might be more efficient to add color to the other plot to illustrate which effects had been previously observed in the DLPFC
#Also, I'm going to change the y-axis to be limited at -0.6 (currently it extends farther due to HTR2B, but that datapoint isn't shared between the two datasets and isn't pictured)

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_RNASeqMeta_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC", ylim=c(-0.6, 0.6), col="gray50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_SCZ.P.Value<0.05,], col=1, ylim=c(-0.6, 0.6))
abline(TempLine, col=2, lwd=3)
dev.off()


pdf("Plot_BP_LogFC_qPCR_vs_Gandal_RNASeqMeta.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_BD.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.247480 -0.049599  0.005109  0.058883  0.168962 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.045260   0.008501   5.324 7.01e-07 ***
#   Gandal_RNASeqMeta_BD.logFC 0.536928   0.099003   5.423 4.61e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08187 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2403,	Adjusted R-squared:  0.2321 
# F-statistic: 29.41 on 1 and 93 DF,  p-value: 4.606e-07


#Filtering by p-value:

pdf("Plot_BP_LogFC_qPCR_vs_Gandal_RNASeqMeta_P05.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_BD.P.Value<0.05,], main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_BD.P.Value<0.05,])
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_BD.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_BD.P.Value < 
#                                                                0.05, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.15205 -0.06127  0.01910  0.04207  0.12997 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.02246    0.02069   1.086 0.293780    
# Gandal_RNASeqMeta_BD.logFC  0.58524    0.12053   4.855 0.000175 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08221 on 16 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.5957,	Adjusted R-squared:  0.5704 
# F-statistic: 23.58 on 1 and 16 DF,  p-value: 0.0001754


#It might be more efficient to add color to the other plot to illustrate which effects had been previously observed in the DLPFC
#For the sake of visual comparison, I plotted using the same y-axis limits as for Schizophrenia (BP has smaller effects in the qPCR dataset)

pdf("Plot_BP_LogFC_qPCR_vs_Gandal_RNASeqMeta_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC", col="gray50", ylim=c(-0.6, 0.6))
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_BD.P.Value<0.05,], col=1, ylim=c(-0.6, 0.6))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_BD.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.247480 -0.049599  0.005109  0.058883  0.168962 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.045260   0.008501   5.324 7.01e-07 ***
#   Gandal_RNASeqMeta_BD.logFC 0.536928   0.099003   5.423 4.61e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.08187 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2403,	Adjusted R-squared:  0.2321 
# F-statistic: 29.41 on 1 and 93 DF,  p-value: 4.606e-07

#Comparing different diagnoses:

pdf("Plot_BP_LogFC_qPCR_vs_SCHIZ_Gandal_RNASeqMeta_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP vs. SCHIZ", ylab="Frontal Pole qPCR: BP Log2FC", xlab="Gandal RNA-Seq Meta-Analysis: SCHIZ Log2FC", col="gray50", ylim=c(-0.6, 0.6))
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_SCZ.P.Value<0.05,], col=1, ylim=c(-0.6, 0.6))
abline(TempLine, col=2, lwd=3)
dev.off()

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_SCZ.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.249466 -0.032257  0.007692  0.052588  0.187639 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.053712   0.009135   5.880 6.38e-08 ***
#   Gandal_RNASeqMeta_SCZ.logFC 0.490145   0.099153   4.943 3.39e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.08359 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2081,	Adjusted R-squared:  0.1996 
# F-statistic: 24.44 on 1 and 93 DF,  p-value: 3.389e-06

pdf("Plot_SCHIZ_LogFC_qPCR_vs_BP_Gandal_RNASeqMeta_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ vs. BP", ylab="Frontal Pole qPCR: Schiz Log2FC", xlab="Gandal RNA-Seq Meta-Analysis: BP Log2FC", col="gray50", ylim=c(-0.6, 0.6))
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_RNASeqMeta_BD.P.Value<0.05,], col=1, ylim=c(-0.6, 0.6))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_BD.logFC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.23882 -0.03972 -0.00478  0.05255  0.36565 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.018049   0.009476   1.905   0.0599 .  
# Gandal_RNASeqMeta_BD.logFC 0.729346   0.110354   6.609 2.39e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09126 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.3196,	Adjusted R-squared:  0.3123 
# F-statistic: 43.68 on 1 and 93 DF,  p-value: 2.394e-09

#Ignore: It turns out that the Gandal Microarray Meta-analysis for Schiz actually includes frontal pole tissue so this comparison is circular.
# I might output it anyway just for the sake of parallel structure:

# #Comparison with the Gandal Microarray Meta-Analysis:
# 
pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_MicroarrayMeta.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()


pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_MicroarrayMeta_wGandalP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC", col="gray50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_SCZ.P.value<0.05,], col=1, ylim=c(-0.6, 0.6))
dev.off()

# 
# summary.lm(TempLine)
# # Call:
# #   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_SCZ.beta_log2FC, 
# #      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.97105 -0.02472  0.02761  0.07970  0.41917 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)  
# # (Intercept)                           -0.005619   0.014812  -0.379   0.7052  
# # Gandal_MicroarrayMeta_SCZ.beta_log2FC  0.335926   0.143769   2.337   0.0214 *
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 0.1518 on 103 degrees of freedom
# # (6 observations deleted due to missingness)
# # Multiple R-squared:  0.05034,	Adjusted R-squared:  0.04112 
# # F-statistic:  5.46 on 1 and 103 DF,  p-value: 0.0214
# 
# 
# #Filtering by p-value:
# 
# pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_MicroarrayMeta_P05.pdf", height=5.5, width=5)
# plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_SCZ.P.value<0.05,], main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC")
# TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_SCZ.P.value<0.05,])
# abline(TempLine, col=2, lwd=3)
# dev.off()
# 
# summary.lm(TempLine)
# 
# # Call:
# #   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_SCZ.beta_log2FC, 
# #      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_SCZ.P.value < 
# #                                                                0.05, ])
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.32388 -0.03970  0.01841  0.06435  0.22188 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)   
# # (Intercept)                           -0.02997    0.01553  -1.929  0.06017 . 
# # Gandal_MicroarrayMeta_SCZ.beta_log2FC  0.30916    0.10153   3.045  0.00392 **
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 0.1053 on 44 degrees of freedom
# # (6 observations deleted due to missingness)
# # Multiple R-squared:  0.174,	Adjusted R-squared:  0.1553 
# # F-statistic: 9.271 on 1 and 44 DF,  p-value: 0.003922
# 
# 
pdf("Plot_BP_LogFC_qPCR_vs_Gandal_MicroarrayMeta.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()


pdf("Plot_BP_LogFC_qPCR_vs_Gandal_MicroarrayMeta_wGandalP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC", col="gray50")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
points(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_BD.P.value<0.05,], col=1, ylim=c(-0.6, 0.6))
dev.off()

#
# summary.lm(TempLine)
# # Call:
# #   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_BD.beta_log2FC,
# #      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# #
# # Residuals:
# #   Min       1Q   Median       3Q      Max
# # -0.70451 -0.05159  0.02711  0.08346  0.23162
# #
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)
# # (Intercept)                          0.023764   0.012949   1.835   0.0694 .
# # Gandal_MicroarrayMeta_BD.beta_log2FC 0.005151   0.123012   0.042   0.9667
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# #
# # Residual standard error: 0.1311 on 103 degrees of freedom
# # (6 observations deleted due to missingness)
# # Multiple R-squared:  1.702e-05,	Adjusted R-squared:  -0.009692
# # F-statistic: 0.001753 on 1 and 103 DF,  p-value: 0.9667

sqrt(1.702e-05)
#[1] 0.00412553
#
# #Filtered by p-value:
#
# pdf("Plot_BP_LogFC_qPCR_vs_Gandal_MicroarrayMeta_P05.pdf", height=5.5, width=5)
# plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_BD.P.value<0.05,], main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC")
# TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_BD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_BD.P.value<0.05,])
# abline(TempLine, col=2, lwd=3)
# dev.off()
#
# summary.lm(TempLine)
#
# # Call:
# #   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_BD.beta_log2FC,
# #      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_BD.P.value <
# #                                                                0.05, ])
# #
# # Residuals:
# #   Min       1Q   Median       3Q      Max
# # -0.45577 -0.03609  0.01874  0.08210  0.15216
# #
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)
# # (Intercept)                          -0.01103    0.02138  -0.516    0.609
# # Gandal_MicroarrayMeta_BD.beta_log2FC  0.01859    0.12812   0.145    0.885
# #
# # Residual standard error: 0.1283 on 35 degrees of freedom
# # (6 observations deleted due to missingness)
# # Multiple R-squared:  0.0006015,	Adjusted R-squared:  -0.02795
# # F-statistic: 0.02106 on 1 and 35 DF,  p-value: 0.8854




pdf("Plot_BP_LogFC_qPCR_vs_Gandal_MicroarrayMeta_MDD_wGandalP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP vs. MDD", ylab="Frontal Pole qPCR: BP Log2FC", xlab="Gandal Microarray Meta-Analysis: MDD Log2FC", col="gray50")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
points(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_MDD.P.value<0.05,], col=1, ylim=c(-0.6, 0.6))
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_MDD.beta_log2FC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.67687 -0.04283  0.01680  0.08026  0.24550 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                            0.02809    0.01230   2.284  0.02432 * 
#   Gandal_MicroarrayMeta_MDD.beta_log2FC  0.69445    0.24481   2.837  0.00545 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.126 on 107 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.06994,	Adjusted R-squared:  0.06125 
# F-statistic: 8.047 on 1 and 107 DF,  p-value: 0.005453


pdf("Plot_SCHIZ_LogFC_qPCR_vs_Gandal_MicroarrayMeta_MDD_wGandalP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ vs. MDD", ylab="Frontal Pole qPCR: Schiz Log2FC", xlab="Gandal Microarray Meta-Analysis: MDD Log2FC", col="gray50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gandal_MicroarrayMeta_MDD.P.value<0.05,], col=1, ylim=c(-0.6, 0.6))
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_MDD.beta_log2FC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.93215 -0.04472  0.02413  0.07221  0.46565 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.002084   0.014066   0.148 0.882472    
# Gandal_MicroarrayMeta_MDD.beta_log2FC 1.074401   0.280068   3.836 0.000212 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1442 on 107 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.1209,	Adjusted R-squared:  0.1127 
# F-statistic: 14.72 on 1 and 107 DF,  p-value: 0.0002116


#Here's a version of AAD to match our other Figure:

pdf("GandalAAD_Vs_PritzkerDiagnosis_allTargets.pdf", width=5, height=5.5)
plot(Gandal_MicroarrayMeta_AAD.beta_log2FC~Diagnosis_BP_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, col=3, pch=16, ylab="Gandal Microarray Meta-Analysis: AAD Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Gandal_MicroarrayMeta_AAD.beta_log2FC~Diagnosis_BP_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
abline(BP_Line, col=3, lwd=3)

points(Gandal_MicroarrayMeta_AAD.beta_log2FC~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, col=4, pch=16)
Schiz_Line<-lm(Gandal_MicroarrayMeta_AAD.beta_log2FC~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()



#Comparison with Lanz grey matter microarray:

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Lanz_Microarray.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Lanz: Microarray Log2FC")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_Schiz_logFC_M2, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94942 -0.04124  0.02701  0.07826  0.43447 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                          -0.009695   0.014464   -0.67    0.504
# FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_Schiz_logFC_M2  0.164038   0.108657    1.51    0.134
# 
# Residual standard error: 0.1515 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.02067,	Adjusted R-squared:  0.0116 
# F-statistic: 2.279 on 1 and 108 DF,  p-value: 0.134

sqrt(0.02067)
#[1] 0.1437707

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Lanz_Microarray_wLanzP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Lanz_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Lanz: Microarray Log2FC", col="grey50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Lanz_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Lanz_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_Schiz_P.Value_M2<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()




pdf("Plot_BP_LogFC_qPCR_vs_Lanz_Microarray.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_BP_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Lanz: Microarray Log2FC")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_BP_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_BP_logFC_M2, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.69496 -0.04205  0.02947  0.07966  0.22204 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                                        0.02118    0.01237   1.712   0.0898 .
# FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_BP_logFC_M2  0.08600    0.11569   0.743   0.4589  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1297 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.00509,	Adjusted R-squared:  -0.004122 
# F-statistic: 0.5526 on 1 and 108 DF,  p-value: 0.4589

sqrt(0.00509)
#0.07134424


pdf("Plot_BP_LogFC_qPCR_vs_Lanz_Microarray_wLanzP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Lanz_BP_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Lanz: Microarray Log2FC", col="grey50")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Lanz_BP_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_BP_PostHocSummary_MLM_Beta~Lanz_BP_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_BP_P.Value_M2<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()




#Comparison with Narayan grey matter microarray:

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Narayan_Microarray.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Narayan: Microarray Log2FC")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_logFC_M2, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.97999 -0.05590  0.02657  0.07456  0.41694 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                             -0.006051   0.014503  -0.417    0.677
# FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_logFC_M2  0.167709   0.101833   1.647    0.102
# 
# Residual standard error: 0.1512 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.0245,	Adjusted R-squared:  0.01547 
# F-statistic: 2.712 on 1 and 108 DF,  p-value: 0.1025

sqrt(0.0245)
#[1] 0.1565248

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Narayan_Microarray_wNarayanP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Narayan_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Narayan: Microarray Log2FC", col="grey50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Narayan_Schiz_logFC_M2, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_P.Value_M2<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()




#Comparison with an average of Lanz and Narayan grey matter microarray:

FrontalPolePCRResults_vsDLPFCResults_forGraphing$LanzNarayan_Schiz_AVELog2FC<-apply(cbind(FrontalPolePCRResults_vsDLPFCResults_forGraphing$Lanz_Schiz_logFC_M2,FrontalPolePCRResults_vsDLPFCResults_forGraphing$Narayan_Schiz_logFC_M2), 1, function(y) median(y, na.rm=TRUE))


pdf("Plot_SCHIZ_LogFC_qPCR_vs_LanzNarayan_Microarray.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$LanzNarayan_Schiz_AVELog2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Grey Matter: Microarray Log2FC")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$LanzNarayan_Schiz_AVELog2FC, data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ FrontalPolePCRResults_vsDLPFCResults_forGraphing$LanzNarayan_Schiz_AVELog2FC, 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.96380 -0.04781  0.02796  0.07599  0.42321 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                                                  -0.007809   0.014422  -0.541   0.5893  
# FrontalPolePCRResults_vsDLPFCResults_forGraphing$LanzNarayan_Schiz_AVELog2FC  0.184744   0.110719   1.669   0.0981 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1512 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.02513,	Adjusted R-squared:  0.0161 
# F-statistic: 2.784 on 1 and 108 DF,  p-value: 0.0981


##################################

#Charts comparing the results for other variables:

#Age

pdf("Plot_Age_LogFC_qPCR_vs_Pritzker_Microarray.pdf", height=5.5, width=5)
plot(Age_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="Age", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC")
TempLine<-lm(Age_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

pdf("Plot_Age_LogFC_qPCR_vs_Pritzker_Microarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Age_PostHocSummary_MLM_Beta~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="Age", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Age_PostHocSummary_MLM_Beta~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Age_PostHocSummary_MLM_Beta~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Pval_Nominal_Age_Centered<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Age_PostHocSummary_MLM_Beta ~ Pritzker_Beta_Age_Centered..per.increase.of.1.year., 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0120105 -0.0029901  0.0001503  0.0029699  0.0130543 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                         -0.0035895  0.0004685  -7.661 1.10e-11 ***
#   Pritzker_Beta_Age_Centered..per.increase.of.1.year.  0.6107912  0.0780446   7.826 4.85e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.00466 on 102 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.3752,	Adjusted R-squared:  0.3691 
# F-statistic: 61.25 on 1 and 102 DF,  p-value: 4.852e-12



pdf("Plot_BrainPH_LogFC_qPCR_vs_Pritzker_Microarray.pdf", height=5.5, width=5)
plot(FrontalPolePCRResults_vsDLPFCResults_forGraphing$pH_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="Brain pH", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC")
TempLine<-lm(FrontalPolePCRResults_vsDLPFCResults_forGraphing$pH_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

pdf("Plot_BrainPH_LogFC_qPCR_vs_Pritzker_Microarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(pH_PostHocSummary_MLM_Beta~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="Brain pH", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(pH_PostHocSummary_MLM_Beta~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(pH_PostHocSummary_MLM_Beta~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Pval_Nominal_BrainPH_Centered<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = pH_PostHocSummary_MLM_Beta ~ Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.52013 -0.14041 -0.01175  0.17461  1.84649 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                              0.07677    0.03867   1.985   0.0498 *  
#   Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit.  0.80600    0.11330   7.114  1.6e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3943 on 102 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.3316,	Adjusted R-squared:  0.3251 
# F-statistic: 50.61 on 1 and 102 DF,  p-value: 1.602e-10


pdf("Plot_PMI_LogFC_qPCR_vs_Pritzker_Microarray.pdf", height=5.5, width=5)
plot(FrontalPolePCRResults_vsDLPFCResults_forGraphing$PMI_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="PMI", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC")
TempLine<-lm(FrontalPolePCRResults_vsDLPFCResults_forGraphing$PMI_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_PMI..per.increase.of.1.hour.,  data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()


pdf("Plot_PMI_LogFC_qPCR_vs_Pritzker_Microarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(PMI_PostHocSummary_MLM_Beta~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="PMI", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(PMI_PostHocSummary_MLM_Beta~Pritzker_Beta_PMI..per.increase.of.1.hour.,  data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(PMI_PostHocSummary_MLM_Beta~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Pval_Nominal_PMI<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary(TempLine)

# Call:
#   lm(formula = PMI_PostHocSummary_MLM_Beta ~ Pritzker_Beta_PMI..per.increase.of.1.hour., 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.018383 -0.004136  0.000043  0.004146  0.031752 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                0.0033482  0.0007403   4.523 1.65e-05 ***
#   Pritzker_Beta_PMI..per.increase.of.1.hour. 1.0652399  0.1334826   7.980 2.25e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.007471 on 102 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.3844,	Adjusted R-squared:  0.3783 
# F-statistic: 63.69 on 1 and 102 DF,  p-value: 2.253e-12


pdf("Plot_Gender_LogFC_qPCR_vs_Pritzker_Microarray.pdf", height=5.5, width=5)
plot(FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gender_F_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="Gender", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC")
TempLine<-lm(FrontalPolePCRResults_vsDLPFCResults_forGraphing$Gender_F_PostHocSummary_MLM_Beta~FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
abline(TempLine, col=2, lwd=3)
dev.off()

pdf("Plot_Gender_LogFC_qPCR_vs_Pritzker_Microarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Gender_F_PostHocSummary_MLM_Beta~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing, main="Gender", ylab="Frontal Pole: qPCR Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Gender_F_PostHocSummary_MLM_Beta~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing)
points(Gender_F_PostHocSummary_MLM_Beta~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPolePCRResults_vsDLPFCResults_forGraphing[FrontalPolePCRResults_vsDLPFCResults_forGraphing$Pritzker_Pval_Nominal_Gender_F<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Gender_F_PostHocSummary_MLM_Beta ~ Pritzker_Beta_Gender..Female.vs..a.Male.baseline., 
#      data = FrontalPolePCRResults_vsDLPFCResults_forGraphing)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41632 -0.08952  0.01113  0.07731  0.49075 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                       -0.04709    0.01431  -3.290  0.00138 **
#   Pritzker_Beta_Gender..Female.vs..a.Male.baseline.  0.07492    0.21941   0.341  0.73344   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1457 on 102 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.001142,	Adjusted R-squared:  -0.008651 
# F-statistic: 0.1166 on 1 and 102 DF,  p-value: 0.7334


##################################

##I should add in here:
#1) Graphs comparing results with frontal pole microarray meta-analysis output using the same format (to match for paper)

setwd("~/Documents/Microarray Gen/FrontalPole/MaycoxAndIwamoto")

#Model 3 meta-analysis results:
FrontalPole_Schiz_MetaAnalysis<-read.csv("RMAOutputForIwamotoVsMaycoxSchizResults_correctedAgain.csv", header=TRUE, stringsAsFactors = FALSE)

str(FrontalPole_Schiz_MetaAnalysis)
# 'data.frame':	8398 obs. of  7 variables:
#   $ X           : int  2864 6860 6996 6396 5283 4134 138 2035 5901 5736 ...
# $ b           : num  0.166 -0.149 -0.124 -0.228 -0.253 ...
# $ se          : num  0.0389 0.035 0.0301 0.0559 0.0621 ...
# $ pval        : num  2.13e-05 2.16e-05 3.86e-05 4.58e-05 4.63e-05 ...
# $ BH_adjPval  : num  0.0778 0.0778 0.0778 0.0778 0.0778 ...
# $ EntrezGeneID: chr  "27101" "7923" "81605" "7075" ...
# $ GeneSymbol  : chr  "CACYBP" "HSD17B8" "URM1" "TIE1" ...

#Results from all model output:
FrontalPole_Microarray_AllResults<-read.csv("IwamotoVsMaycox_AllResults.csv", header=TRUE, stringsAsFactors = FALSE)

str(FrontalPole_Microarray_AllResults)
#Let's cut that down to the content that is most useful - Model 3 output (the other model output was confounded by gender and contained a lot of sex chromosome genes)
colnames(FrontalPole_Microarray_AllResults)

FrontalPole_Microarray_AllResults<-FrontalPole_Microarray_AllResults[,c(104:107,57:103,150:192)]

str(FrontalPole_Microarray_AllResults)
'data.frame':	20575 obs. of  94 variables:
  
colnames(FrontalPole_Microarray_AllResults)
# [1] "X"                                                          
# [2] "ProbesetID"                                                 
# [3] "EntrezGeneID"                                               
# [4] "GeneSymbol"                                                 
# [5] "Iwamoto_Model3_A"                                           
# [6] "Iwamoto_Model3_Coef..Intercept."                            
# [7] "Iwamoto_Model3_Coef.DiagnosisFactorBipolar"                 
# [8] "Iwamoto_Model3_Coef.DiagnosisFactorDepression"              
# [9] "Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia"           
# [10] "Iwamoto_Model3_Coef.BrainpH"                                
# [11] "Iwamoto_Model3_Coef.RNADegradPerSample"                     
# [12] "Iwamoto_Model3_Coef.RateofDeathFactorOn.respirator"         
# [13] "Iwamoto_Model3_Coef.RateofDeathFactorPossible.anoxia"       
# [14] "Iwamoto_Model3_Coef.Age"                                    
# [15] "Iwamoto_Model3_Coef.PMI"                                    
# [16] "Iwamoto_Model3_Coef.GenderFactorF"                          
# [17] "Iwamoto_Model3_t..Intercept."                               
# [18] "Iwamoto_Model3_t.DiagnosisFactorBipolar"                    
# [19] "Iwamoto_Model3_t.DiagnosisFactorDepression"                 
# [20] "Iwamoto_Model3_t.DiagnosisFactorSchizophrenia"              
# [21] "Iwamoto_Model3_t.BrainpH"                                   
# [22] "Iwamoto_Model3_t.RNADegradPerSample"                        
# [23] "Iwamoto_Model3_t.RateofDeathFactorOn.respirator"            
# [24] "Iwamoto_Model3_t.RateofDeathFactorPossible.anoxia"          
# [25] "Iwamoto_Model3_t.Age"                                       
# [26] "Iwamoto_Model3_t.PMI"                                       
# [27] "Iwamoto_Model3_t.GenderFactorF"                             
# [28] "Iwamoto_Model3_p.value..Intercept."                         
# [29] "Iwamoto_Model3_p.value.DiagnosisFactorBipolar"              
# [30] "Iwamoto_Model3_p.value.DiagnosisFactorDepression"           
# [31] "Iwamoto_Model3_p.value.DiagnosisFactorSchizophrenia"        
# [32] "Iwamoto_Model3_p.value.BrainpH"                             
# [33] "Iwamoto_Model3_p.value.RNADegradPerSample"                  
# [34] "Iwamoto_Model3_p.value.RateofDeathFactorOn.respirator"      
# [35] "Iwamoto_Model3_p.value.RateofDeathFactorPossible.anoxia"    
# [36] "Iwamoto_Model3_p.value.Age"                                 
# [37] "Iwamoto_Model3_p.value.PMI"                                 
# [38] "Iwamoto_Model3_p.value.GenderFactorF"                       
# [39] "Iwamoto_Model3_p.value.adj..Intercept."                     
# [40] "Iwamoto_Model3_p.value.adj.DiagnosisFactorBipolar"          
# [41] "Iwamoto_Model3_p.value.adj.DiagnosisFactorDepression"       
# [42] "Iwamoto_Model3_p.value.adj.DiagnosisFactorSchizophrenia"    
# [43] "Iwamoto_Model3_p.value.adj.BrainpH"                         
# [44] "Iwamoto_Model3_p.value.adj.RNADegradPerSample"              
# [45] "Iwamoto_Model3_p.value.adj.RateofDeathFactorOn.respirator"  
# [46] "Iwamoto_Model3_p.value.adj.RateofDeathFactorPossible.anoxia"
# [47] "Iwamoto_Model3_p.value.adj.Age"                             
# [48] "Iwamoto_Model3_p.value.adj.PMI"                             
# [49] "Iwamoto_Model3_p.value.adj.GenderFactorF"                   
# [50] "Iwamoto_Model3_F"                                           
# [51] "Iwamoto_Model3_F.p.value"                                   
# [52] "Maycox_Model3_A"                                            
# [53] "Maycox_Model3_Coef..Intercept."                             
# [54] "Maycox_Model3_Coef.DiagnosisFactorScz"                      
# [55] "Maycox_Model3_Coef.BrainpH"                                 
# [56] "Maycox_Model3_Coef.Age"                                     
# [57] "Maycox_Model3_Coef.RNADegradPerSample"                      
# [58] "Maycox_Model3_Coef.PMI"                                     
# [59] "Maycox_Model3_Coef.GenderFactorFemale"                      
# [60] "Maycox_Model3_Coef.GenderFactorNot.Available"               
# [61] "Maycox_Model3_Coef.ScanDateDayOnly01.21.04"                 
# [62] "Maycox_Model3_Coef.ScanDateDayOnly01.29.04"                 
# [63] "Maycox_Model3_t..Intercept."                                
# [64] "Maycox_Model3_t.DiagnosisFactorScz"                         
# [65] "Maycox_Model3_t.BrainpH"                                    
# [66] "Maycox_Model3_t.Age"                                        
# [67] "Maycox_Model3_t.RNADegradPerSample"                         
# [68] "Maycox_Model3_t.PMI"                                        
# [69] "Maycox_Model3_t.GenderFactorFemale"                         
# [70] "Maycox_Model3_t.GenderFactorNot.Available"                  
# [71] "Maycox_Model3_t.ScanDateDayOnly01.21.04"                    
# [72] "Maycox_Model3_t.ScanDateDayOnly01.29.04"                    
# [73] "Maycox_Model3_p.value..Intercept."                          
# [74] "Maycox_Model3_p.value.DiagnosisFactorScz"                   
# [75] "Maycox_Model3_p.value.BrainpH"                              
# [76] "Maycox_Model3_p.value.Age"                                  
# [77] "Maycox_Model3_p.value.RNADegradPerSample"                   
# [78] "Maycox_Model3_p.value.PMI"                                  
# [79] "Maycox_Model3_p.value.GenderFactorFemale"                   
# [80] "Maycox_Model3_p.value.GenderFactorNot.Available"            
# [81] "Maycox_Model3_p.value.ScanDateDayOnly01.21.04"              
# [82] "Maycox_Model3_p.value.ScanDateDayOnly01.29.04"              
# [83] "Maycox_Model3_p.value.adj..Intercept."                      
# [84] "Maycox_Model3_p.value.adj.DiagnosisFactorScz"               
# [85] "Maycox_Model3_p.value.adj.BrainpH"                          
# [86] "Maycox_Model3_p.value.adj.Age"                              
# [87] "Maycox_Model3_p.value.adj.RNADegradPerSample"               
# [88] "Maycox_Model3_p.value.adj.PMI"                              
# [89] "Maycox_Model3_p.value.adj.GenderFactorFemale"               
# [90] "Maycox_Model3_p.value.adj.GenderFactorNot.Available"        
# [91] "Maycox_Model3_p.value.adj.ScanDateDayOnly01.21.04"          
# [92] "Maycox_Model3_p.value.adj.ScanDateDayOnly01.29.04"          
# [93] "Maycox_Model3_F"                                            
# [94] "Maycox_Model3_F.p.value"    

#Which annotation should I use for joining?

sum(is.na(FrontalPole_Microarray_AllResults$EntrezGeneID))
#[1] 0

sum(is.na(FrontalPole_Microarray_AllResults$GeneSymbol))
#[1] 94

sum(is.na(FrontalPole_Schiz_MetaAnalysis$EntrezGeneID))
#[1] 0

sum(is.na(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$ENTREZID))
#[1] 9541

sum(is.na(DLPFC_LanzNarayanPrizker_wSymbols$ENTREZID))
#[1] 0

#So the NAs shouldn't cause trouble for the Entrez ID because only the RNA-Seq data includes them.

colnames(FrontalPole_Microarray_AllResults)[3]<-"ENTREZID"
colnames(FrontalPole_Schiz_MetaAnalysis)[6]<-"ENTREZID"

FrontalPole_Microarray_vsDLPFC<-join_all(dfs=list(FrontalPole_Microarray_AllResults, FrontalPole_Schiz_MetaAnalysis, Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols, DLPFC_LanzNarayanPrizker_wSymbols), by="ENTREZID", type="left", match="all")
str(FrontalPole_Microarray_vsDLPFC)
#'data.frame':	20875 obs. of  199 variables:


FrontalPole_qPCRGenes_MasterAnnotation_DiagnosisEffects<-join_all(dfs=list(FrontalPole_qPCRGenes_MasterAnnotation, CellTypeSpecific_DiagnosisEffects, Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols, DLPFC_LanzNarayanPrizker_wSymbols), by="SYMBOL", type="left", match="all")

str(FrontalPole_Microarray_vsDLPFC)
#'data.frame':	20875 obs. of  199 variables:

#So 300 duplicated mappings were added during the joining process

sum(duplicated(FrontalPole_Microarray_vsDLPFC$ENTREZID))
#[1] 300
  
FrontalPole_Microarray_vsDLPFC$GeneSymbol[duplicated(FrontalPole_Microarray_vsDLPFC$ENTREZID)]
#It's a mixture of genes - no NA, only a few genes have 5+ entries, so not likely to dominate analyses but we should double check.

write.csv(FrontalPole_Microarray_vsDLPFC, "FrontalPole_Microarray_vsDLPFC.csv")

colnames(FrontalPolePCRResults_vsDLPFCResults_forGraphing)

colnames(FrontalPole_Microarray_vsDLPFC)[4]<-"Gene.Symbol"

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC<-join(FrontalPolePCRResults_vsDLPFCResults_forGraphing[,c(1:98)], FrontalPole_Microarray_vsDLPFC, by="Gene.Symbol")

dim(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
#[1] 111 293

write.csv(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, "FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC.csv")

##########################


pdf("Plot_SCHIZ_LogFC_qPCR_vs_FrontalPoleMicroarrayMeta_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~b, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Frontal Pole: Microarray Meta-Analysis Log2FC", col="gray50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~b, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~b, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$pval<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ b, data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94519 -0.05279  0.02872  0.07819  0.45594 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -0.01650    0.01614  -1.022   0.3093  
# b            0.36937    0.15269   2.419   0.0175 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1564 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.0592,	Adjusted R-squared:  0.04908 
# F-statistic: 5.852 on 1 and 93 DF,  p-value: 0.01751


pdf("Plot_SCHIZ_LogFC_qPCR_vs_Iwamoto_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Iwamoto: Microarray Log2FC", col="gray50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_p.value.DiagnosisFactorSchizophrenia<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94746 -0.06219  0.02288  0.07806  0.45810 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                      -0.00757    0.01618  -0.468   0.6409  
# Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia  0.32883    0.13638   2.411   0.0179 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1565 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.05883,	Adjusted R-squared:  0.04871 
# F-statistic: 5.813 on 1 and 93 DF,  p-value: 0.01787

pdf("Plot_SCHIZ_LogFC_qPCR_vs_Maycox_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.DiagnosisFactorScz, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Frontal Pole: qPCR Log2FC", xlab="Maycox: Microarray Log2FC", col="gray50")
TempLine<-lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.DiagnosisFactorScz, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.DiagnosisFactorScz, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Maycox_Model3_p.value.DiagnosisFactorScz<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Maycox_Model3_Coef.DiagnosisFactorScz, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.96171 -0.05292  0.02559  0.07713  0.45903 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                           -0.01234    0.01504  -0.820    0.414
# Maycox_Model3_Coef.DiagnosisFactorScz  0.09317    0.09687   0.962    0.338
# 
# Residual standard error: 0.1525 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.008493,	Adjusted R-squared:  -0.0006872 
# F-statistic: 0.9251 on 1 and 108 DF,  p-value: 0.3383



pdf("Plot_BP_LogFC_qPCR_vs_Iwamoto_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Diagnosis_BP_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorBipolar, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="BP", ylab="Frontal Pole: qPCR Log2FC", xlab="Iwamoto: Microarray Log2FC", col="gray50")
TempLine<-lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorBipolar, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Diagnosis_BP_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorBipolar, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_p.value.DiagnosisFactorBipolar<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.DiagnosisFactorBipolar, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.69850 -0.05066  0.03005  0.08516  0.20702 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                 0.01861    0.01408   1.322    0.190
# Iwamoto_Model3_Coef.DiagnosisFactorBipolar  0.18466    0.12756   1.448    0.151
# 
# Residual standard error: 0.1362 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.02204,	Adjusted R-squared:  0.01152 
# F-statistic: 2.096 on 1 and 93 DF,  p-value: 0.1511

#Cross-diagnosis correlation stats:
summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.67678 -0.04518  0.01553  0.07255  0.24542 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                       0.02102    0.01363   1.542  0.12644   
# Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia  0.33496    0.11490   2.915  0.00445 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1318 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.08373,	Adjusted R-squared:  0.07388 
# F-statistic: 8.498 on 1 and 93 DF,  p-value: 0.004454

summary.lm(lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorBipolar, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.DiagnosisFactorBipolar, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.96920 -0.05877  0.03264  0.08074  0.41187 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                -0.009282   0.016462  -0.564    0.574
# Iwamoto_Model3_Coef.DiagnosisFactorBipolar  0.229838   0.149178   1.541    0.127
# 
# Residual standard error: 0.1593 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.02489,	Adjusted R-squared:  0.0144 
# F-statistic: 2.374 on 1 and 93 DF,  p-value: 0.1268

summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~b, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ b, data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.67033 -0.05389  0.02233  0.06701  0.24390 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.01114    0.01334   0.835 0.405659    
# b            0.44583    0.12622   3.532 0.000643 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1293 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.1183,	Adjusted R-squared:  0.1088 
# F-statistic: 12.48 on 1 and 93 DF,  p-value: 0.0006435

summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.DiagnosisFactorScz, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Maycox_Model3_Coef.DiagnosisFactorScz, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.68588 -0.04292  0.01706  0.06501  0.25271 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                            0.01488    0.01259   1.181   0.2401  
# Maycox_Model3_Coef.DiagnosisFactorScz  0.16440    0.08111   2.027   0.0451 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1277 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.03665,	Adjusted R-squared:  0.02773 
# F-statistic: 4.109 on 1 and 108 DF,  p-value: 0.04513


summary.lm(lm(Lanz_Schiz_logFC_M2~Lanz_BP_logFC_M2, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Lanz_Schiz_logFC_M2 ~ Lanz_BP_logFC_M2, data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.260642 -0.048265 -0.000186  0.045416  0.212775 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.003662   0.007303   0.501    0.617    
# Lanz_BP_logFC_M2 1.021284   0.068282  14.957   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.07657 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.6744,	Adjusted R-squared:  0.6714 
# F-statistic: 223.7 on 1 and 108 DF,  p-value: < 2.2e-16

sqrt(0.6744)
#[1] 0.8212186

#Here's a version to match our other Figure:

pdf("IwamotoMDD_Vs_PritzkerDiagnosis_allTargets.pdf", width=5, height=5.5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Diagnosis_BP_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, col=3, pch=16, ylab="Iwamoto Microarray: MDD Log2FC", xlab="Pritzker qPCR: Diagnosis Log2FC")
BP_Line<-lm(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Diagnosis_BP_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
abline(BP_Line, col=3, lwd=3)

points(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, col=4, pch=16)
Schiz_Line<-lm(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
abline(Schiz_Line, col=4, lwd=3)
dev.off()


pdf("Plot_Age_LogFC_qPCR_vs_Iwamoto_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Age_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.Age, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="Age", ylab="Frontal Pole: qPCR Log2FC", xlab="Iwamoto: Microarray Log2FC", col="gray50")
TempLine<-lm(Age_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.Age, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Age_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.Age, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_p.value.Age<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Age_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.Age, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0155569 -0.0031658  0.0004118  0.0032944  0.0091418 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -0.0036454  0.0005514  -6.611 2.38e-09 ***
#   Iwamoto_Model3_Coef.Age  0.6275774  0.1170682   5.361 6.00e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.00511 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2361,	Adjusted R-squared:  0.2278 
# F-statistic: 28.74 on 1 and 93 DF,  p-value: 6.004e-07

pdf("Plot_pH_LogFC_qPCR_vs_Iwamoto_wMetaP05Filled.pdf", height=5.5, width=5)
plot(pH_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.BrainpH, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="pH", ylab="Frontal Pole: qPCR Log2FC", xlab="Iwamoto: Microarray Log2FC", col="gray50")
TempLine<-lm(pH_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.BrainpH, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(pH_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.BrainpH, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_p.value.BrainpH<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = pH_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.BrainpH, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.52138 -0.17022  0.03257  0.22558  1.81088 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  0.11216    0.04443   2.524   0.0133 *  
#   Iwamoto_Model3_Coef.BrainpH  0.46974    0.08123   5.783 9.78e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4281 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2645,	Adjusted R-squared:  0.2566 
# F-statistic: 33.44 on 1 and 93 DF,  p-value: 9.777e-08

pdf("Plot_RIN_LogFC_qPCR_vs_Iwamoto_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Integrity_RIN_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.RNADegradPerSample, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="RNA Integrity", ylab="Frontal Pole: qPCR Log2FC", xlab="Iwamoto: Microarray Log2FC", col="gray50")
TempLine<-lm(Integrity_RIN_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.RNADegradPerSample, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Integrity_RIN_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.RNADegradPerSample, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_p.value.RNADegradPerSample<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Integrity_RIN_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.RNADegradPerSample, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41990 -0.08105 -0.00872  0.06910  0.55734 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                            0.006513   0.015417   0.422    0.674
# Iwamoto_Model3_Coef.RNADegradPerSample 0.017960   0.011129   1.614    0.110
# 
# Residual standard error: 0.1472 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.02724,	Adjusted R-squared:  0.01678 
# F-statistic: 2.604 on 1 and 93 DF,  p-value: 0.11


pdf("Plot_Age_LogFC_qPCR_vs_Maycox_wMaycoxP05Filled.pdf", height=5.5, width=5)
plot(Age_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.Age, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="Age", ylab="Frontal Pole: qPCR Log2FC", xlab="Maycox: Microarray Log2FC", col="gray50")
TempLine<-lm(Age_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.Age, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Age_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.Age, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Maycox_Model3_p.value.Age<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Age_PostHocSummary_MLM_Beta ~ Maycox_Model3_Coef.Age, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0116310 -0.0025524  0.0003132  0.0032897  0.0105909 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -0.0032978  0.0004573  -7.212 7.97e-11 ***
#   Maycox_Model3_Coef.Age  0.8363294  0.1018115   8.214 5.04e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.004556 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3845,	Adjusted R-squared:  0.3788 
# F-statistic: 67.48 on 1 and 108 DF,  p-value: 5.041e-13

pdf("Plot_pH_LogFC_qPCR_vs_Maycox_wMetaP05Filled.pdf", height=5.5, width=5)
plot(pH_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.BrainpH, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="pH", ylab="Frontal Pole: qPCR Log2FC", xlab="Maycox: Microarray Log2FC", col="gray50")
TempLine<-lm(pH_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.BrainpH, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(pH_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.BrainpH, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Maycox_Model3_p.value.BrainpH<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = pH_PostHocSummary_MLM_Beta ~ Maycox_Model3_Coef.BrainpH, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.52119 -0.15903  0.01555  0.19045  1.97472 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.03179    0.04221   0.753    0.453    
# Maycox_Model3_Coef.BrainpH  0.52945    0.10661   4.966 2.57e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4268 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.1859,	Adjusted R-squared:  0.1784 
# F-statistic: 24.66 on 1 and 108 DF,  p-value: 2.569e-06


pdf("Plot_RNADeg_LogFC_qPCR_vs_Maycox_wMetaP05Filled.pdf", height=5.5, width=5)
plot(Integrity_RIN_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.RNADegradPerSample, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, main="RNA Integrity", ylab="Frontal Pole: qPCR Log2FC", xlab="Maycox: Microarray Log2FC", col="gray50")
TempLine<-lm(Integrity_RIN_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.RNADegradPerSample, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
points(Integrity_RIN_PostHocSummary_MLM_Beta~Maycox_Model3_Coef.RNADegradPerSample, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Maycox_Model3_p.value.RNADegradPerSample<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Integrity_RIN_PostHocSummary_MLM_Beta ~ Maycox_Model3_Coef.RNADegradPerSample, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38609 -0.08081 -0.00560  0.06249  0.57957 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.005169   0.013549   0.382 0.703553    
# Maycox_Model3_Coef.RNADegradPerSample 0.018952   0.005036   3.764 0.000273 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1371 on 108 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.1159,	Adjusted R-squared:  0.1078 
# F-statistic: 14.16 on 1 and 108 DF,  p-value: 0.0002728


#Conclusion:
#Some of the Schiz effects replicate in qPCR vs. microarray meta-analysis
#The effects of variables with larger effect sizes (age, pH) replicate better.


#For simplicity sake, I wonder if I could output a correlation matrix:
#(did I already do this?)

FrontalPole_qPCRGenes_Log2FC_Diagnosis<-cbind.data.frame(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Diagnosis_BP_PostHocSummary_MLM_Beta,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Diagnosis_Schiz_PostHocSummary_MLM_Beta,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_Coef.DiagnosisFactorBipolar,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_Coef.DiagnosisFactorDepression,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Maycox_Model3_Coef.DiagnosisFactorScz,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$b,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.logFC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_BD.logFC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_MicroarrayMeta_SCZ.beta_log2FC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_MicroarrayMeta_BD.beta_log2FC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_MicroarrayMeta_MDD.beta_log2FC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_MicroarrayMeta_AAD.beta_log2FC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_MicroarrayMeta_ASD.beta_log2FC,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Lanz_BP_logFC_M2,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Lanz_Schiz_logFC_M2,
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Narayan_Schiz_logFC_M2)

CorMatrix_FrontalPole_qPCRGenes_Log2FC_Diagnosis<-cor(as.matrix(FrontalPole_qPCRGenes_Log2FC_Diagnosis), use="pairwise.complete.obs")
colnames(CorMatrix_FrontalPole_qPCRGenes_Log2FC_Diagnosis)<-colnames(FrontalPole_qPCRGenes_Log2FC_Diagnosis)
row.names(CorMatrix_FrontalPole_qPCRGenes_Log2FC_Diagnosis)<-colnames(FrontalPole_qPCRGenes_Log2FC_Diagnosis)

write.csv(CorMatrix_FrontalPole_qPCRGenes_Log2FC_Diagnosis, "CorMatrix_FrontalPole_qPCRGenes_Log2FC_Diagnosis.csv")

#I need to output some stats to go with that:
summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Diagnosis_Schiz_PostHocSummary_MLM_Beta, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Diagnosis_Schiz_PostHocSummary_MLM_Beta, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.250664 -0.042162  0.007351  0.036918  0.214745 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                             0.027284   0.007014    3.89 0.000173 ***
#   Diagnosis_Schiz_PostHocSummary_MLM_Beta 0.698617   0.046322   15.08  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.07376 on 109 degrees of freedom
# Multiple R-squared:  0.676,	Adjusted R-squared:  0.6731 
# F-statistic: 227.5 on 1 and 109 DF,  p-value: < 2.2e-16

summary.lm(lm(Iwamoto_Model3_Coef.DiagnosisFactorBipolar~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorBipolar ~ Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.21724 -0.05395  0.00268  0.05460  0.34099 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      -0.007612   0.010329  -0.737    0.463    
# Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia  0.400730   0.087091   4.601 1.32e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09991 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.1854,	Adjusted R-squared:  0.1767 
# F-statistic: 21.17 on 1 and 93 DF,  p-value: 1.325e-05

summary.lm(lm(Iwamoto_Model3_Coef.DiagnosisFactorBipolar~Iwamoto_Model3_Coef.DiagnosisFactorDepression, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorBipolar ~ Iwamoto_Model3_Coef.DiagnosisFactorDepression, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.31416 -0.03498  0.01531  0.05618  0.25022 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   -0.009897   0.009695  -1.021     0.31    
# Iwamoto_Model3_Coef.DiagnosisFactorDepression  0.555582   0.093741   5.927 5.19e-08 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09431 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2742,	Adjusted R-squared:  0.2664 
# F-statistic: 35.13 on 1 and 93 DF,  p-value: 5.192e-08

summary.lm(lm(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Iwamoto_Model3_Coef.DiagnosisFactorDepression, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia ~ 
#        Iwamoto_Model3_Coef.DiagnosisFactorDepression, data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.31171 -0.05618  0.01785  0.06148  0.27581 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   -0.01132    0.01092  -1.037    0.303    
# Iwamoto_Model3_Coef.DiagnosisFactorDepression  0.51393    0.10555   4.869 4.57e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1062 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2032,	Adjusted R-squared:  0.1946 
# F-statistic: 23.71 on 1 and 93 DF,  p-value: 4.573e-06

summary.lm(lm(Gandal_RNASeqMeta_BD.logFC~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Gandal_RNASeqMeta_BD.logFC ~ Gandal_RNASeqMeta_SCZ.logFC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.185384 -0.025977  0.006596  0.028609  0.105777 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.011673   0.005623   2.076   0.0407 *  
#   Gandal_RNASeqMeta_SCZ.logFC 0.784705   0.061039  12.856   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.05146 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.6399,	Adjusted R-squared:  0.636 
# F-statistic: 165.3 on 1 and 93 DF,  p-value: < 2.2e-16

summary.lm(lm(Gandal_MicroarrayMeta_BD.beta_log2FC~Gandal_MicroarrayMeta_SCZ.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Gandal_MicroarrayMeta_BD.beta_log2FC ~ Gandal_MicroarrayMeta_SCZ.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.241543 -0.025757  0.005723  0.028029  0.089209 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           -0.014530   0.004539  -3.201  0.00182 ** 
#   Gandal_MicroarrayMeta_SCZ.beta_log2FC  0.905440   0.044057  20.551  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04651 on 103 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.8039,	Adjusted R-squared:  0.802 
# F-statistic: 422.4 on 1 and 103 DF,  p-value: < 2.2e-16

summary.lm(lm(Gandal_MicroarrayMeta_BD.beta_log2FC~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Gandal_MicroarrayMeta_BD.beta_log2FC ~ Gandal_MicroarrayMeta_MDD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.20020 -0.06030 -0.00141  0.02988  0.41072 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           -0.007389   0.009584  -0.771    0.443    
# Gandal_MicroarrayMeta_MDD.beta_log2FC  0.871881   0.195255   4.465 2.05e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09614 on 103 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.1622,	Adjusted R-squared:  0.1541 
# F-statistic: 19.94 on 1 and 103 DF,  p-value: 2.052e-05

summary.lm(lm(Gandal_MicroarrayMeta_SCZ.beta_log2FC~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Gandal_MicroarrayMeta_SCZ.beta_log2FC ~ Gandal_MicroarrayMeta_MDD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.16374 -0.04883 -0.01773  0.02207  0.44751 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.007934   0.009253   0.857    0.393    
# Gandal_MicroarrayMeta_MDD.beta_log2FC 0.967603   0.188504   5.133 1.35e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09281 on 103 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.2037,	Adjusted R-squared:  0.196 
# F-statistic: 26.35 on 1 and 103 DF,  p-value: 1.351e-06

summary.lm(lm(Gandal_MicroarrayMeta_BD.beta_log2FC~Gandal_MicroarrayMeta_ASD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Gandal_MicroarrayMeta_BD.beta_log2FC ~ Gandal_MicroarrayMeta_ASD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.223053 -0.042194 -0.001518  0.037982  0.213986 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.002401   0.007800   0.308    0.759    
# Gandal_MicroarrayMeta_ASD.beta_log2FC 0.388433   0.042995   9.034 1.63e-14 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.07375 on 97 degrees of freedom
# (12 observations deleted due to missingness)
# Multiple R-squared:  0.4569,	Adjusted R-squared:  0.4514 
# F-statistic: 81.62 on 1 and 97 DF,  p-value: 1.63e-14

summary.lm(lm(Gandal_MicroarrayMeta_SCZ.beta_log2FC~Gandal_MicroarrayMeta_ASD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Gandal_MicroarrayMeta_SCZ.beta_log2FC ~ Gandal_MicroarrayMeta_ASD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.257240 -0.037728 -0.006516  0.039829  0.239748 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.018026   0.007432   2.425   0.0171 *  
#   Gandal_MicroarrayMeta_ASD.beta_log2FC 0.403580   0.040969   9.851 2.81e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.07028 on 97 degrees of freedom
# (12 observations deleted due to missingness)
# Multiple R-squared:  0.5001,	Adjusted R-squared:  0.495 
# F-statistic: 97.04 on 1 and 97 DF,  p-value: 2.813e-16

#Cross diagnosis:
summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_SCZ.logFC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.249466 -0.032257  0.007692  0.052588  0.187639 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.053712   0.009135   5.880 6.38e-08 ***
#   Gandal_RNASeqMeta_SCZ.logFC 0.490145   0.099153   4.943 3.39e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.08359 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.2081,	Adjusted R-squared:  0.1996 
# F-statistic: 24.44 on 1 and 93 DF,  p-value: 3.389e-06

summary.lm(lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_RNASeqMeta_BD.logFC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_RNASeqMeta_BD.logFC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.23882 -0.03972 -0.00478  0.05255  0.36565 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.018049   0.009476   1.905   0.0599 .  
# Gandal_RNASeqMeta_BD.logFC 0.729346   0.110354   6.609 2.39e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.09126 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.3196,	Adjusted R-squared:  0.3123 
# F-statistic: 43.68 on 1 and 93 DF,  p-value: 2.394e-09

summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorDepression, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.DiagnosisFactorDepression, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.63825 -0.05538  0.02235  0.08478  0.22059 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                    0.01879    0.01344   1.398  0.16543   
# Iwamoto_Model3_Coef.DiagnosisFactorDepression  0.41480    0.12994   3.192  0.00193 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1307 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.09876,	Adjusted R-squared:  0.08907 
# F-statistic: 10.19 on 1 and 93 DF,  p-value: 0.001926

summary.lm(lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Iwamoto_Model3_Coef.DiagnosisFactorDepression, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Iwamoto_Model3_Coef.DiagnosisFactorDepression, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.88656 -0.06425  0.01767  0.08000  0.31363 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   -0.008708   0.015409  -0.565 0.573363    
# Iwamoto_Model3_Coef.DiagnosisFactorDepression  0.570320   0.148989   3.828 0.000234 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1499 on 93 degrees of freedom
# (16 observations deleted due to missingness)
# Multiple R-squared:  0.1361,	Adjusted R-squared:  0.1268 
# F-statistic: 14.65 on 1 and 93 DF,  p-value: 0.0002341

summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_MDD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.67687 -0.04283  0.01680  0.08026  0.24550 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                            0.02809    0.01230   2.284  0.02432 * 
#   Gandal_MicroarrayMeta_MDD.beta_log2FC  0.69445    0.24481   2.837  0.00545 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.126 on 107 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.06994,	Adjusted R-squared:  0.06125 
# F-statistic: 8.047 on 1 and 107 DF,  p-value: 0.005453

summary.lm(lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_MDD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.93215 -0.04472  0.02413  0.07221  0.46565 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.002084   0.014066   0.148 0.882472    
# Gandal_MicroarrayMeta_MDD.beta_log2FC 1.074401   0.280068   3.836 0.000212 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1442 on 107 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.1209,	Adjusted R-squared:  0.1127 
# F-statistic: 14.72 on 1 and 107 DF,  p-value: 0.0002116

summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_ASD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_ASD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.70741 -0.03962  0.02780  0.08554  0.19042 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                            0.02694    0.01389   1.940   0.0551 .
# Gandal_MicroarrayMeta_ASD.beta_log2FC  0.09913    0.07726   1.283   0.2025  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1332 on 100 degrees of freedom
# (9 observations deleted due to missingness)
# Multiple R-squared:  0.01619,	Adjusted R-squared:  0.006355 
# F-statistic: 1.646 on 1 and 100 DF,  p-value: 0.2025

summary.lm(lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_ASD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_ASD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.98181 -0.04134  0.02725  0.07404  0.35432 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                           0.002901   0.016021   0.181   0.8567  
# Gandal_MicroarrayMeta_ASD.beta_log2FC 0.203672   0.089148   2.285   0.0244 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1537 on 100 degrees of freedom
# (9 observations deleted due to missingness)
# Multiple R-squared:  0.04961,	Adjusted R-squared:  0.0401 
# F-statistic:  5.22 on 1 and 100 DF,  p-value: 0.02445

summary.lm(lm(Diagnosis_BP_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_AAD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))
# Call:
#   lm(formula = Diagnosis_BP_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_AAD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.68748 -0.04005  0.02089  0.07967  0.18906 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                            0.02509    0.01296   1.936   0.0556 .
# Gandal_MicroarrayMeta_AAD.beta_log2FC  0.10911    0.05803   1.880   0.0629 .
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.131 on 102 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.0335,	Adjusted R-squared:  0.02402 
# F-statistic: 3.535 on 1 and 102 DF,  p-value: 0.06292

summary.lm(lm(Diagnosis_Schiz_PostHocSummary_MLM_Beta~Gandal_MicroarrayMeta_AAD.beta_log2FC, data=FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC))

# Call:
#   lm(formula = Diagnosis_Schiz_PostHocSummary_MLM_Beta ~ Gandal_MicroarrayMeta_AAD.beta_log2FC, 
#      data = FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94740 -0.04092  0.02255  0.08141  0.32599 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                           -0.003772   0.015032  -0.251   0.8024  
# Gandal_MicroarrayMeta_AAD.beta_log2FC  0.168471   0.067301   2.503   0.0139 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1519 on 102 degrees of freedom
# (7 observations deleted due to missingness)
# Multiple R-squared:  0.05788,	Adjusted R-squared:  0.04864 
# F-statistic: 6.266 on 1 and 102 DF,  p-value: 0.01389

##########################

#Making plots of frontal pole microarray vs. DLPFC:

pdf("Plot_SCHIZ_FrontalPole_MicroarrayMetaAnalysis_vs_RNASeqMetaAnalysis_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(b~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Frontal Pole: Microarray Meta-Analysis Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC", col="grey50")
TempLine<-lm(b~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC)
points(b~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.P.Value<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()
#Wow. That's pretty abyssmal. I guess the frontal pole meta-analysis is underpowered and based on old microarray data, but still. Ouch.

summary.lm(TempLine)

# Call:
#   lm(formula = b ~ Gandal_RNASeqMeta_SCZ.logFC, data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.39252 -0.04652 -0.00320  0.04492  0.83974 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.0070372  0.0009756   7.213 6.07e-13 ***
#   Gandal_RNASeqMeta_SCZ.logFC 0.0194573  0.0125894   1.546    0.122    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07981 on 6701 degrees of freedom
# (14172 observations deleted due to missingness)
# Multiple R-squared:  0.0003563,	Adjusted R-squared:  0.0002072 
# F-statistic: 2.389 on 1 and 6701 DF,  p-value: 0.1223


summary.lm(lm(b~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.P.Value<0.05,]))

# Call:
#   lm(formula = b ~ Gandal_RNASeqMeta_SCZ.logFC, data = FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.P.Value < 
#                                                                                         0.05, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.34865 -0.04960 -0.00380  0.04374  0.53871 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.005654   0.001910   2.960  0.00312 ** 
#   Gandal_RNASeqMeta_SCZ.logFC 0.074585   0.014989   4.976 7.15e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07923 on 1720 degrees of freedom
# (9035 observations deleted due to missingness)
# Multiple R-squared:  0.01419,	Adjusted R-squared:  0.01362 
# F-statistic: 24.76 on 1 and 1720 DF,  p-value: 7.145e-07

#If we filter by DLPFC p-value, there is a slight positive relationship between the two datasets.


#Plotting for the individual datasets.

#This one had a weird y-axis limits (lots of blank space) so I set it artificially

pdf("Plot_SCHIZ_FrontalPole_IwamotoMicroarray_vs_RNASeqMetaAnalysis_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.P.Value<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia ~ 
#        Gandal_RNASeqMeta_SCZ.logFC, data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.79365 -0.05254  0.00827  0.06227  0.53048 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 -0.006841   0.001313  -5.208 1.96e-07 ***
#   Gandal_RNASeqMeta_SCZ.logFC  0.050775   0.016879   3.008  0.00264 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1079 on 6761 degrees of freedom
# (14112 observations deleted due to missingness)
# Multiple R-squared:  0.001337,	Adjusted R-squared:  0.001189 
# F-statistic: 9.049 on 1 and 6761 DF,  p-value: 0.002638

#Also had a weird ylim, so I set it artificially:

pdf("Plot_BP_FrontalPole_IwamotoMicroarray_vs_RNASeqMetaAnalysis_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorBipolar~Gandal_RNASeqMeta_BD.logFC, data=FrontalPole_Microarray_vsDLPFC, main="BP", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Iwamoto_Model3_Coef.DiagnosisFactorBipolar~Gandal_RNASeqMeta_BD.logFC, data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.DiagnosisFactorBipolar~Gandal_RNASeqMeta_BD.logFC, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_BD.P.Value<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorBipolar ~ Gandal_RNASeqMeta_BD.logFC, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.59575 -0.05230  0.00129  0.05501  0.74305 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                -0.003730   0.001119  -3.334  0.00086 ***
#   Gandal_RNASeqMeta_BD.logFC  0.156248   0.017350   9.005  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09114 on 6761 degrees of freedom
# (14112 observations deleted due to missingness)
# Multiple R-squared:  0.01185,	Adjusted R-squared:  0.01171 
# F-statistic:  81.1 on 1 and 6761 DF,  p-value: < 2.2e-16

pdf("Plot_SCHIZ_FrontalPole_MaycoxMicroarray_vs_RNASeqMetaAnalysis_wRNASeqP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.DiagnosisFactorScz~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Gandal: RNA-Seq Meta-Analysis Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Maycox_Model3_Coef.DiagnosisFactorScz~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.DiagnosisFactorScz~Gandal_RNASeqMeta_SCZ.logFC, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.P.Value<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.DiagnosisFactorScz ~ Gandal_RNASeqMeta_SCZ.logFC, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.73467 -0.06284 -0.00852  0.05556  0.91106 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  0.019013   0.000883  21.533   <2e-16 ***
#   Gandal_RNASeqMeta_SCZ.logFC -0.027316   0.012083  -2.261   0.0238 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1023 on 13446 degrees of freedom
# (7427 observations deleted due to missingness)
# Multiple R-squared:  0.00038,	Adjusted R-squared:  0.0003056 
# F-statistic: 5.111 on 1 and 13446 DF,  p-value: 0.02379


#I wonder if there would be a stronger correlation with DLPFC microarray experiments:

pdf("Plot_SCHIZ_FrontalPole_IwamotoMicroarray_vs_NarayanMicroarray_wNarayanP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Narayan_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Narayan: Microarray Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Narayan_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Narayan_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Narayan_Schiz_P.Value_M2<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia ~ 
#        Narayan_Schiz_logFC_M2, data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.81513 -0.04948  0.00732  0.05784  0.56202 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -0.001264   0.001106  -1.143    0.253    
# Narayan_Schiz_logFC_M2 -0.133348   0.013311 -10.018   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1019 on 8485 degrees of freedom
# (12388 observations deleted due to missingness)
# Multiple R-squared:  0.01169,	Adjusted R-squared:  0.01157 
# F-statistic: 100.4 on 1 and 8485 DF,  p-value: < 2.2e-16

#Negative correlation!


pdf("Plot_SCHIZ_FrontalPole_IwamotoMicroarray_vs_LanzMicroarray_wNarayanP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Lanz_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Lanz: Microarray Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Lanz_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia~Lanz_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Lanz_Schiz_P.Value_M2<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia ~ 
#        Lanz_Schiz_logFC_M2, data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.80193 -0.04906  0.00819  0.05873  0.49503 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -0.001296   0.001112  -1.165    0.244    
# Lanz_Schiz_logFC_M2  0.053961   0.012102   4.459 8.34e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1024 on 8485 degrees of freedom
# (12388 observations deleted due to missingness)
# Multiple R-squared:  0.002338,	Adjusted R-squared:  0.00222 
# F-statistic: 19.88 on 1 and 8485 DF,  p-value: 8.345e-06

#Positive correlation

pdf("Plot_SCHIZ_FrontalPole_MaycoxMicroarray_vs_NarayanMicroarray_wNarayanP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.DiagnosisFactorScz~Narayan_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Narayan: Microarray Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Maycox_Model3_Coef.DiagnosisFactorScz~Narayan_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.DiagnosisFactorScz~Narayan_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Narayan_Schiz_P.Value_M2<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.DiagnosisFactorScz ~ Narayan_Schiz_logFC_M2, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.64807 -0.05840 -0.01115  0.05048  1.16087 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -0.004754   0.000671  -7.085 1.44e-12 ***
#   Narayan_Schiz_logFC_M2 -0.346561   0.009202 -37.660  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09399 on 19634 degrees of freedom
# (1239 observations deleted due to missingness)
# Multiple R-squared:  0.06737,	Adjusted R-squared:  0.06732 
# F-statistic:  1418 on 1 and 19634 DF,  p-value: < 2.2e-16

#That's a pretty strong *negative* correlation


pdf("Plot_SCHIZ_FrontalPole_MaycoxMicroarray_vs_LanzMicroarray_wNarayanP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.DiagnosisFactorScz~Lanz_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Lanz: Microarray Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Maycox_Model3_Coef.DiagnosisFactorScz~Lanz_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.DiagnosisFactorScz~Lanz_Schiz_logFC_M2, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Lanz_Schiz_P.Value_M2<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.DiagnosisFactorScz ~ Lanz_Schiz_logFC_M2, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.70826 -0.05893 -0.01201  0.05083  0.90823 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -0.0041009  0.0006953  -5.898 3.73e-09 ***
#   Lanz_Schiz_logFC_M2  0.0046784  0.0080812   0.579    0.563    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09733 on 19634 degrees of freedom
# (1239 observations deleted due to missingness)
# Multiple R-squared:  1.707e-05,	Adjusted R-squared:  -3.386e-05 
# F-statistic: 0.3352 on 1 and 19634 DF,  p-value: 0.5626

#Also pretty underwhelming. I think these individual datasets may just be too weakly-powered.


#Other variables:

pdf("Plot_MDD_FrontalPole_IwamotoMicroarray_vs_GandalMicroarrayMeta_wGandalP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_Microarray_vsDLPFC, main="MDD", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Gandal: Microarray Meta-Analysis Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Gandal_MicroarrayMeta_MDD.beta_log2FC, data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Gandal_MicroarrayMeta_MDD.P.value<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorDepression ~ 
#        Gandal_MicroarrayMeta_MDD.beta_log2FC, data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.43545 -0.04591 -0.00404  0.04155  0.58290 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                           -0.0007863  0.0008431  -0.933    0.351
# Gandal_MicroarrayMeta_MDD.beta_log2FC  0.0246786  0.0200303   1.232    0.218
# 
# Residual standard error: 0.07424 on 8195 degrees of freedom
# (12678 observations deleted due to missingness)
# Multiple R-squared:  0.0001852,	Adjusted R-squared:  6.319e-05 
# F-statistic: 1.518 on 1 and 8195 DF,  p-value: 0.218

#Not much relationship. 


pdf("Plot_MDD_FrontalPole_IwamotoMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Pritzker_Beta_Diagnosis..MDD.vs..Control.baseline., data=FrontalPole_Microarray_vsDLPFC, main="MDD", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50", ylim=c(-1,1))
TempLine<-lm(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Pritzker_Beta_Diagnosis..MDD.vs..Control.baseline., data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.DiagnosisFactorDepression~Pritzker_Beta_Diagnosis..MDD.vs..Control.baseline., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_Diagnosis_MD<0.05,], col=1, ylim=c(-1,1))
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.DiagnosisFactorDepression ~ 
#        Pritzker_Beta_Diagnosis..MDD.vs..Control.baseline., data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.46678 -0.04536 -0.00425  0.04225  0.57029 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                        -0.0008414  0.0008215  -1.024    0.306    
# Pritzker_Beta_Diagnosis..MDD.vs..Control.baseline. -0.2134216  0.0182736 -11.679   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07455 on 8234 degrees of freedom
# (12639 observations deleted due to missingness)
# Multiple R-squared:  0.0163,	Adjusted R-squared:  0.01618 
# F-statistic: 136.4 on 1 and 8234 DF,  p-value: < 2.2e-16

#That's a pretty strong *negative* relationship.  Good thing we aren't hanging our hat on those results - I would have no idea what to make of them.

#Hmm... what about other variables that typically have stronger signals in microarray data?  

pdf("Plot_Age_FrontalPole_IwamotoMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.Age~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPole_Microarray_vsDLPFC, main="Age", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Iwamoto_Model3_Coef.Age~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.Age~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_Age_Centered<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.Age ~ Pritzker_Beta_Age_Centered..per.increase.of.1.year., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0162493 -0.0015329  0.0002184  0.0016580  0.0292297 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                         8.680e-06  3.151e-05   0.275    0.783    
# Pritzker_Beta_Age_Centered..per.increase.of.1.year. 3.968e-01  1.347e-02  29.465   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.00286 on 8234 degrees of freedom
# (12639 observations deleted due to missingness)
# Multiple R-squared:  0.09538,	Adjusted R-squared:  0.09527 
# F-statistic: 868.2 on 1 and 8234 DF,  p-value: < 2.2e-16

pdf("Plot_pH_FrontalPole_IwamotoMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.BrainpH~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPole_Microarray_vsDLPFC, main="pH", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Iwamoto_Model3_Coef.BrainpH~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.BrainpH~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_BrainPH_Centered<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.BrainpH ~ Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.69354 -0.14218  0.02242  0.16276  1.45133 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                             0.008391   0.002924    2.87  0.00411 ** 
#   Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit. 1.100465   0.018492   59.51  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2653 on 8234 degrees of freedom
# (12639 observations deleted due to missingness)
# Multiple R-squared:  0.3007,	Adjusted R-squared:  0.3007 
# F-statistic:  3541 on 1 and 8234 DF,  p-value: < 2.2e-16


pdf("Plot_PMI_FrontalPole_IwamotoMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.PMI~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPole_Microarray_vsDLPFC, main="PMI", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Iwamoto_Model3_Coef.PMI~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.PMI~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_PMI<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.PMI ~ Pritzker_Beta_PMI..per.increase.of.1.hour., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.012826 -0.002245 -0.000410  0.001931  0.033754 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                1.173e-04  3.565e-05   3.291    0.001 ** 
#   Pritzker_Beta_PMI..per.increase.of.1.hour. 6.380e-02  9.872e-03   6.463 1.09e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.003232 on 8234 degrees of freedom
# (12639 observations deleted due to missingness)
# Multiple R-squared:  0.005047,	Adjusted R-squared:  0.004926 
# F-statistic: 41.77 on 1 and 8234 DF,  p-value: 1.085e-10

#Less correlation - maybe because Iwamoto includes RIN and Pritzker does not?


pdf("Plot_Gender_FrontalPole_IwamotoMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Iwamoto_Model3_Coef.GenderFactorF~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPole_Microarray_vsDLPFC, main="Gender", ylab="Iwamoto: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Iwamoto_Model3_Coef.GenderFactorF~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPole_Microarray_vsDLPFC)
points(Iwamoto_Model3_Coef.GenderFactorF~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_Gender_F<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

#LOL - that's the sex chromosome genes. Not particularly informative.

summary.lm(TempLine)

# Call:
#   lm(formula = Iwamoto_Model3_Coef.GenderFactorF ~ Pritzker_Beta_Gender..Female.vs..a.Male.baseline., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.89698 -0.04051  0.00515  0.04752  0.34691 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                       -0.0030607  0.0009525  -3.213  0.00132 ** 
#   Pritzker_Beta_Gender..Female.vs..a.Male.baseline.  0.9113272  0.0143476  63.518  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08644 on 8234 degrees of freedom
# (12639 observations deleted due to missingness)
# Multiple R-squared:  0.3288,	Adjusted R-squared:  0.3288 
# F-statistic:  4034 on 1 and 8234 DF,  p-value: < 2.2e-16


##I should add in here:
#1) Graphs comparing results with frontal pole microarray meta-analysis output using the same format (to match for paper)
#2) Comparison with drug treatment or stress results???  Are there any studies with great enough power?


pdf("Plot_Age_FrontalPole_MaycoxMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.Age~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPole_Microarray_vsDLPFC, main="Age", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Maycox_Model3_Coef.Age~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.Age~Pritzker_Beta_Age_Centered..per.increase.of.1.year., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_Age_Centered<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)


# Call:
#   lm(formula = Maycox_Model3_Coef.Age ~ Pritzker_Beta_Age_Centered..per.increase.of.1.year., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0226551 -0.0014246 -0.0002056  0.0012395  0.0276365 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                         1.947e-04  2.287e-05   8.514   <2e-16 ***
#   Pritzker_Beta_Age_Centered..per.increase.of.1.year. 6.867e-01  1.035e-02  66.378   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.002497 on 11918 degrees of freedom
# (8955 observations deleted due to missingness)
# Multiple R-squared:  0.2699,	Adjusted R-squared:  0.2699 
# F-statistic:  4406 on 1 and 11918 DF,  p-value: < 2.2e-16

pdf("Plot_pH_FrontalPole_MaycoxMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.BrainpH~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPole_Microarray_vsDLPFC, main="pH", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Maycox_Model3_Coef.BrainpH~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.BrainpH~Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_BrainPH_Centered<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.BrainpH ~ Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.29250 -0.11231 -0.01585  0.10449  1.45349 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                             0.012617   0.001643   7.681  1.7e-14 ***
#   Pritzker_Beta_BrainPH_Centered..per.increase.of.1.unit. 1.060547   0.010865  97.612  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1793 on 11918 degrees of freedom
# (8955 observations deleted due to missingness)
# Multiple R-squared:  0.4443,	Adjusted R-squared:  0.4442 
# F-statistic:  9528 on 1 and 11918 DF,  p-value: < 2.2e-16

pdf("Plot_PMI_FrontalPole_MaycoxMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.PMI~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPole_Microarray_vsDLPFC, main="PMI", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Maycox_Model3_Coef.PMI~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.PMI~Pritzker_Beta_PMI..per.increase.of.1.hour., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_PMI<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.PMI ~ Pritzker_Beta_PMI..per.increase.of.1.hour., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.062038 -0.002940  0.000561  0.003447  0.034697 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                -2.407e-04  5.237e-05  -4.595 4.37e-06 ***
#   Pritzker_Beta_PMI..per.increase.of.1.hour.  1.739e-01  1.540e-02  11.290  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.005718 on 11918 degrees of freedom
# (8955 observations deleted due to missingness)
# Multiple R-squared:  0.01058,	Adjusted R-squared:  0.0105 
# F-statistic: 127.5 on 1 and 11918 DF,  p-value: < 2.2e-16

pdf("Plot_Gender_FrontalPole_MaycoxMicroarray_vs_PritzkerMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.GenderFactorFemale~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPole_Microarray_vsDLPFC, main="Gender", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Pritzker: Microarray Log2FC", col="grey50")
TempLine<-lm(Maycox_Model3_Coef.GenderFactorFemale~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPole_Microarray_vsDLPFC)
points(Maycox_Model3_Coef.GenderFactorFemale~Pritzker_Beta_Gender..Female.vs..a.Male.baseline., data=FrontalPole_Microarray_vsDLPFC[FrontalPole_Microarray_vsDLPFC$Pritzker_Pval_Nominal_Gender_F<0.05,], col=1)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.GenderFactorFemale ~ Pritzker_Beta_Gender..Female.vs..a.Male.baseline., 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.08232 -0.03358 -0.00140  0.03279  2.88086 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                       0.0042795  0.0007266    5.89 3.96e-09 ***
#   Pritzker_Beta_Gender..Female.vs..a.Male.baseline. 0.6928116  0.0122040   56.77  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07932 on 11918 degrees of freedom
# (8955 observations deleted due to missingness)
# Multiple R-squared:  0.2129,	Adjusted R-squared:  0.2128 
# F-statistic:  3223 on 1 and 11918 DF,  p-value: < 2.2e-16



####Comparing Frontal Pole with Frontal Pole:

#Has a weird x-axis - artificially set it at -1 to 1
pdf("Plot_SCHIZ_FrontalPole_MaycoxMicroarray_vs_IwamotoMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.DiagnosisFactorScz~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_Microarray_vsDLPFC, main="SCHIZ", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Iwamoto: Frontal Pole Microarray Log2FC", col=1, xlim=c(-1,1))
TempLine<-lm(Maycox_Model3_Coef.DiagnosisFactorScz~Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, data=FrontalPole_Microarray_vsDLPFC)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.DiagnosisFactorScz ~ Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.73157 -0.06336 -0.00974  0.05884  0.87839 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      0.004702   0.001126   4.177 2.98e-05 ***
#   Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia 0.105668   0.010986   9.618  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.104 on 8535 degrees of freedom
# (12338 observations deleted due to missingness)
# Multiple R-squared:  0.01072,	Adjusted R-squared:  0.01061 
# F-statistic: 92.51 on 1 and 8535 DF,  p-value: < 2.2e-16


pdf("Plot_Age_FrontalPole_MaycoxMicroarray_vs_IwamotoMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.Age~Iwamoto_Model3_Coef.Age, data=FrontalPole_Microarray_vsDLPFC, main="Age", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Iwamoto: Frontal Pole Microarray Log2FC", col=1)
TempLine<-lm(Maycox_Model3_Coef.Age~Iwamoto_Model3_Coef.Age, data=FrontalPole_Microarray_vsDLPFC)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.Age ~ Iwamoto_Model3_Coef.Age, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0239107 -0.0016758 -0.0002746  0.0013611  0.0315767 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             2.208e-04  3.218e-05   6.859 7.41e-12 ***
#   Iwamoto_Model3_Coef.Age 2.118e-01  1.071e-02  19.779  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.002974 on 8535 degrees of freedom
# (12338 observations deleted due to missingness)
# Multiple R-squared:  0.04383,	Adjusted R-squared:  0.04371 
# F-statistic: 391.2 on 1 and 8535 DF,  p-value: < 2.2e-16


pdf("Plot_pH_FrontalPole_MaycoxMicroarray_vs_IwamotoMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.BrainpH~Iwamoto_Model3_Coef.BrainpH, data=FrontalPole_Microarray_vsDLPFC, main="pH", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Iwamoto: Frontal Pole Microarray Log2FC", col=1)
TempLine<-lm(Maycox_Model3_Coef.BrainpH~Iwamoto_Model3_Coef.BrainpH, data=FrontalPole_Microarray_vsDLPFC)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.BrainpH ~ Iwamoto_Model3_Coef.BrainpH, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.6217 -0.1502 -0.0244  0.1298  1.3931 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.013227   0.002490   5.312 1.11e-07 ***
#   Iwamoto_Model3_Coef.BrainpH 0.304216   0.007905  38.485  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.23 on 8535 degrees of freedom
# (12338 observations deleted due to missingness)
# Multiple R-squared:  0.1479,	Adjusted R-squared:  0.1478 
# F-statistic:  1481 on 1 and 8535 DF,  p-value: < 2.2e-16


pdf("Plot_PMI_FrontalPole_MaycoxMicroarray_vs_IwamotoMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.PMI~Iwamoto_Model3_Coef.PMI, data=FrontalPole_Microarray_vsDLPFC, main="PMI", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Iwamoto: Frontal Pole Microarray Log2FC", col=1)
TempLine<-lm(Maycox_Model3_Coef.PMI~Iwamoto_Model3_Coef.PMI, data=FrontalPole_Microarray_vsDLPFC)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.PMI ~ Iwamoto_Model3_Coef.PMI, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.062915 -0.002692  0.000486  0.003349  0.028758 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -3.658e-04  6.152e-05  -5.947 2.84e-09 ***
#   Iwamoto_Model3_Coef.PMI -4.737e-01  1.902e-02 -24.903  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.005682 on 8535 degrees of freedom
# (12338 observations deleted due to missingness)
# Multiple R-squared:  0.06774,	Adjusted R-squared:  0.06763 
# F-statistic: 620.2 on 1 and 8535 DF,  p-value: < 2.2e-16

#Negative relationship. Huh. 

pdf("Plot_Gender_FrontalPole_MaycoxMicroarray_vs_IwamotoMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.GenderFactorFemale~Iwamoto_Model3_Coef.GenderFactorF, data=FrontalPole_Microarray_vsDLPFC, main="Gender", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Iwamoto: Frontal Pole Microarray Log2FC", col=1)
TempLine<-lm(Maycox_Model3_Coef.GenderFactorFemale~Iwamoto_Model3_Coef.GenderFactorF, data=FrontalPole_Microarray_vsDLPFC)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.GenderFactorFemale ~ Iwamoto_Model3_Coef.GenderFactorF, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.1129 -0.0398 -0.0009  0.0378  3.4848 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       0.0051666  0.0009421   5.484 4.28e-08 ***
#   Iwamoto_Model3_Coef.GenderFactorF 0.4895468  0.0089266  54.841  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08703 on 8535 degrees of freedom
# (12338 observations deleted due to missingness)
# Multiple R-squared:  0.2606,	Adjusted R-squared:  0.2605 
# F-statistic:  3008 on 1 and 8535 DF,  p-value: < 2.2e-16


pdf("Plot_RNADeg_FrontalPole_MaycoxMicroarray_vs_IwamotoMicroarray_wPritzkerP05Filled.pdf", height=5.5, width=5)
plot(Maycox_Model3_Coef.RNADegradPerSample~Iwamoto_Model3_Coef.RNADegradPerSample, data=FrontalPole_Microarray_vsDLPFC, main="RNA Degradation", ylab="Maycox: Frontal Pole Microarray Log2FC", xlab="Iwamoto: Frontal Pole Microarray Log2FC", col=1)
TempLine<-lm(Maycox_Model3_Coef.RNADegradPerSample~Iwamoto_Model3_Coef.RNADegradPerSample, data=FrontalPole_Microarray_vsDLPFC)
abline(TempLine, col=2, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = Maycox_Model3_Coef.RNADegradPerSample ~ Iwamoto_Model3_Coef.RNADegradPerSample, 
#      data = FrontalPole_Microarray_vsDLPFC)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -18.0683  -0.9605   0.0215   0.9074  13.9477 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            -0.02486    0.01962  -1.267    0.205    
# Iwamoto_Model3_Coef.RNADegradPerSample  0.72854    0.01801  40.448   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.812 on 8535 degrees of freedom
# (12338 observations deleted due to missingness)
# Multiple R-squared:  0.1609,	Adjusted R-squared:  0.1608 
# F-statistic:  1636 on 1 and 8535 DF,  p-value: < 2.2e-16


#Conclusion:
#As would be expected, the effects of variables that produce large effect sizes show greater replication across datasets, including the small frontal pole datsets.
#There is only a little bit of replication of SCHIZ signal in the frontal pole datasets.
#There is not strong replication of SCHIZ signal between frontal pole and DLPFC unless you first filter based on which genes had previously shown SCHIZ effects in the DLPFC.

#########################

#After making all of those individual plots, perhaps I should have made some heatmaps first?

setwd("~/Documents/Microarray Gen/FrontalPole/CompareWDLPFC")

SCHIZ_LogFCMatrix<-cbind(FrontalPole_Microarray_vsDLPFC$b,  FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, FrontalPole_Microarray_vsDLPFC$Maycox_Model3_Coef.DiagnosisFactorScz, FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.logFC, FrontalPole_Microarray_vsDLPFC$Lanz_Schiz_logFC_M2, FrontalPole_Microarray_vsDLPFC$Narayan_Schiz_logFC_M2)

colnames(SCHIZ_LogFCMatrix)<-c("FP_MicroarrayMetaAnalysis", "FP_Iwamoto", "FP_Maycox", "Gandal_RNASeqMetaAnalysis", "Lanz_GreyMatter_Microarray", "Narayan_GreyMatter_Microarray")


library(lattice)

pdf("ScatterplotMatrix_SCHIZ_LogFCMatrix.pdf", height=8, width=8)
splom(data.frame(SCHIZ_LogFCMatrix), col=1)
dev.off()

#That's pretty ugly and overwhelming. The names are too long, and it would be nice if it had the correlation matrix heatmap as part of it. hmmm...

heatmap(cor(SCHIZ_LogFCMatrix, use="pairwise.complete"))
#that's sort of nice... except the color values are hard to interpret. And it would be useful to have the qPCR values included too.

dim(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
#[1] 111 293
#This matrix only includes the qPCR genes

qPCR_SCHIZ_LogFCMatrix<-cbind(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Diagnosis_Schiz_PostHocSummary_MLM_Beta, FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$b,  FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Iwamoto_Model3_Coef.DiagnosisFactorSchizophrenia, FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Maycox_Model3_Coef.DiagnosisFactorScz, FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gandal_RNASeqMeta_SCZ.logFC, FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Lanz_Schiz_logFC_M2, FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Narayan_Schiz_logFC_M2)

colnames(qPCR_SCHIZ_LogFCMatrix)<-c("FP_qPCR", "FP_MicroarrayMetaAnalysis", "FP_Iwamoto", "FP_Maycox", "Gandal_RNASeqMetaAnalysis", "Lanz_GreyMatter_Microarray", "Narayan_GreyMatter_Microarray")


heatmap(cor(qPCR_SCHIZ_LogFCMatrix, use="pairwise.complete"))



colnames(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)

FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[,c(42:48, 104:113, 151:156, 192, 198, 202, 214,223,226,249,255,261, 267:271, 273)]

pdf("Heatmap_FrontalPole_qPCRGenes_vs_OtherDatasets_LogFC_Everything.pdf", height=10, width=10)
heatmap(cor(FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything, use="pairwise.complete"), margins=c(10,10))
dev.off()

#Ok, that's pretty overwhelming, LOL.

write.csv(cor(FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything, use="pairwise.complete"), "CorMatrix_FrontalPole_qPCRGenes_vs_OtherDatasets_LogFC_Everything.csv")

cor(cbind(FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything$Age_PostHocSummary_MLM_Beta, FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything$Iwamoto_Model3_Coef.Age, FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything$Maycox_Model3_Coef.Age, FrontalPole_qPCRGenes_vs_FrontalPole_Microarray_vsDLPFC_LogFC_Everything$Pritzker_Beta_Age_Centered..per.increase.of.1.year.), use="pairwise.complete")

#That's actually pretty interesting - the Maycox Frontal Pole dataset definitely stands out as the odd one out, and seems to have a schiz signature that correlates with pH signature. This dataset was confounded by pH - I wonder if our correction for that confound using pH as a covariate in the model wasn't a strong enough correction.


##################################


#I wonder how AverageCq relates to RNA-Seq expression values:

#CMC RNA-Seq methods from Hagenauer et al. 2018
#The dataset was further filtered using an expression threshold (CPM>1 in at least 50 individuals)
#Output from limmaVoom method (log2CPM)

str(GabaGlu_Cq_AllSubjects)
# 'data.frame':	96 obs. of  3 variables:
#   $ X                        : chr  "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# $ GabaGlu_AverageCq_PerGene: num  21.1 25.5 24.1 27.2 23.2 ...
# $ GabaGlu_NumberNA_PerGene : int  0 0 0 0 0 0 0 0 0 0 ...

str(DA5HT_Cq_AllSubjects)

# 'data.frame':	48 obs. of  3 variables:
#   $ X                      : chr  "ADRB1" "ADRB2" "COMT" "DBH" ...
# $ DA5HT_AverageCq_PerGene: num  24.1 24.9 24.5 28.3 33 ...
# $ DA5HT_NumberNA_PerGene : int  0 0 0 0 0 0 0 0 0 0 ...

#This is left over from older code:
# length(c(GabaGlu_AverageCq_PerGene,DA5HT_AverageCq_PerGene))
# DA5HT_AverageCq_NumberNA_PerGene<-cbind(DA5HT_AverageCq_PerGene, DA5HT_NumberNA_PerGene)
# GabaGlu_AverageCq_NumberNA_PerGene<-cbind(GabaGlu_AverageCq_PerGene, GabaGlu_NumberNA_PerGene)

DA5HT_AverageCq_NumberNA_PerGene<-DA5HT_Cq_AllSubjects
colnames(DA5HT_AverageCq_NumberNA_PerGene)<-c("Gene.Symbol", "AverageCq_PerGene", "NumberNA_PerGene")
GabaGlu_AverageCq_NumberNA_PerGene<-GabaGlu_Cq_AllSubjects
colnames(GabaGlu_AverageCq_NumberNA_PerGene)<-c("Gene.Symbol", "AverageCq_PerGene", "NumberNA_PerGene")

AverageCq_NumberNA_PerGene<-rbind(GabaGlu_AverageCq_NumberNA_PerGene, DA5HT_AverageCq_NumberNA_PerGene)

AverageCq_NumberNA_PerGene_DF<-data.frame(AverageCq_NumberNA_PerGene, stringsAsFactors=FALSE)
dim(AverageCq_NumberNA_PerGene_DF)
#[1] 144   3

#Old code:
# colnames(AverageCq_NumberNA_PerGene_DF)[2]<-"AverageCq_PerGene"
# colnames(AverageCq_NumberNA_PerGene_DF)[3]<-"NumberNA_PerGene"

str(Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols)

AverageExpr_RNASeq<-data.frame(Gene.Symbol=Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$SYMBOL, CMC_AveExpr=Gandal_Microarray_RNASeq_W_CMCExpression_wSymbols$CMC_Schiz_AveExpr_M2, stringsAsFactors=FALSE)
str(AverageExpr_RNASeq)
# 'data.frame':	32050 obs. of  2 variables:
#   $ Gene.Symbol: chr  "TSPAN6" "DPM1" "SCYL3" "C1orf112" ...
# $ CMC_AveExpr: num  3.58 4.99 4.13 2.93 1.54 ...

library(plyr)

AverageExpr_FrontalPoleVsRNASeq<-join(AverageCq_NumberNA_PerGene_DF, AverageExpr_RNASeq, by="Gene.Symbol", type="left", match="all")
str(AverageExpr_FrontalPoleVsRNASeq)
# 'data.frame':	144 obs. of  4 variables:
#   $ Gene.Symbol      : chr  "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# $ AverageCq_PerGene: num  21.1 25.5 24.1 27.2 23.2 ...
# $ NumberNA_PerGene : int  0 0 0 0 0 0 0 0 0 0 ...
# $ CMC_AveExpr      : num  8.21 2.93 5.12 NA 6.67 ...

setwd("~/Documents/Microarray Gen/FrontalPole/Output/FinalOutput_AcrossDatasets")
write.csv(AverageExpr_FrontalPoleVsRNASeq, "AverageExpr_FrontalPoleVsRNASeq.csv")

#Oooh - that's pretty cool:

pdf("RNASeq_Log2CPM_vs_qPCR_AverageCq.pdf", width=5, height=5)
plot(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr~AverageExpr_FrontalPoleVsRNASeq$AverageCq_PerGene, xlab="qPCR: Average Cq per Gene", ylab="RNA-Seq: Average log2CPM per Gene", xlim=c(16, 37))
abline(a=0, b=0, col="grey")
dev.off()

sum(is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr))
#[1] 13

AverageExpr_FrontalPoleVsRNASeq$Gene.Symbol[is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr)]
#[1] ADORA2A 18S     GABRR1  GABRR2  DDC     DRD2    DRD3    DRD4    HTR3A   SLC18A1 SLC18A2 SLC6A3  TH

#18S - by definition (since they probably used dT primers for RT)
#So basically 12

#We also were not able to reliably measure 4 of these genes: 18S, GABRR1,  SLC18A1, SLC6A3:
data.frame(AverageExpr_FrontalPoleVsRNASeq$Gene.Symbol[is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr)], AverageExpr_FrontalPoleVsRNASeq$NumberNA_PerGene[is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr)])

table(is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr), AverageExpr_FrontalPoleVsRNASeq$NumberNA_PerGene)
#         0   1   2   4   6  20  24  27  60  68
# FALSE 123   1   4   1   1   0   0   0   1   0
# TRUE    9   0   0   0   0   1   1   1   0   1

#So everything that we couldn't reliably measure, RNA-Seq couldn't reliably measure either, but there were still 9 genes that we were able to measure that RNA-Seq could not.

#Note: 139 is the actual number of transcripts (144 double counts the 5 HK genes included in both datasets)
#139-4=135, the number of transcripts that we were actually able to measure reliably in our dataset
sum(is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr))/139
#[1] 0.09027778
#9% of these transcripts were unmeasurable in RNASeq in the DLPFC (either because of RNA-Seq filtering or because of region or because RNA-Seq is less sensitive than PCR)
(sum(is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr))-4)/(139-4)
#0.06666667
#6.7% of the transcripts (9) were not measurable in RNASeq, but were measureable using qPCR.

#I would add that there seems to be another 4 that are above 30 Cq and don't look like their measurements in DLPFC are reliable:
#In general, it seems like the relationship between Cq and log2CPM falls apart above 30. What are those transcripts?

AverageExpr_FrontalPoleVsRNASeq$Gene.Symbol[AverageExpr_FrontalPoleVsRNASeq$AverageCq_PerGene>30]
# [1] "GABRR1"  "GABRR2"  "GRM6"    "DDC"     "DRD3"    "DRD4"    "HTR2B"   "HTR6"    "HTR7"    "SLC18A1" "SLC18A2" "SLC6A3" 
# [13] "SLC6A4"  "TH" 
#14 genes total.
#Huh. But 5 of these genes were measurable in the RNA-Seq dataset. 14-5=9, 13-9=4.
#That means that at least 4 of the NAs in the RNAseq data came from genes with higher expression. 
#18S (duh), ADORA2A (Cq 27), DRD2 (Cq 29), HTR3A (Cq 29.7)
#interesting.

AverageExpr_FrontalPoleVsRNASeq$Gene.Symbol[AverageExpr_FrontalPoleVsRNASeq$AverageCq_PerGene>30 & is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr)==FALSE]
#[1] "GRM6"   "HTR2B"  "HTR6"   "HTR7"   "SLC6A4"
#These genes were detectable in previous RNA-Seq studies despite low-level expression.


#That's pretty fun too:

pdf("qPCR_AverageCq_Vs_MeasurableInRNASeq.pdf", width=5, height=5)
boxplot(AverageExpr_FrontalPoleVsRNASeq$AverageCq_PerGene~is.na(AverageExpr_FrontalPoleVsRNASeq$CMC_AveExpr)==FALSE, ylab="qPCR: Average Cq per Gene", xlab="Included in RNASeq Results", ylim=c(17, 37))
dev.off()
#Note that 18s is still included in this plot - I just set the ylim to exclude it, since it is inapplicable.

#So qPCR is able to reliably capture low-expression transcripts in a manner that is an improvement on RNA-Seq. Not exactly rocket science, but still fun to see in action.
#Also, some of the gene families that are more measurable with PCR than RNA-Seq are the ones that we are most interested in relationship to Schiz/BP in the frontal cortex: DA system, 5HT system.  Actually, seems like receptors in general seem to tend to have lower expression.

table(AverageCq_NumberNA_PerGene_DF$NumberNA_PerGene, trunc(AverageCq_NumberNA_PerGene_DF$AverageCq_PerGene))

pdf("qPCR_NumberNA_vsCq_perTranscript.pdf", width=5, height=5)
plot(AverageCq_NumberNA_PerGene_DF$NumberNA_PerGene~AverageCq_NumberNA_PerGene_DF$AverageCq_PerGene, ylab="% NA for Transcript", xlab="Average Cq for Transcript", pch=(as.numeric(AverageCq_NumberNA_PerGene_DF$NumberNA_PerGene>10)+0)*16, xlim=c(7,38))
dev.off()

####

##It might be interesting to run this comparison with the frontal pole microarray data too.

#I guess we want the full frontal pole microarray data for this correlation (i.e., not just the target genes)

AverageExpr_FrontalPoleMicroarray<-data.frame(Gene.Symbol=FrontalPole_Microarray_AllResults$GeneSymbol, Iwamoto_AveExpr=FrontalPole_Microarray_AllResults$Iwamoto_Model3_A, Maycox_AveExpr=FrontalPole_Microarray_AllResults$Maycox_Model3_A, stringsAsFactors=FALSE)
str(AverageExpr_FrontalPoleMicroarray)
# 'data.frame':	20575 obs. of  3 variables:
#   $ Gene.Symbol    : chr  "NAT2" "ADA" "CDH2" "AKT3" ...
# $ Iwamoto_AveExpr: num  2.91 6.09 6.41 6.63 3.45 ...
# $ Maycox_AveExpr : num  7.09 7.17 9.59 8.11 7.41 ...

library(plyr)

str(AverageCq_NumberNA_PerGene_DF)
# 'data.frame':	144 obs. of  3 variables:
#   $ Gene.Symbol      : chr  "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# $ AverageCq_PerGene: num  21.1 25.5 24.1 27.2 23.2 ...
# $ NumberNA_PerGene : int  0 0 0 0 0 0 0 0 0 0 ...

AverageExprPCR_vs_FrontalPoleMicroarray<-join(AverageCq_NumberNA_PerGene_DF, AverageExpr_FrontalPoleMicroarray, by="Gene.Symbol", type="left", match="all")
str(AverageExprPCR_vs_FrontalPoleMicroarray)

# 'data.frame':	144 obs. of  5 variables:
#   $ Gene.Symbol      : chr  "ABAT" "ADCY7" "ADORA1" "ADORA2A" ...
# $ AverageCq_PerGene: num  21.1 25.5 24.1 27.2 23.2 ...
# $ NumberNA_PerGene : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Iwamoto_AveExpr  : num  6.44 3.87 8.3 NA 7.88 ...
# $ Maycox_AveExpr   : num  10.86 6.41 8.83 NA 9.09 ...

setwd("~/Documents/Microarray Gen/FrontalPole/Output/FinalOutput_AcrossDatasets")
write.csv(AverageExprPCR_vs_FrontalPoleMicroarray, "AverageExprPCR_vs_FrontalPoleMicroarray.csv")

#Oooh - that's pretty cool:

pdf("Iwamoto_AVEExpress_vs_qPCR_AverageCq.pdf", width=5, height=5)
plot(AverageExprPCR_vs_FrontalPoleMicroarray$Iwamoto_AveExpr~AverageExprPCR_vs_FrontalPoleMicroarray$AverageCq_PerGene, xlab="qPCR: Average Cq per Gene", ylab="Iwamoto Microarray: Average Log2 Signal", xlim=c(16, 37))
dev.off()

pdf("Maycox_AVEExpress_vs_qPCR_AverageCq.pdf", width=5, height=5)
plot(AverageExprPCR_vs_FrontalPoleMicroarray$Maycox_AveExpr~AverageExprPCR_vs_FrontalPoleMicroarray$AverageCq_PerGene, xlab="qPCR: Average Cq per Gene", ylab="Maycox Microarray: Average Log2 Signal", xlim=c(16, 37))
dev.off()

pdf("Maycox_AVEExpress_vs_Iwamoto_AVEExpr.pdf", width=5, height=5)
plot(AverageExprPCR_vs_FrontalPoleMicroarray$Maycox_AveExpr~AverageExprPCR_vs_FrontalPoleMicroarray$Iwamoto_AveExpr, xlab="Iwamoto Microarray: Average log2 Signal", ylab="Maycox Microarray: Average Log2 Signal")
dev.off()

#compare to RNA-Seq?

############################

#Determining which of the genes in our dataset have previous indications that there may be genetic variation associated with Schiz

#Genetic info from SZDB_2
#http://www.szdb.org/integrative2.php

setwd("~/Documents/Microarray Gen/FrontalPole/SZDB_2")

list.files()
# [1] "~$FrontalPole_SZDB2.xlsx"  "CNV_genes.txt"             "Diff-methy"                "Diff-methy.zip"           
# [5] "ExomeSequencing_Genes.txt" "FrontalPole_SZDB2.xlsx"    "GWAS-Genes"                "GWAS-Genes.zip"           
# [9] "integrative"               "integrative.zip"  

#Data from "Contribution of copy number variants to schizophrenia from a genome-wide study of 41,321 subjects" pmid:27869829
#I had to fix some formatting issues in the column names before reading this in
SCHIZ_CNV_genes<-read.csv("CNV_genes_forR.csv", header=TRUE, stringsAsFactors = FALSE)

str(SCHIZ_CNV_genes)

# 'data.frame':	541 obs. of  16 variables:
#   $ gene              : chr  "CA15P1" "DGCR2" "Y_RNA" "DGCR11" ...
# $ genechr           : chr  "chr22" "chr22" "chr22" "chr22" ...
# $ genestart         : int  19019077 19023795 19032769 19033675 19043314 19109042 19111822 19112392 19117792 19118321 ...
# $ geneend           : int  19022077 19109967 19032880 19035888 19044025 19112016 19115962 19113159 19132197 19120134 ...
# $ gene_type         : chr  "pseudogene" "protein_coding" "misc_RNA" "sense_intronic" ...
# $ cytoband          : chr  "22q11.21" "22q11.21" "22q11.21" "22q11.21" ...
# $ cnvchr            : chr  "chr22" "chr22" "chr22" "chr22" ...
# $ cnvstart          : int  19020000 19020000 19020000 19020000 19020000 19020000 19020000 19020000 19020000 19020000 ...
# $ cnvend            : int  21420000 21420000 21420000 21420000 21420000 21420000 21420000 21420000 21420000 21420000 ...
# $ Putative_mechanism: chr  "NAHR" "NAHR" "NAHR" "NAHR" ...
# $ CNV_test          : chr  "Loss" "Loss" "Loss" "Loss" ...
# $ Direction         : chr  "Risk" "Risk" "Risk" "Risk" ...
# $ case_number       : chr  "64/21094" "64/21094" "64/21094" "64/21094" ...
# $ control_number    : chr  "1/20227" "1/20227" "1/20227" "1/20227" ...
# $ Pregional         : num  5.7e-18 5.7e-18 5.7e-18 5.7e-18 5.7e-18 ...
# $ OR                : chr  "67.7 (9.3-492.8)" "67.7 (9.3-492.8)" "67.7 (9.3-492.8)" "67.7 (9.3-492.8)" ...

#Total # of genes implicated:

length(unique(SCHIZ_CNV_genes$gene))
#[1] 408


dim(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC)
#[1] 111 293
#Target genes that survived QC for our qPCR experiment.
sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_CNV_genes$gene)
#[1] 1
#COMT
#As an aside: the reference gene TFRC was also implicated by a SCHIZ CNV study

SCHIZ_ExomeSequencing_Genes<-read.delim("ExomeSequencing_Genes.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

str(SCHIZ_ExomeSequencing_Genes)

# 'data.frame':	9365 obs. of  18 variables:
#   $ Author              : chr  "Xu" "Xu" "Xu" "Xu" ...
# $ Year                : int  2011 2011 2011 2011 2011 2011 2011 2011 2011 2011 ...
# $ PMID                : int  21822266 21822266 21822266 21822266 21822266 21822266 21822266 21822266 21822266 21822266 ...
# $ Exome.Capture.Kit   : chr  "Agilent SureSelect Human All Exon Target Enrichment System" "Agilent SureSelect Human All Exon Target Enrichment System" "Agilent SureSelect Human All Exon Target Enrichment System" "Agilent SureSelect Human All Exon Target Enrichment System" ...
# $ Sequence.Platform   : chr  "Illumina paired-end" "Illumina paired-end" "Illumina paired-end" "Illumina paired-end" ...
# $ Case.Number         : chr  "53 trios" "53 trios" "53 trios" "53 trios" ...
# $ Control.Number      : chr  "22 trios" "22 trios" "22 trios" "22 trios" ...
# $ Sample.Description  : chr  "proband diagnosed with schizophrenia with unaffected parents" "proband diagnosed with schizophrenia with unaffected parents" "proband diagnosed with schizophrenia with unaffected parents" "proband diagnosed with schizophrenia with unaffected parents" ...
# $ Familial.or.Sporadic: chr  "Sporadic family" "Sporadic family" "Sporadic family" "Sporadic family" ...
# $ Ancestry            : chr  "European" "European" "European" "European" ...
# $ Amino.Acid.Change   : chr  "p.C30S" "p.R1081H" "p.G539R" "p.G202R" ...
# $ Mutation.Type       : chr  "SNV" "SNV" "SNV" "SNV" ...
# $ Annotation          : chr  "missense" "missense" "missense" "missense" ...
# $ De.novo.            : chr  "de novo" "de novo" "de novo" "de novo" ...
# $ PoluPhen.2          : chr  "Probably damaging" "Probably damaging" "Probably damaging" "Probably damaging" ...
# $ Gene                : chr  "PLCL2" "WDR11" "DPYD" "OR4C46" ...
# $ Position            : chr  "chr3:17051253; absent in 1000 Genomes Project" "chr10:122664879; absent in 1000 Genomes Project" "chr1:97981407; absent in 1000 Genomes Project" "chr11:51515885; absent in 1000 Genomes Project" ...
# $ SIFT                : chr  "Tolerated" "Damaging" "Damaging" "Damaging" ...

# of unique genes implicated:
length(unique(SCHIZ_ExomeSequencing_Genes$Gene))
#[1] 5830

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_ExomeSequencing_Genes$Gene)
#[1] 45
45/111
#[1] 0.4054054
#!!! 40% of the genes in our target list (that survived QC)

setwd("~/Documents/Microarray Gen/FrontalPole/SZDB_2/GWAS-Genes")

list.files()
#[1] "CLOZUK_GENE.txt" "PGC2_GENE.txt" 

SCHIZ_CLOZUK_GENE<-read.delim("CLOZUK_GENE.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
str(SCHIZ_CLOZUK_GENE)
# 'data.frame':	480 obs. of  11 variables:
#   $ Gene       : chr  "MHC" "CACNA1C" "CACNA1C-AS4" "CACNA1C-IT3" ...
# $ Biotype    : chr  "-" "protein_coding" "antisense" "sense_intronic" ...
# $ Gene_start : chr  "28477797" "2079952" "2329703" "2378942" ...
# $ Gene_end   : chr  "33448354" "2802108" "2332645" "2398103" ...
# $ Index_SNP  : chr  "rs3130820" "rs2007044" "rs2007044" "rs2007044" ...
# $ P.value    : num  2.12e-44 5.63e-20 5.63e-20 5.63e-20 1.10e-18 ...
# $ CLOZUK_loci: int  1 2 2 2 3 3 4 4 4 4 ...
# $ chr        : chr  "6" "12" "12" "12" ...
# $ Loci_start : int  24988105 2321868 2321868 2321868 1877502 1877502 98341152 98341152 98341152 98341152 ...
# $ Loci_end   : int  33842877 2523772 2523772 2523772 2190100 2190100 98559093 98559093 98559093 98559093 ...
# $ GWS.in.PGC2: chr  "Yes" "Yes" "Yes" "Yes" ...

length(unique(SCHIZ_CLOZUK_GENE$Gene))
#[1] 475

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_CLOZUK_GENE$Gene)
#[1] 6
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_CLOZUK_GENE$Gene]
#[1] "GABBR2"  "GRIN2A"  "GRM3"    "SLC32A1" "SRR"     "DRD2"  

SCHIZ_PGC2_GENE<-read.delim("PGC2_GENE.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
str(SCHIZ_PGC2_GENE)
# 'data.frame':	350 obs. of  11 variables:
#   $ Gene         : chr  "MHC" "DPYD" "MIR137" "ARL3" ...
# $ Biotype      : chr  "-" "protein_coding" "miRNA" "protein_coding" ...
# $ Gene_start   : chr  "-" "97543299" "98511626" "104433488" ...
# $ Gene_end     : chr  "-" "98386605" "98511727" "104474164" ...
# $ Index_SNP    : chr  "rs1233578" "rs2660304" "rs2660304" "rs12416331" ...
# $ P.value      : num  3.48e-31 3.36e-19 3.36e-19 6.20e-19 6.20e-19 ...
# $ CLOZUK_loci  : int  1 2 2 3 3 3 3 3 3 3 ...
# $ chr          : chr  "6" "1" "1" "10" ...
# $ Loci_start   : int  28303247 97792625 97792625 104423800 104423800 104423800 104423800 104423800 104423800 104423800 ...
# $ Loci_end     : int  28712247 98559084 98559084 105165583 105165583 105165583 105165583 105165583 105165583 105165583 ...
# $ GWS.in.CLOZUK: chr  "Yes" "Yes" "Yes" "Yes" ...

length(unique(SCHIZ_PGC2_GENE$Gene))
#[1] 348

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_PGC2_GENE$Gene)
#[1] 6
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_PGC2_GENE$Gene]
#[1] "GRIA1"   "GRIN2A"  "GRM3"    "SLC32A1" "SRR"     "DRD2" 
#5 are the same as what was found in CLOZUK

#Total of our targets sig in either CLOZUK or PGC2
7/111
#[1] 0.06306306

length(unique(c(SCHIZ_CLOZUK_GENE$Gene, SCHIZ_PGC2_GENE$Gene)))
#[1] 571

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_PGC2_GENE<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_PGC2_GENE$Gene

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_CLOZUK_GENE<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_CLOZUK_GENE$Gene

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_ExomeSequencing_Genes<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_ExomeSequencing_Genes$Gene

#Maybe it would be good to note how many times a gene pops up for each of these types of evidence?
table(SCHIZ_ExomeSequencing_Genes$Gene)
#Some genes definitely pop up more than others... but were these studies targeted? (vs. unbiased?)
length(unique(SCHIZ_ExomeSequencing_Genes$Gene))
#[1] 5830

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_ExomeSequencing_Genes)
#[1] 45

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_CNV_Genes<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_CNV_genes$gene

table(SCHIZ_CNV_genes$gene)
length(unique(SCHIZ_CNV_genes$gene))
#[1] 408

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_CNV_Genes)
#[1] 1

#Getting into the predictive/integrative analyses (GWAS+eQTL) - there are a lot of files for these:
#I only downloaded the SMR output (not the sherlock) for each combination of GWAS + Brain/Blood eQTL database.

setwd("~/Documents/Microarray Gen/FrontalPole/SZDB_2/integrative")
list.files()
# [1] "SMR_CLOZUK_Brain_Amygdala.txt"                        "SMR_CLOZUK_Brain_Anterior_cingulate_cortex_BA24.txt" 
# [3] "SMR_CLOZUK_Brain_Caudate_basal_ganglia.txt"           "SMR_CLOZUK_Brain_Cerebellar_Hemisphere.txt"          
# [5] "SMR_CLOZUK_Brain_Cerebellum.txt"                      "SMR_CLOZUK_Brain_Cortex.txt"                         
# [7] "SMR_CLOZUK_Brain_eMeta.txt"                           "SMR_CLOZUK_Brain_Frontal_Cortex_BA9.txt"             
# [9] "SMR_CLOZUK_Brain_Hippocampus.txt"                     "SMR_CLOZUK_Brain_Hypothalamus.txt"                   
# [11] "SMR_CLOZUK_Brain_mMeta.txt"                           "SMR_CLOZUK_Brain_Nucleus_accumbens_basal_ganglia.txt"
# [13] "SMR_CLOZUK_Brain_Putamen_basal_ganglia.txt"           "SMR_CLOZUK_Brain_Spinal_cord_cervical.txt"           
# [15] "SMR_CLOZUK_Brain_Substantia_nigra.txt"                "SMR_CLOZUK_cage.txt"                                 
# [17] "SMR_CLOZUK_GTEx_brain.txt"                            "SMR_CLOZUK_Hannon_blood1.txt"                        
# [19] "SMR_CLOZUK_Hannon_blood2.txt"                         "SMR_CLOZUK_Hannon_FetalBrain.txt"                    
# [21] "SMR_CLOZUK_LBC_BSGS.txt"                              "SMR_CLOZUK_Pituitary.txt"                            
# [23] "SMR_CLOZUK_psychENCODE.txt"                           "SMR_CLOZUK_westra.txt"                               
# [25] "SMR_PGC2_Brain_Amygdala.txt"                          "SMR_PGC2_Brain_Anterior_cingulate_cortex_BA24.txt"   
# [27] "SMR_PGC2_Brain_Caudate_basal_ganglia.txt"             "SMR_PGC2_Brain_Cerebellar_Hemisphere.txt"            
# [29] "SMR_PGC2_Brain_Cerebellum.txt"                        "SMR_PGC2_Brain_Cortex.txt"                           
# [31] "SMR_PGC2_Brain_eMeta.txt"                             "SMR_PGC2_Brain_Frontal_Cortex_BA9.txt"               
# [33] "SMR_PGC2_Brain_Hippocampus.txt"                       "SMR_PGC2_Brain_Hypothalamus.txt"                     
# [35] "SMR_PGC2_Brain_mMeta.txt"                             "SMR_PGC2_Brain_Nucleus_accumbens_basal_ganglia.txt"  
# [37] "SMR_PGC2_Brain_Putamen_basal_ganglia.txt"             "SMR_PGC2_Brain_Spinal_cord_cervical.txt"             
# [39] "SMR_PGC2_Brain_Substantia_nigra.txt"                  "SMR_PGC2_cage.txt"                                   
# [41] "SMR_PGC2_GTEx_brain.txt"                              "SMR_PGC2_Hannon_blood1.txt"                          
# [43] "SMR_PGC2_Hannon_blood2.txt"                           "SMR_PGC2_Hannon_FetalBrain.txt"                      
# [45] "SMR_PGC2_LBC_BSGS.txt"                                "SMR_PGC2_Pituitary.txt"                              
# [47] "SMR_PGC2_psychENCODE.txt"                             "SMR_PGC2_westra.txt"     


SCHIZ_GWAS_eQTL_Integrative_Files<-list.files()

SCHIZ_GWAS_eQTL_Integrative_Results<-read.csv("SMR_CLOZUK_Brain_Amygdala.txt", header=FALSE, sep="\t", stringsAsFactors = FALSE) 

str(SCHIZ_GWAS_eQTL_Integrative_Results)

#Doublechecking if they all have the same structure.
for(i in c(1: length(SCHIZ_GWAS_eQTL_Integrative_Files))){
  temp<-read.csv(SCHIZ_GWAS_eQTL_Integrative_Files[i], header=FALSE, sep="\t", stringsAsFactors = FALSE) 
  
  print(dim(temp))
  
}
#... and they don't. :( 
#some have 10 columns, some have 11. Let's compare and see if there is a pattern:

colnames(SCHIZ_GWAS_eQTL_Integrative_Results)
#[1] "V1"  "V2"  "V3"  "V4"  "V5"  "V6"  "V7"  "V8"  "V9"  "V10" "V11"
colnames(read.csv(SCHIZ_GWAS_eQTL_Integrative_Files[2], header=FALSE, sep="\t", stringsAsFactors = FALSE) )
#[1] "V1"  "V2"  "V3"  "V4"  "V5"  "V6"  "V7"  "V8"  "V9"  "V10"
#Well that's informative, LOL. Duh...
SCHIZ_GWAS_eQTL_Integrative_Results[1,]
# V1                V2           V3         V4           V5       V6       V7       V8       V9 V10 V11
# 1 RP1-140A9.1 ENSG00000231050.1 chr1:1823503 rs28423370 chr1:1855450 1.10-e03 7.64-e13 3.06-e03 3.98-e01  10  NA
(read.csv(SCHIZ_GWAS_eQTL_Integrative_Files[2], header=FALSE, sep="\t", stringsAsFactors = FALSE) )[1,]
# V1                V2            V3        V4            V5       V6       V7       V8       V9 V10
# 1 CLCNKA ENSG00000186510.7 chr1:16352957 rs3738633 chr1:16359827 3.19-e02 4.56-e10 4.27-e02 8.17-e02   9

#Looks like the final column in the files with 11 columns might be artifact.

SCHIZ_GWAS_eQTL_Integrative_Results[,11]
# [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
# [39] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
# [77] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
# [115] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#Yep.

SCHIZ_GWAS_eQTL_Integrative_Results<-SCHIZ_GWAS_eQTL_Integrative_Results[,c(1:10)]

for(i in c(3: length(SCHIZ_GWAS_eQTL_Integrative_Files))){
  temp<-read.csv(SCHIZ_GWAS_eQTL_Integrative_Files[i], header=FALSE, sep="\t", stringsAsFactors = FALSE) 
  
  SCHIZ_GWAS_eQTL_Integrative_Results<-rbind.data.frame(SCHIZ_GWAS_eQTL_Integrative_Results, temp[,c(1:10)])  
  
}

dim(SCHIZ_GWAS_eQTL_Integrative_Results)
#[1] 98067    10

setwd("~/Documents/Microarray Gen/FrontalPole/SZDB_2")

write.csv(SCHIZ_GWAS_eQTL_Integrative_Results, "SCHIZ_GWAS_eQTL_Integrative_Results.csv")

SCHIZ_GWAS_eQTL_Integrative_Results[1,]

colnames(SCHIZ_GWAS_eQTL_Integrative_Results)[1]<-"gene"
  
length(unique(SCHIZ_GWAS_eQTL_Integrative_Results$gene))
#[1] 8578
#hmmm... that's definitely likely to contain some noise. 

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_GWAS_eQTL_Integrative<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%SCHIZ_GWAS_eQTL_Integrative_Results$gene

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$IN_SCHIZ_GWAS_eQTL_Integrative)
#[1] 65
65/111
#[1] 0.5855856
#So approximately half. Maybe not much better than random chance?


#I would like to add in the linkage and methylation data, but I'm having trouble downloading it in bulk from the SZDB.org website.
#I downloaded the results for our genes by hand - that will provide annotation, but not context with how to interpret it.

LinkageAssociation_JustOurGenes<-read.csv("LinkageAssociation_JustOurGenes.csv", header=TRUE, stringsAsFactors = FALSE)
str(LinkageAssociation_JustOurGenes)
# 'data.frame':	83 obs. of  8 variables:
#   $ Gene      : chr  "BDNF" "BDNF" "BDNF" "COMT" ...
# $ Biotype   : chr  "protein-coding" "protein-coding" "protein-coding" "protein-coding" ...
# $ Locationa : chr  "chr11:27676440..27743605" "chr11:27676440..27743605" "chr11:27676440..27743605" "chr22:19929263..19957498" ...
# $ Full.Name : chr  "brain derived neurotrophic factor" "brain derived neurotrophic factor" "brain derived neurotrophic factor" "catechol-O-methyltransferase" ...
# $ Study.Type: chr  "linkage" "linkage" "Association" "Association" ...
# $ Authors   : chr  "Lewis et al." "Ng et al." "Sun et al." "Allen et al." ...
# $ Year      : int  2003 2009 2008 2008 2003 2009 2008 2009 2008 2003 ...
# $ PubMed    : int  12802786 19349958 18361404 18583979 12802786 19349958 18361404 19349958 18583979 12802786 ...

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_SCHIZ_LinkageAssociation<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%LinkageAssociation_JustOurGenes$Gene

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_SCHIZ_LinkageAssociation)
#[1] 30

DiffMethylation_JustOurGenes<-read.csv("DiffMethylation_JustOurGenes.csv", header=TRUE, stringsAsFactors = FALSE)
str(DiffMethylation_JustOurGenes)

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_SCHIZ_DiffMethylation<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%DiffMethylation_JustOurGenes$Gene

sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_SCHIZ_DiffMethylation)
#[1] 6

head(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[,c(294:300)])

table(apply(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[,c(294:300)], 1, sum))
# 0  1  2  3  4  5 
# 25 35 34 14  1  2 

#Very few of our target genes lack annotation linking them genetically to SCHIZ.

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[apply(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[,c(294:300)], 1, sum)<1]
# [1] "ABAT"    "ADORA2A" "SST"     "GABRE"   "GABRG1"  "GABRG3"  "GABRQ"   "GNAQ"    "GRIA2"   "GRIN2C"  "GRM2"   
# [12] "HOMER1"  "AQP4"    "GJA1"    "SLC17A7" "SLC17A8" "SLC1A6"  "SLC6A11" "ADRB2"   "DBH"     "HTR1B"   "HTR1D"  
# [23] "HTR2C"   "HTR3B"   "HTR6"
#... and interestingly the genes lacking any annotation suggesting a relationship with diagnosis are some of the top genes altered in our dataset (and other datasets) - ABAT, SST, AQP4. Huh.  Maybe those genes are more related to stress/substance use/medication?

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[apply(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[,c(294:300)], 1, sum)>2]
# [1] "BDNF"    "CACNA1A" "GABBR1"  "GABRA1"  "GRIA1"   "GRIN1"   "GRIN2A"  "GRIN2B"  "GRM3"    "SLC32A1" "SRR"    
# [12] "COMT"    "DRD1"    "DRD2"    "DRD3"    "HTR2A"   "SLC18A2"
#These genes have the most different types of genetic evidence backing up their association with schiz.

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[apply(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC[,c(294:300)], 1, sum)>3]
#[1] "GRIN2A" "GRM3"   "DRD2" 
#... and especially these.

setwd("~/Documents/Microarray Gen/FrontalPole/SZDB_2")

write.csv(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, "FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_wGenetics.csv")

################################

#Hmm... let's see if we can build a case for some of these genes showing expression differences due to stress instead.

#Genes identified in association with stress in Flati 2020 meta-analysis of mice RNA-Seq data from across the brain.
setwd("~/Documents/Literature/Flati_2020_StressRNASeq_MetaAnalysis")
list.files()
# [1] "Flati_2020_41597_2020_772_MOESM1_ESM.xlsx"                  
# [2] "Flati_2020_41597_2020_772_MOESM2_ESM.xlsx"                  
# [3] "Flati_2020_41597_2020_772_MOESM3_ESM.xlsx"                  
# [4] "Flati_2020_41597_2020_772_MOESM4_ESM.xlsx"                  
# [5] "Flati_2020_41597_2020_772_MOESM5_ESM.xlsx"                  
# [6] "Flati_2020_41597_2020_772_MOESM6_ESM.xlsx"                  
# [7] "Flati_2020_41597_2020_772_MOESM7_ESM.xlsx"                  
# [8] "Flati_2020_41597_2020_772_MOESM8_ESM.pdf"                   
# [9] "Flati_2020_StressRNASeq_MetaAnalysis_s41597-020-00772-z.pdf"
# [10] "Stress_Mice_CS.csv"                                         
# [11] "Stress_Mice_HIPP.csv"                                       
# [12] "Stress_Mice_SUPP1.csv"                                      
# [13] "Stress_Mice_SUPP2.csv"                                      
# [14] "Stress_Mice_T3.csv" 

Stress_Mice_T3<-read.csv("Stress_Mice_T3.csv", header=FALSE, stringsAsFactors = FALSE)
str(Stress_Mice_T3)
# 'data.frame':	109 obs. of  1 variable:
#   $ V1: chr  "Islr2" "Btg2" "Col11a2" "Dusp1" ...

Stress_Mice_SUPP1<-read.csv("Stress_Mice_SUPP1.csv", header=FALSE, stringsAsFactors = FALSE)
str(Stress_Mice_SUPP1)
# 'data.frame':	364 obs. of  1 variable:
#   $ V1: chr  "1110002L01Rik" "1110032F04Rik" "1700020I14Rik" "3110035E14Rik" ...

Stress_Mice_SUPP2<-read.csv("Stress_Mice_SUPP2.csv", header=FALSE, stringsAsFactors = FALSE)
str(Stress_Mice_SUPP2)
# 'data.frame':	1452 obs. of  1 variable:
#   $ V1: chr  "Supplementary Table 2: DEGs shared among 2 BioProjects" "1110008P14Rik" "1500015O10Rik" "1700001L05Rik" ...

#Digging up orthologs
setwd("~/Documents/Microarray Gen/GeneOrthology_HumanRatMouse")
Orthologs_Mice_Humans<-read.delim("Orthologs_Mice_Humans_Symbols_NoNA_InformaticsJax_OrthologDatabase_20210228.txt", header=TRUE, sep=",", stringsAsFactors = FALSE)
str(Orthologs_Mice_Humans)

Stress_Mice<-data.frame(Symbol_Mouse=c(Stress_Mice_T3[,1],Stress_Mice_SUPP1[,1], Stress_Mice_SUPP2[,1]), StressEvidence=c(rep("MoreThan4", length(Stress_Mice_T3[,1])), rep("3studies", length(Stress_Mice_SUPP1[,1])),  rep("2studies", length(Stress_Mice_SUPP2[,1]))))

library(plyr)

Stress_Mice<-join(Stress_Mice, Orthologs_Mice_Humans, by="Symbol_Mouse", type="left", match="all")
dim(Stress_Mice)
#[1] 1931    4
length(unique(Stress_Mice$Symbol_Mouse))
#[1] 1925

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_Stress_MetaAnalysis<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol%in%(c(toupper(Stress_Mice$Symbol_Mouse),Stress_Mice$Symbol_Human))
  
sum(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_Stress_MetaAnalysis)
#[1] 27
27/111
#[1] 0.2432432
#So probably not just noise.
FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$In_Stress_MetaAnalysis]
# [1] "ADCY7"   "ADORA1"  "ADORA2A" "SST"     "PVALB"   "BDNF"    "GABRA4"  "GABRD"   "GABRG1"  "GABRG3"  "GABRQ"  
# [12] "GRIN2A"  "GRIN2B"  "GRM2"    "HOMER2"  "PHGDH"   "SHANK2"  "SLC17A6" "SLC17A7" "SLC6A11" "SLC6A13" "SLC7A11"
# [23] "DRD1"    "DRD2"    "HTR1A"   "HTR2A"   "HTR2C" 

#Interesting - SST and PVALB. That might explain why they have such strong effects across datasets but no strong genetic linkage to the disorder.
#Looking at the tables, both have shown up as stress-related in 3 other studies.

#DRD2 and BDNF look like they might be a point of convergence between genetic and stress.

4/27
#[1] 0.1481481
#15% of the stress associated target genes are arguably differentially expressed in the frontal pole (either FDR<0.1 or p<0.05 in at least 2 frontal cortex datasets)
12/111
#[1] 0.1081081
#11% of the target genes are DE (by above definition)
#So not a clear enrichment.


library(VennDiagram)

pdf("VennDiagram_OurTargets_vs_StressMicePortal.pdf", height=4.5, width=4.5)
grid.newpage()
draw.pairwise.venn(111, 1925, 27, category = c("Our qPCR Targets", "Stress Associated"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cat.fontfamily = rep("sans", 2))
dev.off()


setwd("~/Documents/Microarray Gen/FrontalPole/CompareWDLPFC")

write.csv(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, "FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_wGeneticsStress.csv")

#####################

#How about drug exposure?

#Referencing Martin 2015 for effects of typical and atypical antipsychotics and PCP:
setwd("~/Documents/Literature/Martin_2015_AntiPsychotics")

Martin_Suppl_Table1<-read.csv("Martin_Suppl_Table1.csv", header=TRUE, stringsAsFactors = FALSE)
Martin_Suppl_Table2<-read.csv("Martin_Suppl_Table2.csv", header=TRUE, stringsAsFactors = FALSE)

str(Martin_Suppl_Table1)

colnames(Martin_Suppl_Table1)<-paste(colnames(Martin_Suppl_Table1), "S1", sep="_")
colnames(Martin_Suppl_Table2)<-paste(colnames(Martin_Suppl_Table2), "S2", sep="_")

colnames(Martin_Suppl_Table1)[2]<-"Gene.Symbol"
colnames(Martin_Suppl_Table2)[2]<-"Gene.Symbol"

#Let's trim these down a little bit to make this less overwhelming. 
#Acute treatment isn't very interesting for our applications.

Martin_Suppl_Table1_forJoin<-Martin_Suppl_Table1[,c(2,5,7,20, 22,23, 26, 28,29,32,34,35, 38, 40,41)]
Martin_Suppl_Table2_forJoin<-Martin_Suppl_Table2[,c(2,5,7,20, 22,23, 26, 28,29,32,34,35, 38, 40,41)]

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin<-join_all(list(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC, Martin_Suppl_Table1_forJoin, Martin_Suppl_Table2_forJoin), by="Gene.Symbol", type="left", match="all")

dim(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin)
#[1] 195 329
#Definitely some multimapped genes.
length(unique(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Gene.Symbol[is.na(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$p.value.Treatment._S1)==FALSE]))
#[1] 78
#There are 78 of our genes included in S1

length(unique(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Gene.Symbol[is.na(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$p.value.Treatment._S2 )==FALSE]))
#[1] 78
#There are 78 of our genes included in S2



#I wonder how many of the genes have both Table 1 and Table 2 representation:
sum((is.na(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$p.value.Treatment._S1)==FALSE) & (is.na(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$p.value.Treatment._S2 )==FALSE))
#[1] 162

length(unique(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Gene.Symbol[(is.na(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$p.value.Treatment._S1)==FALSE) & (is.na(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$p.value.Treatment._S2 )==FALSE)]))
#[1] 78
#Looks like the same ones

#How well do the values correlate between the two studies?

setwd("~/Documents/Microarray Gen/FrontalPole/Drugs")

pdf("Scatterplot_Martin2015_PCP_OurGenes_S1vsS2.pdf", height=5, width=5)
plot(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH._S1~FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH._S2, xlab="S2", ylab="S1", main="PCP Ratio")
dev.off()

cor(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH._S1,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH._S2, use="pairwise.complete")
#[1] 0.3389428

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH_AVE<-apply(cbind(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH._S1,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH._S2), 1, function(y) mean(y,na.rm=TRUE))

pdf("Scatterplot_Martin2015_HAL_OurGenes_S1vsS2.pdf", height=5, width=5)
plot(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH._S1~FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH._S2, xlab="S2", ylab="S1", main="HAL Ratio")
dev.off()

cor(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH._S1,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH._S2, use="pairwise.complete")
#[1] 0.4320719

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH_AVE<-apply(cbind(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH._S1,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH._S2), 1, function(y) mean(y,na.rm=TRUE))


pdf("Scatterplot_Martin2015_OLANZ_OurGenes_S1vsS2.pdf", height=5, width=5)
plot(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLZ.vs..VEH._S1~FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLANZ.vs..VEH._S2, xlab="S2", ylab="S1", main="OLANZ Ratio")
dev.off()

cor(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLZ.vs..VEH._S1,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLANZ.vs..VEH._S2, use="pairwise.complete")
#[1] 0.2548856

FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLANZ.vs..VEH_AVE<-apply(cbind(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLZ.vs..VEH._S1,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLANZ.vs..VEH._S2), 1, function(y) mean(y,na.rm=TRUE))

#How do they compare to each other?

pdf("Scatterplot_Martin2015_PCPvsHAL_OurGenes.pdf", height=5, width=5)
plot(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH_AVE~FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH_AVE, xlab="PCP", ylab="HAL", main="Gene Expression Ratio")
dev.off()
#Mostly a positive relationship, oddly enough.
cor(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.HAL.vs..VEH_AVE,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH_AVE, use="pairwise.complete")
#[1] 0.4850871

pdf("Scatterplot_Martin2015_PCPvsOLZ_OurGenes.pdf", height=5, width=5)
plot(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLANZ.vs..VEH_AVE~FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH_AVE, xlab="PCP", ylab="OLZ", main="Gene Expression Ratio")
dev.off()
#Also mostly a positive relationship, oddly enough.
cor(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.OLANZ.vs..VEH_AVE,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin$Ratio.PCP.vs..VEH_AVE, use="pairwise.complete")
#[1] 0.4124871

write.csv(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin, "FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC_vsMartin.csv")

#Hmm... some of those targets are represented quite a few times.  
#I think our best bet might be to come up with some indication as to whether the effect met even nominal significance in the Martin study
#And then average by target to get a sense of likely overall effect 

#####################

#Taking a peek at the functional pathways associated with the top genes:
#This is currently set up with rat annotation, but it looks like it doesn't matter for our genes (because they are so well-known/old, they have the same symbols across species)
#This also has a whole bunch of other databases added into it (from Gemma/GeneWeaver) for our NACC and HC experiment - ignore for now.

setwd("~/Documents/Microarray Gen/FrontalPole/fGSEA")

FunctionalOntology_ReferenceDatabase<-read.delim("c5withBrainCellTypesFunctionGemmaGeneWeaver_RatOrtholog.gmt.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)

str(FunctionalOntology_ReferenceDatabase)

FunctionalOntology_ReferenceDatabase<-unname(t(FunctionalOntology_ReferenceDatabase))

str(FunctionalOntology_ReferenceDatabase)

head(FunctionalOntology_ReferenceDatabase[,2])

FunctionalOntology_ReferenceDatabase[1,apply(FunctionalOntology_ReferenceDatabase, 2, function(y) "Htr2b"%in%y)]
#Nice. 
#There are a lot of them.  We should probably see how many of them show up repeatedly for our genes.

FunctionalAnnotation_ForqPCRGenes<-matrix(nrow=length(FunctionalOntology_ReferenceDatabase[1,]), ncol=length(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol))

colnames(FunctionalAnnotation_ForqPCRGenes)<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol
row.names(FunctionalAnnotation_ForqPCRGenes)<-FunctionalOntology_ReferenceDatabase[1,]

for(i in c(1:length(FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol))){
  for(j in c(1:length(FunctionalOntology_ReferenceDatabase[1,]))){
    FunctionalAnnotation_ForqPCRGenes[j,i]<-FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Gene.Symbol[i]%in%toupper(FunctionalOntology_ReferenceDatabase[,j])
  }
}

head(FunctionalAnnotation_ForqPCRGenes)
#Neat
ContainsGenesFromQPCR_FDR10<-apply(FunctionalAnnotation_ForqPCRGenes[,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Diagnosis_ANOVA_MLM_FDR<0.10], 1, sum)
table(ContainsGenesFromQPCR_FDR10)
#     0     1     2     3     4     5 
# 14748   985   120    49     7     5

ContainsGenesFromQPCR_P05<-apply(FunctionalAnnotation_ForqPCRGenes[,FrontalPole_qPCRGenes_MasterAnnotation_vs_FrontalPole_Microarray_vsDLPFC$Diagnosis_ANOVA_MLM_NominalPval<0.05], 1, sum)
table(ContainsGenesFromQPCR_P05)
# 0     1     2     3     4     5     6     7     8     9    10    11    13    14    15 
# 13601  1506   415   176    94    49    24    18    15     7     2     2     3     1     1 

write.csv(cbind.data.frame(row.names(FunctionalAnnotation_ForqPCRGenes), FunctionalAnnotation_ForqPCRGenes, ContainsGenesFromQPCR_FDR10,ContainsGenesFromQPCR_P05), "FunctionalAnnotation_ForqPCRGenes.csv") 
#Lots of uninteresting findings (i.e.default related to synapses...)
#Most interesting: 3 genes related to MAPK signalling, including DRD4 and HTR2B, 3 genes related to dopamine signalling.
#Also interesting is the overlap in the P05 set with genes previously identified as responsive to glucorticoids in particular cell types (Gemma)


################################

#I didn't re-do any of this yet:

##################################

#What about regional and or cell type specificity?

setwd("~/Documents/ABA_Amygdala")

ABA_Expression<-read.csv("ABA_AllExpressionsNoEmptyRow.csv", header=T, row.names=1, stringsAsFactors = F)
str(ABA_Expression)
ABA_Expression<-ABA_Expression[,-1]
ABA_Columns<-read.csv("ABA_Columns.csv", header=T, stringsAsFactors = F)
str(ABA_Columns)
ABA_Probes<-read.csv("ABA_AllProbesNoEmptyRow.csv", header=T, stringsAsFactors = F)

sum(ABA_Probes$gene.symbol%in%AverageExpr_FrontalPoleVsRNASeq$GeneSymbol)
#[1] 563
#Probably includes NAs
#Hmm... I should probably average by gene symbol. 

#But before doing that, let's remove the all NA columns because they're causing trouble.

ABA_Expression_NAbyCol<-apply(ABA_Expression, 2, function(y) sum(is.na(y)))
head(ABA_Expression_NAbyCol)
hist(ABA_Expression_NAbyCol)
#The columns either are all NA or have no NA
table(ABA_Expression_NAbyCol)
# 0 30000 
# 841   173 

dim(ABA_Expression)
#[1] 30000  1014
ABA_Expression_NoNA<-ABA_Expression[, which(ABA_Expression_NAbyCol<30000)]
dim(ABA_Expression_NoNA)
#[1] 30000   841
head(ABA_Expression_NoNA)

dim(ABA_Columns)
#[1] 1014   17
ABA_Columns_NoNA<-ABA_Columns[which(ABA_Expression_NAbyCol<30000),]
dim(ABA_Columns_NoNA)
#[1] 841  17

ABA_Expression_AveByGeneSymbol<-matrix(0, ncol=ncol(ABA_Expression_NoNA), nrow=length(names(table(ABA_Probes$gene.symbol))))
row.names(ABA_Expression_AveByGeneSymbol)<-names(table(ABA_Probes$gene.symbol))

for(i in c(1:ncol(ABA_Expression_NoNA))){
  ABA_Expression_AveByGeneSymbol[,i]<-tapply(ABA_Expression_NoNA[,i], ABA_Probes$gene.symbol, function(y) mean(y, na.rm=TRUE))
}

head(ABA_Expression_AveByGeneSymbol)
str(ABA_Expression_AveByGeneSymbol)
#side note: I probably could have used the collapseRows function for this instead of tapply...

write.csv(ABA_Expression_AveByGeneSymbol, "ABA_Expression_AveByGeneSymbol.csv")
write.csv(ABA_Columns_NoNA, "ABA_Columns_NoNA.csv")

#Cleaning up the workspace of large objects:
rm(ABA_Expression)
rm(ABA_Expression_NoNA)

str(ABA_Columns_NoNA)
table(ABA_Columns_NoNA$structure_name[ABA_Columns_NoNA$top_level_structure_abbreviation=="FL"])
# anterior orbital gyrus         frontal operculum                            frontal pole 
# 4                                       6                                       2 
# gyrus rectus  inferior frontal gyrus, opercular part    inferior frontal gyrus, orbital part 
# 6                                       4                                       6 
# inferior frontal gyrus, triangular part  inferior rostral gyrus                   lateral orbital gyrus 
# 6                                       5                                       6 
# medial orbital gyrus                    middle frontal gyrus       paracentral lobule, anterior part 
# 6                                       6                                       6 
# paraterminal gyrus                       parolfactory gyri                 posterior orbital gyrus 
# 1                                       5                                       5 
# precentral gyrus                  superior frontal gyrus                  superior rostral gyrus 
# 6                                       6                                       5 

#The Frontal pole only has two samples. :(

#What are these areas in BA?
#https://atlas.brain-map.org
#They call frontal pole area 10
#DLPFC is area 8, 9, 46
#None of this terminology fits their column labels...
#Ah - looks like the terminology is under "gyral"
#Immediately caudal to the frontal pole is the middle frontal gyrus (area 46, roughly) and superior frontal gyrus (area 9, roughly)

ABA_Columns_NoNA_JustFrontalPoleDLPFC<-ABA_Columns_NoNA[ABA_Columns_NoNA$structure_name%in%c("frontal pole", "middle frontal gyrus", "superior frontal gyrus"),]

ABA_Expression_AveByGeneSymbol_JustFrontalPoleDLPFC<-ABA_Expression_AveByGeneSymbol[,ABA_Columns_NoNA$structure_name%in%c("frontal pole", "middle frontal gyrus", "superior frontal gyrus")]

head(ABA_Expression_AveByGeneSymbol_JustFrontalPoleDLPFC)

colnames(ABA_Expression_AveByGeneSymbol_JustFrontalPoleDLPFC)<-paste(ABA_Columns_NoNA_JustFrontalPoleDLPFC$structure_name, ABA_Columns_NoNA_JustFrontalPoleDLPFC$donor_id)

#Since the Frontal Pole only has samples from two subjects, would it make sense to only compare within those two subjects?
ABA_Columns_NoNA[ABA_Columns_NoNA$structure_name=="frontal pole",]
#     donor_id donor_name donor_age donor_color sample_well sample_polygon sample_mri_0 sample_mri_1 sample_mri_2 structure_id
# 13     9861 H0351.2001  24 years      EC891D          97          28028           85           96           19         4888
# 14    10021 H0351.2002  39 years      0079C0        4561         980065          107           78           18         4888
# structure_name structure_abbreviation structure_color top_level_structure_id top_level_structure_name
# 13   frontal pole                     FP          E8C359                   4009             frontal lobe
# 14   frontal pole                     FP          E8C359                   4009             frontal lobe
# top_level_structure_abbreviation top_level_structure_color
# 13                               FL                    E8CD59
# 14                               FL                    E8CD59


ABA_Expression_PCRGenes<-ABA_Expression_AveByGeneSymbol_JustFrontalPoleDLPFC[which(row.names(ABA_Expression_AveByGeneSymbol_JustFrontalPoleDLPFC)%in%AverageExpr_FrontalPoleVsRNASeq$GeneSymbol),]
dim(ABA_Expression_PCRGenes)
#[1]  104  14
str(ABA_Expression_PCRGenes)
head(ABA_Expression_PCRGenes)

write.csv(ABA_Expression_PCRGenes, "ABA_Expression_PCRGenes.csv")

#Hard to look at, because all of the values differ on average by gene. Let's z-score it again:
#hmmm... that may be misleading - it emphasizes small differences, and de-emphasizes large differences. Maybe just center?
ABA_Expression_PCRGenes_Centered<-t(scale(t(ABA_Expression_PCRGenes), center=TRUE, scale=FALSE))
write.csv(ABA_Expression_PCRGenes_Centered, "ABA_Expression_PCRGenes_Centered.csv")

#Looks like there are big differences per region by person - maybe I should center by subject instead?
#This only works for the subjects with all three regions:
ABA_Expression_PCRGenes_Just2subjects<-ABA_Expression_PCRGenes[,ABA_Columns_NoNA_JustFrontalPoleDLPFC$donor_id%in%c(9861,10021)]
head(ABA_Expression_PCRGenes_Just2subjects)
ABA_Expression_PCRGenes_Just2subjects_Centered<-cbind(t(scale(t(ABA_Expression_PCRGenes_Just2subjects[,c(1,3,5)]), center=TRUE, scale=FALSE)), t(scale(t(ABA_Expression_PCRGenes_Just2subjects[,c(2,4,6)]), center=TRUE, scale=FALSE)))

write.csv(ABA_Expression_PCRGenes_Just2subjects_Centered, "ABA_Expression_PCRGenes_Just2subjects_Centered.csv")
#These results do not suggest that our targest are more expressed in the frontal pole than in the neighboring PFC regions - baseline expression for these molecules does not make the frontal pole special.

#oooh - look - there is a bioconductor package that does these sorts of analyses automatically!  I guess I didn't need to code it...
https://bioconductor.org/packages/release/bioc/vignettes/ABAEnrichment/inst/doc/ABAEnrichment.html

#Here's a meta-search engine for identifying spatiotemporal patterns and co-expression...but it doesn't include frontal pole:
http://www.brainexp.org/Search.aspx

#Many of our more interesting genes (e.g., HTR2B) have too low of expression to be present in many common datasets.


############################
#Additional notes:

#Depending on author/reviewer preferences, we may need to make a version of the boxplot figures that averages the replicate measurements for each subject
#Depending on author/reviewer preferences, we may also need to re-run the analysis relating gender and diagnosis on the original subject sample (vs. measurements sampled)


colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2)

colnames(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,1:82])

GabaGlu_NegDeltaDeltaCq_GeneGeneCorMatrix<-cor(GabaGlu_NegDeltaDeltaCq_AllSubjects_QCed2[,1:82], use="pairwise.complete.obs")

pdf("GabaGlu_NegDeltaDeltaCq_GeneGeneCorMatrix.pdf", height=12, width=12)
heatmap(GabaGlu_NegDeltaDeltaCq_GeneGeneCorMatrix)
dev.off()



colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2)
colnames(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,c(31:46)])

DA5HT_NegDeltaDeltaCq_GeneGeneCorMatrix<-cor(DA5HT_NegDeltaDeltaCq_AllSubjects_QCed2[,-c(31:46)], use="pairwise.complete.obs")
#Why is GAPDH included here? Ah - it is different from GADPH.

pdf("DA5HT_NegDeltaDeltaCq_GeneGeneCorMatrix.pdf", height=8, width=8)
heatmap(DA5HT_NegDeltaDeltaCq_GeneGeneCorMatrix)
dev.off()


