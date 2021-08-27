#This is some sample code to adapt for running the analysis determining how many of the drugs that our subjects were exposed to are things that are known to effect expression levels for our top differentially expressed genes
#Megan Hagenauer, 8-27-2021
#To be adapted by Evelyn Richardson


#Set working directory
#Can be done in the GUI (user interface - drop down menus) in Rstudio
#Or can be done by code using (replace YourDirectory with your actual directory...):
setwd(YourDirectory)
#If you use the GUI to set working directory, you may still want to copy and paste the code that appears in your console after you do that so that the wd is documented.

#This code tells you what files are in your working directory (prints them on your console):
list.files()

#Reading in Files and making them objects in the R environment
#Before doing this, save your files as single worksheets in .csv format
#It would be helpful to make sure that the worksheets have the same column name for the column that contains the list of drugs. 

#ExampleCode (header means column names)

ListOfDrugsForPritzker<-read.csv("YourFileName_ForDrugsInPritzker.csv", header=TRUE, stringsAsFactors=FALSE)

#Double-checking it read in properly:

#Tells you the structure - dimensions and datatypes - for your data.frame
str(ListOfDrugsForPritzker)

#Shows you a glimpse of the top of your data.frame
head(ListOfDrugsForPritzker)

ATopGene_DrugBudger_Database1<-read.csv("ATopGene_DrugBudger_Database1.csv", header=TRUE, stringsAsFactors=FALSE)


#Joining your list of drugs for Pritzker subjects with the drug budger output for each top gene for each drug budger database:
#This could be looped, but you can also do each by hand if you like (18 iterations? 6 top genes * 3 DrugBudger databases).

#You may need to install the plyr package
#If so, it may be easiest to do that with the GUI (drop down menus) in Rstudio

#load the plyr package:
library(plyr)

#Before joining, we will need a column in each of our data.frames that is equivalent - both the same column name and similar formatting for the drug names

str(ListOfDrugsForPritzker)
str(ATopGene_DrugBudger_Database1)

#To make formatting similar, we can force all of the drug names to be lower case using this code
#You will need to do this for all of your data.frames (19 iterations?), so could be looped or just done by hand (copy, paste, change data.frame name)
#To get this code to work, you will need to make the column name match what you see using the str() command

ATopGene_DrugBudger_Database1$'Drug Name'<-tolower(ATopGene_DrugBudger_Database1$'Drug Name') 

#The command for joining things is:
#By references the shared column name in the two databases

ListOfDrugsForPritzker_vs_ATopGene_DrugBudger_Database1<-join(ListOfDrugsForPritzker, ATopGene_DrugBudger_Database1, by='Drug Name', type="left", match="all")

#You will need to do this for all of your data.frames (18 iterations?), so could be looped or just done by hand (copy, paste, change data.frame name)

#Let's stop here, upload the results to Dropbox, and we'll figure out whether it is worth it to make a master database (vs. just read over the results)

##############################

#In the end, you can make a master database by joining all of our individual databases using a similar method:

dfs<-list(ListOfDrugsForPritzker, ListOfDrugsForPritzker_vs_ATopGene_DrugBudger_Database1, ListOfDrugsForPritzker_vs_ATopGene_DrugBudger_Database2, ...)

MasterDatabase_TopGenes_DrugBudger<-join_all(dfs, by='Drug Name', type = "left", match = "all")

#On second thought, if there are more than one entry per drug per gene per database, we are going to end up with all combinations across databases, which may increase the number of rows exponentially :(

