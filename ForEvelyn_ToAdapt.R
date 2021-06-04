#Within Excel: 
#Create a new spreadsheet that includes just the columns that you are interested in
#Save as a .csv (comma separated variable)

#Read the .csv into R:
SubjectData<-read.csv("YourFile.csv", header=TRUE, stringsAsFactors=FALSE)

#You can grab particular columns of data using the $
SubjectData$VariableName

#to see the distribution of values in a particular column of categorical data
table(SubjectData$VariableName)

#to see the distribution of two categorical variables in relationship to each other (cross table, similar to a pivot table in Excel):
table(SubjectData$Diagnosis, SubjectData$VariableName)

#if you want to output these tables as a file instead of printing it in the console:
MyTable<-table(SubjectData$Diagnosis, SubjectData$VariableName)
write.csv(MyTable, "FileName_DiagnosisVsMyVariable.csv")

