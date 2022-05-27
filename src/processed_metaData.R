################################################################################  
# Script: processed_metaData.R                                                 #
# Author: Suyeon Kim                                                           #
# Created: April 14, 2022                                                      #
# Last edited: April 14, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to 
#
# Input: 1) Three processed Gene data (metatranscriptomic and host             #
#           metranscriptome folder from ibdmdb database are needed.            #
#        2) sample.node.[NAME WILL BE VARIED! - 'test1' is at time of output   #
#           from the script called 'participant_sample.Node_AND_Edge.R'].csv   #
# Output:                         
# (in neo4j bulk importer format).                                             #
################################################################################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#q()
sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
#outName<- args[4]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

################################################################################
# Read the Input file into R.                                                  #
################################################################################
InFile1<- "/Users/suyeonkim/OneDrive - University of Nebraska at Omaha/Research/Research-Fall-2020/Graph-Database-Project/IBD-project/ibd-rawdata/hmp2_metadata-20180820.csv"
InFile2<- "/Users/suyeonkim/OneDrive - University of Nebraska at Omaha/Research/Research-Fall-2020/Graph-Database-Project/IBD-project/ibd-rawdata/skim_newmeta.tsv"

# Lets read the metadata (downloaded from ibdmdb database) into R 
ihmp_meta<- read.csv(InFile1, header=TRUE, sep=",", stringsAsFactors =FALSE) 

# Lets read the created by author into R 
cust_meta<- read.csv(InFile2, header=FALSE, sep="\t", stringsAsFactors = FALSE, 
                     check.names=FALSE)
colnames(cust_meta)<- c("External.ID","property","data_source","project_name",
                        "data_type", "date_of_collection")



################################################################################
# Main Program: Merge two metadata to find find common sample-Ids              #
# Input: ihmp_meta and cust_meta                                               #
#                                                                              #
################################################################################

new_meta<- ihmp_meta[which(ihmp_meta$External.ID %in% cust_meta$External.ID 
                              & ihmp_meta$data_type %in% cust_meta$data_type),]

# Identify unique line based on "ExternalID" and "data_type" 
uniq_cust_meta<-unique(cust_meta[,c("External.ID", "data_source", "data_type")])

# Merge two metadata 
merged_dat<- merge(new_meta, uniq_cust_meta, by=c("External.ID", "data_type"), 
                   all=FALSE)

table(ihmp_meta$data_type)
table(cust_meta$data_type)
table(uniq_cust_meta$data_type)

# Write the processed metadata Output file to a .csv file 
fileN1= paste("correct_metadata", "csv", sep=".")
write.csv(merged_dat, fileN1, row.names=F, quote=T)

#write.csv(merged_dat, file="correct_metadata.csv", row.names=F, quote=T)