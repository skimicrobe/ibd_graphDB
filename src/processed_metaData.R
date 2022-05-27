################################################################################  
# Script: processed_metaData.R                                                 #
# Author: Suyeon Kim                                                           #
# Created: April 14, 2022                                                      #
# Last edited: April 14, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to identify sample Id's information are  #
#    both appearing in 'hmp2_metadata-20180820.csv' and our created metadata.  #
#    Our metadata is created by extracting sample Ids from downloaded all data #
#     files.                                                                   #
#                                                                              #
# Input: 1) hmp2_metadata-20180820.csv (from IBDMDB)                           #
#        2) skim_newmeta.tsv (our created input dataframe)                     #
# Output: correct_metadata.csv file                                            #
################################################################################

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

################################################################################
# Read the Input file into R.                                                  #
################################################################################
# Lets read the metadata (downloaded from ibdmdb database) into R 
ihmp_meta<- read.csv(InFile1, header=TRUE, sep=",", stringsAsFactors =FALSE) 

# Lets read the created by author into R 
cust_meta<- read.csv(InFile2, header=FALSE, sep="\t", stringsAsFactors = FALSE, 
                     check.names=FALSE)
colnames(cust_meta)<- c("External.ID","property","data_source","project_name",
                        "data_type", "date_of_collection")

################################################################################
# Main Program: Merge two metadata to find find common sample-Ids              #
# Input: metadata from IBDMDB and out created metadata                         #
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
