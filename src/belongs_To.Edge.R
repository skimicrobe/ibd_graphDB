###############################################################################
# Script: BELONG_TO_Edge.R                                                    #
# Author: Suyeon Kim                                                          #
# Created: March 10, 2022                                                     #
# Last edited: March 10, 2022                                                 #
#                                                                             #
# Goal: The purpose of this script is to create edges between 'KeggOrthologue'# 
#       node and 'EC' node as a neo4j bulk importer format                    #  
# Input: KO.node.ibd_dataset.csv, EC.node.ibd_dataset.csv,                    #
#        KEGGOrthology_ECnum.csv                                              #
# Output: belongs_To.edge.ibd_dataset.csv                                     #
###############################################################################

library(stringr)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
InFile3<- args[4]
outName<- args[5]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

###############################################################################
# Read the Input file into R.                                                 #
###############################################################################
# This should result in a dataframe called 'ko_Node'
ko_Node<- read.csv(InFile1, header=TRUE, sep=",", check.names=FALSE)

# This should result in a dataframe called 'ec_Node'
ec_Node<- read.csv(InFile2, header=TRUE, sep=",", check.names=FALSE)

# This should result in a dataframe called 'koEC_info' 
koEc_info<- read.csv(InFile3, header=TRUE, sep=",")

# 
# Pre-processing Data 
#

# We want to subset rows if column 3 (EC_numb) in koEc_info has EC-Ids. 
filt_koEC<- koEc_info[grep("\\EC", koEc_info$EC_number),]
ec_Ids<- data.frame(str_split_fixed(filt_koEC$EC_number, "\\:", 2))
colnames(ec_Ids)<- c("first", "second")

# Let's create a dataframe called 'ko_ec_edge'
ko_ec_edge<- data.frame(cbind("KO"=filt_koEC$KO, "EC"=ec_Ids$second))
ko_ec_edge$EC<- gsub("\\]", "", ko_ec_edge$EC)
 
# Because some KEGG Orthologue entries have multiple EC numbers, they need to be
# split as separate entries and inserted as new rows. 
split_multipleIds<- strsplit(ko_ec_edge$EC, split=" ")
remo_koec<- data.frame(KO=rep(ko_ec_edge$KO, sapply(split_multipleIds,length)), 
                       EC=unlist(split_multipleIds))

###############################################################################
# Main Program: Create a KO-EC edge/relationship                              #
# Input: 'remo_koec' dataframe                                                #         
# Return: ko_ec_Edge dataframe (ko_ID and ec_ID)                              #
###############################################################################
# Merge "KO-EC" Edge and ko_Node" dataframes by 'KO' column. 
matched_KO<- merge(remo_koec,ko_Node,by=c("KO"))
# Merge a created previous dataframe called 'matched_KO' and 
# 'EC_Node' by ' EC' column. 
matched_EC<- merge(matched_KO,ec_Node,by=c("EC"))
# We want to create a 'compound_metabolite' edge. 
ko_ec_Edge<- matched_EC %>% select(c("ko_ID:ID", "ec_ID:ID"))

###############################################################################
# Main Program: Reformat data in neo4j-admin import format and save output    #
# Input: 'ko_ec_Edge' dataframe                                               #                       
# Return: neo4j_ko_ec_Edge                                                    #
# Note: Required source file called "formatting_data_neo4j.R"                 #
###############################################################################
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))

# Process KO-EC edge/relationship in neo4j-import format 
neo4j_ko_ec_Edge<- reformat_Edge_file(ko_ec_Edge, start_col="ko_ID:ID", 
                                    end_col="ec_ID:ID",edge_label="BELONGS_TO")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_ko_ec_Edge[is.na(neo4j_ko_ec_Edge)]<-""

# Write the Compound node Output file to a .csv file 
fileN= paste("belongs_To", "edge", outName, "csv", sep=".")
write.csv(neo4j_ko_ec_Edge, fileN, row.names=F, quote=F)

