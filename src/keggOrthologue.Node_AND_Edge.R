################################################################################
# Script: keggOrthologue.Node_AND_Edge.R                                       #
# Author: Suyeon Kim                                                           #
# Created: April 15, 2022                                                      #
# Last edited: April 15, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to create a KO nodes and sample-KO edges #
#      file in neo4j bulk importer format                                      #
# Input: 1) Two processed KEGG Orthologue data (metagenomic AND meproteomic    #
#           folder) from ibdmdb database are needed.                           #                          
#        2) sample.node.[output argument].csv                                  #
# Output: a KeggOrthologue node and sample_KeggOrthologue edge list            #
################################################################################

library(stringr)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
sample_node<- args[4]
outName<- args[5]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

# In this function, we want to split "KO" column into multiple columns. Input 
# must have a column (called "KO")
processing_Dat<- function(all_ko){
  split_info<- data.frame(str_split_fixed(all_ko$KO, "\\:",2))
  colnames(split_info)<- c("koId","ko_info")
  split_ko_ecId<- data.frame(str_split_fixed(split_info$ko_info, "\\[EC", 2))
  colnames(split_ko_ecId)<- c("KEGGOrthology", "EC_number")
  split_ko_ecId$EC_number<- gsub(":","[EC:", split_ko_ecId$EC_number)
  ko_KEGGinfo_ec<- cbind("KO"=split_info$koId, split_ko_ecId)
  return(ko_KEGGinfo_ec)
}

# In this function, we want to replace the KO name and sample ID with 
# the created sample Ids and KO Ids. 
get_Ids<- function(abundance, entity_Node, entity, sample_info){
  print(dim(abundance))
  # Merge abudance data with KO Node dataframe based on entity
  merged_KOId<- merge(abundance, entity_Node, by=c(entity))
  filter_cols<- merged_KOId %>% select(-c(entity, "KEGGOrthology"))
  matched_sampleIndex<- 
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampleIndex
  # Create abundance table only with created KO-ids and sample-ids
  final_ab<- cbind("ko_ID"=filter_cols$ko_ID, 
                   filter_cols[,-ncol(filter_cols)])
  print("**********Before Reshape final_ab *******")
  print(dim(final_ab))
  
  return(final_ab)
}

# In this function, we convert dataframe from wide to long format.
reshape_dat<- function(new_abundance, by_melt){
  reshape_dat<- melt(new_abundance, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID:ID"=as.character(reshape_dat$variable), 
                              "ko_ID:ID"= reshape_dat$ko_ID, 
                              "Abundance:double"=reshape_dat$value), 
                        check.names=FALSE)
  print("**********After reshape final_ab and entity_Node dimensions *******")
  print(dim(final_ab))
  
  return(final_ab)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
print("Reading Input Data")
# This should result in a dataframe called 'mgx_ko' that contains 
# 7761 observations of 1639 variables. 
mgx_ko<- read.csv(InFile1, header=TRUE, comment.char = "$", sep="\t", 
                  stringsAsFactors = FALSE, check.names=FALSE)
colnames(mgx_ko)[1]<- "KO"
print(dim(mgx_ko))
# This should result in a dataframe called 'mpx_ko' that contains 
# 1823 observations of 451 variables. We want to remove 'UNGROUPED' information
# out in the list. The final dimension of 'mpx_ko' should be 1822 observations 
# of 451 variables. 
mpx_ko<- read.csv(InFile2, header=TRUE, comment.char= "!", check.names=F, 
                  sep="\t", stringsAsFactors = FALSE)
mpx_ko<-mpx_ko[-1,]
print(dim(mpx_ko))

# lets read the 'sample node' file into R  
sample_N<- read.csv(sample_node, header=T, sep=",", check.names=FALSE)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))
print(dim(sample_N))

################################################################################
# Main Program: Create a KEGG Orthologue (KO) Node                             #
# Input: mgx_ko and mpx_ko                                                     #
# Return: KO node dataframe                                                    #
################################################################################
print("Calling formatting_data_neo4j.R")
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))

# First, let's combine each column containing KEGG Orthologue information from 
# both inputs. 
all_ko<- union(mgx_ko$KO, mpx_ko$KO)
all_ko<- data.frame(cbind("KO" = all_ko))

# Pre-processing data 
ko_KEGGinfo_ec<- processing_Dat(all_ko)

# Write the 'KEGG Orthology' with 'EC number' table to a .csv file 
write.csv(ko_KEGGinfo_ec, file="KEGGOrthology_ECnum.csv", row.names=F,quote = T)

# Create a KEGG Orthologue Node. Remove a column ("EC_number") since we don't 
# want this information for creating a KEGG Orthologue Node. 
final_koNode<- ko_KEGGinfo_ec %>% select(-"EC_number")

# We want to create and add 'unique_ids' to a KEGG Orthologue Node. 
ko_uniqIds<- create_Entity_uniqueIds(entity="KEGG Orthologue node", prefix="ko", 
                                     final_propertyTab=final_koNode)
# Add created above unique-Ids to the KEGG Orthologue Node. 
KO_node<- cbind("ko_ID"=ko_uniqIds, final_koNode)
print("We have created KO Node dataframe")

################################################################################
# Main Program: Create a Sample-KEGG Orthologue(KO) edge/relationship          #
# Input: KO_node (KO Node), sample_info, mgx_ko, and mpx_ko                    #
# Return: ko_mgx_reshp and ko_mpx_reshp                                        #
################################################################################
# Pre-process a column containing "KO" in 'mgx_ko' dataframe.   
koInf_mgx<- processing_Dat(mgx_ko)
# Subset MGX samples in sample_info dataframe. 
mgx_sample<- sample_info[sample_info$data_type == "metagenomics",]
# Create a dataframe by using KO_Id column and abundance values from metagenomic
# data 
mgx_ab<- cbind("KO"=koInf_mgx$KO, mgx_ko %>% select(-"KO"))
ko_mgx_abundance<- get_Ids(mgx_ab, entity_Node = KO_node, entity="KO", 
                           sample_info = mgx_sample)
ko_mgx_reshp<- reshape_dat(ko_mgx_abundance, by_melt="ko_ID")

# Pre-process a column coatining "KO" in 'mpx_ko' dataframe.
koInf_mpx<- processing_Dat(mpx_ko)
# Subset MPX samples in sample_info dataframe. 
mpx_sample<- sample_info[sample_info$data_type == "proteomics",]
# Create a dataframe by using KO_Id column and abundance values from  
# metaproteomic data 
mpx_ab<- cbind("KO"= koInf_mpx$KO, mpx_ko %>% select(-"KO"))
ko_mpx_abundance<- get_Ids(mpx_ab, entity_Node = KO_node, entity="KO", 
                           sample_info = mpx_sample)
ko_mpx_reshp<- reshape_dat(ko_mpx_abundance, by_melt="ko_ID")

################################################################################
# Main Program: Reformat data in neo4j-admin import format and Save Output     #
# Input: KO_node, ko_mgx_reshp, and ko_mtx_reshp                               #
# Return: neo4j_KO_Node, neo4j_mgx_Edge, neo4j_mpx_Edge                        #                                     
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
# Process KO_node in neo4j-import format 
neo4j_KO_Node<- reformat_Node_file(KO_node, old_col="ko_ID", 
                                   new_col="ko_ID:ID", 
                                   node_label="KEGG Orthologue")

# Process edge/relationship in neo4j-import format 
neo4j_mgx_Edge<- reformat_Edge_file(ko_mgx_reshp, start_col="s_ID:ID", 
                                    end_col="ko_ID:ID", 
                                    edge_label="HAS_KO_ABUNDANCE_IN_MGX")
neo4j_mpx_Edge<- reformat_Edge_file(ko_mpx_reshp, start_col="s_ID:ID", 
                                    end_col="ko_ID:ID", 
                                    edge_label="HAS_KO_ABUNDANCE_IN_MPX")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_KO_Node[is.na(neo4j_KO_Node)]<-""
neo4j_mgx_Edge[is.na(neo4j_mgx_Edge)]<- ""
neo4j_mpx_Edge[is.na(neo4j_mpx_Edge)]<- ""
print("summary info")
print(dim(neo4j_KO_Node))
print(dim(neo4j_mgx_Edge))
print(dim(neo4j_mpx_Edge))

# Write the KO node Output file to a .csv file 
fileN1= paste("KO", "node", outName, "csv", sep=".")
write.csv(neo4j_KO_Node, fileN1, row.names=F, quote=T, sep="\t")

# Write the sample-KO-mgx file to a .csv file 
fileN2= paste("sample_KO_mgx","edge", outName, "csv", sep=".")
write.csv(neo4j_mgx_Edge, fileN2, row.names=F, quote=F, sep=",")

# Write the sample-KO-mpx file to a .csv file 
fileN3= paste("sample_KO_mpx","edge", outName, "csv", sep=".")
write.csv(neo4j_mpx_Edge, fileN3, row.names=F, quote=F, sep=",")
