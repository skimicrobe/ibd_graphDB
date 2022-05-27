################################################################################
# Script: ec.Node_AND_Edge.R                                                   #
# Author: Suyeon Kim                                                           #
# Created: April 16, 2022                                                      #
# Last edited: April 16, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to create a EC nodes and sample-EC edges #
#       file in neo4j bulk importer format                                     #
# Input: 1) Three processed EC data (metagenomic, metatranscriptomic,          #
#           meproteomic folder) from ibdmdb database are needed.               #
#        2) sample.node.[output argument].csv                                  #
# Output: a EC node and sample_EC edge list                                    #  
################################################################################
library(stringr)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
InFile3<- args[4]
sample_node<- args[5]
outName<- args[6]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

# In this function, we want to split "EC" column into multiple columns. Input 
# must have a column (called "EC")
processing_Dat<- function(all_EC){
  split_info<- data.frame(str_split_fixed(all_EC$EC, "\\:",2))
  colnames(split_info)<- c("EC","Enzyme")
  return(split_info)
}

# In this function, we want to replace the EC name and sample ID with 
# the created sample Ids and ec Ids.
get_Ids<- function(abundance, entity_Node, entity, sample_info){
  print(dim(abundance))
  # Merge abudance data with EC Node dataframe based on entity
  merged_ecId<- merge(abundance, entity_Node, by=c(entity))
  filter_cols<- merged_ecId %>% select(-c(entity, "Enzyme"))
  matched_sampleIndex<- 
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampleIndex
  # Create abundance table only with created ec-ids and sample-ids
  final_ab<- cbind("ec_ID"=filter_cols$ec_ID, 
                   filter_cols[,-ncol(filter_cols)])
  print("**********Before Reshape final_ab *******")
  print(dim(final_ab))
  
  return(final_ab)
}

# In this function, we convert dataframe from wide to long format. 
reshape_dat<- function(new_abundance, by_melt){
  reshape_dat<- melt(new_abundance, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID:ID"=as.character(reshape_dat$variable), 
                              "ec_ID:ID"= reshape_dat$ec_ID, 
                              "Abundance:double"=reshape_dat$value), 
                        check.names=FALSE)
  print("**********After reshape final_ab and entity_Node dimensions *******")
  print(dim(final_ab))
  
  return(final_ab)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
source(paste0(sourceDir,"/","filtering_Input.R"))

# This should result in a dataframe called 'mgx_ec' that contains 
# 2173 observations of 1639 variables. 
mgx_ec<- read.csv(InFile1, header=TRUE, comment.char = "$", sep="\t", 
                  stringsAsFactors = FALSE, check.names=FALSE)

# This should result in a dataframe called 'mtx_ec' that contains 
# 70711 observations of 736 variables. 
mtx_ec<- read.csv(InFile2, header=TRUE, comment.char = "$", sep="\t", 
                  stringsAsFactors = FALSE, check.names=FALSE)

# This should result in a dataframe called 'mpx_ec' that contains 
# 910 observations of 451 variables. 
mpx_ec<- read.csv(InFile3, header=TRUE, comment.char = "$", sep="\t", 
                  stringsAsFactors = FALSE, check.names=FALSE)
print("*************Size of ec-dataframe Input***************")
print(dim(mgx_ec))
print(dim(mtx_ec))
print(dim(mpx_ec))

# Lets read the 'sample node' file into R  
sample_N<- read.csv(sample_node, header=T, sep=",",check.names=FALSE)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))

################################################################################
# Main Program: Create a Enzyme Commission (EC) Node                           #
# Input: mgx_ec, mtx_ec, and mpx_ec                                            #
# Return: ec_node dataframe                                                    #
################################################################################
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))

# The dataframe that has been filtered by 'filt_Input' function. 
filt_mgx<- filt_Input(mgx_ec, entity="EC")
filt_mtx<- filt_Input(mtx_ec, entity="EC")  
filt_mpx<- filt_Input(mpx_ec, entity="EC")

# First, let's combine each column containing EC information from 
# the all inputs. 
all_EC<- union(filt_mpx$EC, union(filt_mgx$EC, filt_mtx$EC))
all_EC<- data.frame(cbind("EC" = all_EC))

# Pre-processed the dataframed called 'all_EC'. 
final_ecNode<- processing_Dat(all_EC)

# We want to create and add 'unique_ids' to a EC Node. 
ec_uniqIds<- create_Entity_uniqueIds(entity="Enzyme node", prefix="ec", 
                                     final_propertyTab=final_ecNode)
# Add created above unique-Ids to the EC Node. 
ec_node<- cbind("ec_ID"=ec_uniqIds, final_ecNode)

################################################################################
# Main Program: Create a Sample-Enzyme Commission (EC) Edge                    #
# Input: ec_node, sample_info, filt_mgx, filt_mtx, filt_mpx                    #
# Return: ec_mgx_reshp, ec_mtx_reshp, and ec_mpx_resh                          #
################################################################################
# Pre-process a column containing "EC" in 'mgx_ec' dataframe.   
ecInf_mgx<- processing_Dat(filt_mgx)
# Subset MGX samples in sample_info dataframe. 
mgx_sample<- sample_info[sample_info$data_type == "metagenomics",]
# Create a dataframe by using EC_Id column and abundance values from metagenomic
# data 
mgx_ab<- cbind("EC"=ecInf_mgx$EC, filt_mgx %>% select(-"EC"))
# Replace sample-Ids and ec info with created sample-Ids and ec-Ids. 
ec_mgx_abundance<- get_Ids(mgx_ab, entity_Node = ec_node, entity="EC", 
                           sample_info = mgx_sample)
# Reshape the created dataframe wide to long format 
ec_mgx_reshp<- reshape_dat(ec_mgx_abundance, by_melt="ec_ID")
  
# Pre-process a column containing "EC" in 'mtx_ec' dataframe.
ecInf_mtx<- processing_Dat(filt_mtx)
# Subset MTX samples in sample_info dataframe.
mtx_sample<- sample_info[sample_info$data_type == "metatranscriptomics",]
# Create a dataframe by using EC_Id column and abundance values from 
# metatranscriptomic data
mtx_ab<- cbind("EC"=ecInf_mtx$EC, filt_mtx %>% select(-"EC"))
# Replace sample-Ids and ec info with created sample-Ids and ec-Ids. 
ec_mtx_abundance<- get_Ids(mtx_ab, entity_Node = ec_node, entity="EC", 
                           sample_info = mtx_sample)
# Reshape the created dataframe wide to long format 
ec_mtx_reshp<- reshape_dat(ec_mtx_abundance, by_melt="ec_ID")

# Pre-process a column containing "EC" in 'mpx_ec' dataframe.
ecInf_mpx<- processing_Dat(filt_mpx)
# Subset MPX samples in sample_info dataframe. 
mpx_sample<- sample_info[sample_info$data_type == "proteomics",]
# Create a dataframe by using EC_Id column and abundance values from 
# metaproteomics data
mpx_ab<- cbind("EC"=ecInf_mpx$EC, filt_mpx %>% select(-"EC"))
# Replace sample-Ids and ec info with created sample-Ids and ec-Ids. 
ec_mpx_abundance<- get_Ids(mpx_ab, entity_Node = ec_node, entity="EC", 
                           sample_info = mpx_sample)
# Reshape the created dataframe wide to long format 
ec_mpx_reshp<- reshape_dat(ec_mpx_abundance, by_melt="ec_ID")

################################################################################
# Main Program: Reformat data in neo4j-admin import format and Save Output     #
# Input: EC_node, ec_mgx_reshp, ec_mtx_reshp, and ec_mpx_resh                  #            
# Return: neo4j_EC_node, noe4j_mgx/mtx/mpx_Edge                                #
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
print("*************Formatting dataframe in Neo4j format*********************")
# Process EC_node in neo4j-import format 
neo4j_EC_Node<- reformat_Node_file(ec_node, old_col="ec_ID", 
                                   new_col="ec_ID:ID", 
                                   node_label="Enzyme")

# Process edge/relationship in neo4j-import format 
neo4j_mgx_Edge<- reformat_Edge_file(ec_mgx_reshp, start_col="s_ID:ID", 
                                    end_col="ec_ID:ID", 
                                    edge_label="HAS_EC_ABUNDANCE_IN_MGX")
neo4j_mtx_Edge<- reformat_Edge_file(ec_mtx_reshp, start_col="s_ID:ID", 
                                    end_col="ec_ID:ID", 
                                    edge_label="HAS_EC_ABUNDANCE_IN_MTX")
neo4j_mpx_Edge<- reformat_Edge_file(ec_mpx_reshp, start_col="s_ID:ID", 
                                    end_col="ec_ID:ID", 
                                    edge_label="HAS_KO_ABUNDANCE_IN_MPX")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_EC_Node[is.na(neo4j_EC_Node)]<-""
neo4j_mgx_Edge[is.na(neo4j_mgx_Edge)]<- ""
neo4j_mtx_Edge[is.na(neo4j_mtx_Edge)]<- ""
neo4j_mpx_Edge[is.na(neo4j_mpx_Edge)]<- ""

# Write the EO node Output file to a .csv file 
fileN= paste("EC", "node", outName, "csv", sep=".")
write.csv(neo4j_EC_Node, fileN, row.names=F, quote=T, sep="\t")

# Write the sample-EC-mgx file to a .csv file 
fileN= paste("sample_EC_mgx","edge", outName, "csv", sep=".")
write.csv(neo4j_mgx_Edge, fileN, row.names=F, quote=F, sep=",")

# Write the sample-EC-mtx file to a .csv file 
fileN= paste("sample_EC_mtx","edge", outName, "csv", sep=".")
write.csv(neo4j_mtx_Edge, fileN, row.names=F, quote=F, sep=",")

# Write the sample-EC-mpx file to a .csv file 
fileN= paste("sample_EC_mpx","edge", outName, "csv", sep=".")
write.csv(neo4j_mpx_Edge, fileN, row.names=F, quote=F, sep=",")
