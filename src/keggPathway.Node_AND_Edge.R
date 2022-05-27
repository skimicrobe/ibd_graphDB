################################################################################
# Script: keggPathway.Node_AND_Edge.R                                          #
# Author: Suyeon Kim                                                           #
# Created: April 15, 2022                                                      #
# Last edited: April 15, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to create a KeggPathway node and sample- #
#       KeggPathway edges file in neo4j bulk importer format                   #
# Input: 1) Three processed keggPathway data (metagenomic, metatranscriptomic, #
#           from ibdmdb database are needed.                                   #
#        2) sample.node.[output argument].csv                                  #
# Output: a keggPathway node and sample_keggPathway edge list                  #
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

# In this function, we want to split "Pathway" column into multiple columns. Input 
# must have a column (called "Pathway")
processing_Dat<- function(all_Pathway){
  split_info<- data.frame(str_split_fixed(all_Pathway$Pathway, "\\:",2))
  colnames(split_info)<- c("PWY","Pathway")
  return(split_info)
}

get_Ids<- function(abundance, entity_Node, entity, sample_info){
  print(dim(abundance))
  # Merge abudance data with Pathway Node dataframe based on entity
  merged_PathwayId<- merge(abundance, entity_Node, by=c(entity))
  filter_cols<- merged_PathwayId %>% select(-c(entity, "PWY"))
  matched_sampleIndex<- 
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampleIndex
  # Create abundance table only with created keggpathway-ids and sample-ids
  final_ab<- cbind("pathway_ID"=filter_cols$path_ID, 
                   filter_cols[,-ncol(filter_cols)])
  print("**********Before Reshape final_ab *******")
  print(dim(final_ab))
  
  return(final_ab)
}

reshape_dat<- function(new_abundance, by_melt){
  reshape_dat<- melt(new_abundance, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID:ID"=as.character(reshape_dat$variable), 
                              "pathway_ID:ID"= reshape_dat$pathway_ID, 
                              "Abundance:double"=reshape_dat$value), 
                        check.names=FALSE)
  print("**********After reshape final_ab and entity_Node dimensions *******")
  print(dim(final_ab))
  
  return(final_ab)
}


################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'mgx_path' that contains 
# 10884 observations of 1639 variables. 
mgx_path<- read.csv(InFile1, header=TRUE, comment.char = "$", sep="\t", 
                    stringsAsFactors = FALSE, check.names=FALSE)

# This should result in a dataframe called 'mtx_path' that contains 
# 6061 observations of 736 variables. 
mtx_path<- read.csv(InFile2, header=TRUE, comment.char = "$", sep="\t", 
                    stringsAsFactors = FALSE, check.names=FALSE)

# Lets read the 'sample node' file into R  
sample_N<- read.csv(sample_node, header=T, sep=",", check.names=FALSE)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))

################################################################################
# Main Program: Create a KEGG Pathway Node                                     #
# Input: mgx_path and mtx_path                                                 #
# Return: Pathway_node                                                         #
################################################################################
source(paste0(sourceDir,"/","filtering_Input.R"))
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))

# The dataframe that has been filtered by 'filt_Input' function. 
filt_mgx<- filt_Input(mgx_path, entity="Pathway")
filt_mtx<- filt_Input(mtx_path, entity="Pathway")  

# First, let's combine each column containing Pathway information from 
# the all inputs. 
all_Pathway<- union(filt_mgx$Pathway, filt_mtx$Pathway)
all_Pathway<- data.frame(cbind("Pathway" = all_Pathway))

# Pre-processed the dataframed called 'all_EC'. 
final_PathNode<- processing_Dat(all_Pathway)

# We want to create and add 'unique_ids' to a EC Node. 
path_uniqIds<- create_Entity_uniqueIds(entity="Pathway node", prefix="path", 
                                       final_propertyTab=final_PathNode)
# Add created above unique-Ids to the EC Node. 
Pathway_node<- cbind("path_ID"=path_uniqIds, final_PathNode)

################################################################################
# Main Program: Create a Sample-Pathway Edge                                   #                         
# Input: Pathway_node, sample_info, filt_mgx, filt_mtx                         #
# Return: path_mgx_reshp, path_mtx_reshp dataframe                             #
################################################################################
# Pre-process a column containing "Pathway" in 'mgx_path' dataframe.   
pInf_mgx<- processing_Dat(filt_mgx)
# Subset mgx samples in sample_info dataframe. 
mgx_sample<- sample_info[sample_info$data_type == "metagenomics",]
# Create a dataframe by using Pathway column and abundance values from 
# metagenomic data 
mgx_ab<- cbind("Pathway"=pInf_mgx$Pathway, filt_mgx %>% select(-"Pathway"))
path_mgx_abundance<- get_Ids(mgx_ab, entity_Node = Pathway_node, 
                             entity="Pathway", sample_info = mgx_sample)
path_mgx_reshp<- reshape_dat(path_mgx_abundance, by_melt="pathway_ID")

# Pre-process a column containing "Pathway" in 'mtx_path' dataframe.
pInf_mtx<- processing_Dat(filt_mtx)
# Subset MTX samples in sample_info dataframe.
mtx_sample<- sample_info[sample_info$data_type == "metatranscriptomics",]
# Create a dataframe by using Pathway column and abundance values from 
# metatranscriptomic data
mtx_ab<- cbind("Pathway"=pInf_mtx$Pathway, filt_mtx %>% select(-"Pathway"))
path_mtx_abundance<- get_Ids(mtx_ab, entity_Node = Pathway_node, 
                             entity="Pathway", sample_info = mtx_sample)
path_mtx_reshp<- reshape_dat(path_mtx_abundance, by_melt="pathway_ID")

################################################################################
# Main Program: Reformat data in neo4j-admin import format and Save Output     #
# Input: Pathway_node, path_mgx_reshp and path_mtx_reshp                       #
# Return: neo4j_Path_Node, neo4j_mgx_Edge, neo4j_mtx_Edge                      #
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
print("*************Formatting dataframe in Neo4j format*********************")
# Process Pathway_node in neo4j-import format 
neo4j_Path_Node<- reformat_Node_file(Pathway_node, old_col="path_ID", 
                                     new_col="pathway_ID:ID", 
                                     node_label="Kegg Pathway")

# Process edge/relationship in neo4j-import format 
neo4j_mgx_Edge<- reformat_Edge_file(path_mgx_reshp, start_col="s_ID:ID", 
                                    end_col="pathway_ID:ID", 
                                    edge_label="HAS_Pathway_ABUNDANCE_IN_MGX")

neo4j_mtx_Edge<- reformat_Edge_file(path_mtx_reshp, start_col="s_ID:ID", 
                                    end_col="pathway_ID:ID", 
                                    edge_label="HAS_Pathway_ABUNDANCE_IN_MTX")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_Path_Node[is.na(neo4j_Path_Node)]<-""
neo4j_mgx_Edge[is.na(neo4j_mgx_Edge)]<- ""
neo4j_mtx_Edge[is.na(neo4j_mtx_Edge)]<- ""

# Write the Pathway node Output file to a .csv file 
fileN= paste("KeggPathway", "node", outName, "csv", sep=".")
write.csv(neo4j_Path_Node, fileN, row.names=F, quote=T, sep="\t")

# Write the sample-Pathway-mgx file to a .csv file 
fileN= paste("sample_KeggPathway_mgx","edge", outName, "csv", sep=".")
write.csv(neo4j_mgx_Edge, fileN, row.names=F, quote=F, sep=",")

# Write the sample-Pathway-mtx file to a .csv file 
fileN= paste("sample_KeggPathway_mtx","edge", outName, "csv", sep=".")
write.csv(neo4j_mtx_Edge, fileN, row.names=F, quote=F, sep=",")
